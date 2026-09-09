//  (C) Copyright Matt Borland 2026.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
//  CUDA backend of the prime sieve. Only compiled by nvcc with BOOST_MATH_ENABLE_CUDA.
//
//  The range is processed in chunks of up to 2^32 integers. Each chunk is sieved by a grid of
//  thread blocks, one block per 32 KiB shared-memory segment (983040 integers) in the same
//  mod-30 byte layout as the CPU engine: primes 7..37 are removed by ANDing periodic tables,
//  primes 41..1024 are handled one warp per prime, and primes up to the segment span one
//  thread per prime with atomicAnd on shared memory. Primes above the segment span (at most
//  one multiple per segment) are applied in a second kernel directly on the chunk's bitset in
//  global memory. Counting uses __popc with a block reduction; generation compacts the set
//  bits into a sorted 64-bit array through a two-level exclusive scan.

#ifndef BOOST_MATH_SF_DETAIL_PRIME_SIEVE_CUDA_HPP
#define BOOST_MATH_SF_DETAIL_PRIME_SIEVE_CUDA_HPP

#include <boost/math/tools/config.hpp>

#if defined(__CUDACC__) && defined(BOOST_MATH_ENABLE_CUDA) && !defined(BOOST_MATH_HAS_NVRTC)

#define BOOST_MATH_HAS_CUDA_PRIME_SIEVE

#include <boost/math/special_functions/detail/prime_sieve/layout.hpp>
#include <boost/math/special_functions/detail/prime_sieve/options.hpp>
#include <boost/math/special_functions/detail/prime_sieve/presieve.hpp>
#include <boost/math/special_functions/detail/prime_sieve/driver.hpp>
#include <boost/math/special_functions/detail/prime_sieve/parallel.hpp>
#include <boost/math/special_functions/detail/prime_sieve/primality.hpp>
#include <cuda_runtime.h>
#include <cstdint>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <string>
#include <stdexcept>
#include <algorithm>

namespace boost::math::detail::prime_sieve {

namespace cuda_detail {

inline constexpr unsigned segment_bytes {32768};                       // shared memory per block
inline constexpr unsigned segment_span {segment_bytes * 30u};          // integers per segment
inline constexpr unsigned block_threads {256};
inline constexpr unsigned warp_prime_limit {1024};                      // primes up to this: one warp per prime
inline constexpr unsigned presieve_limit {37};                          // primes up to this: periodic tables
inline constexpr unsigned segments_per_chunk {4369};                     // about 2^32 integers per device chunk
inline constexpr std::size_t chunk_bytes {static_cast<std::size_t>(segments_per_chunk) * segment_bytes};
inline constexpr std::uint64_t chunk_numbers {static_cast<std::uint64_t>(chunk_bytes) * 30u};   // multiple of 30

// Reports a failed CUDA call. Without host exceptions there is nowhere to return the
// error to, so the message goes to stderr and the process aborts.
inline void cuda_check(cudaError_t e, const char* what)
{
    if (e != cudaSuccess)
    {
#ifndef BOOST_MATH_PRIME_SIEVE_NO_EXCEPTIONS
        throw std::runtime_error(std::string("boost::math::prime_sieve CUDA error in ") + what + ": " + cudaGetErrorString(e));
#else
        std::fprintf(stderr, "boost::math::prime_sieve CUDA error in %s: %s\n", what, cudaGetErrorString(e));
        std::abort();
#endif
    }
}

template <class T>
class device_buffer
{
public:
    device_buffer() = default;

    explicit device_buffer(std::size_t n)
    {
        allocate(n);
    }

    ~device_buffer()
    {
        release();
    }

    device_buffer(const device_buffer&) = delete;
    device_buffer& operator=(const device_buffer&) = delete;

    void allocate(std::size_t n)
    {
        release();
        if (n != 0)
        {
            void* p {nullptr};
            cuda_check(cudaMalloc(&p, n * sizeof(T)), "cudaMalloc");
            ptr_ = static_cast<T*>(p);
            size_ = n;
        }
    }

    void release() noexcept
    {
        if (ptr_ != nullptr)
        {
            cudaFree(ptr_);
            ptr_ = nullptr;
            size_ = 0;
        }
    }

    void upload(const T* host, std::size_t n)
    {
        cuda_check(cudaMemcpy(ptr_, host, n * sizeof(T), cudaMemcpyHostToDevice), "cudaMemcpy H2D");
    }

    void download(T* host, std::size_t n) const
    {
        cuda_check(cudaMemcpy(host, ptr_, n * sizeof(T), cudaMemcpyDeviceToHost), "cudaMemcpy D2H");
    }

    T* get() noexcept
    {
        return ptr_;
    }

    const T* get() const noexcept
    {
        return ptr_;
    }

    std::size_t size() const noexcept
    {
        return size_;
    }

private:
    T* ptr_ {nullptr};
    std::size_t size_ {0};
};

// Page-locked host memory for fast device to host copies; falls back to pageable memory.
class pinned_buffer
{
public:
    pinned_buffer() = default;

    ~pinned_buffer()
    {
        release();
    }

    pinned_buffer(const pinned_buffer&) = delete;
    pinned_buffer& operator=(const pinned_buffer&) = delete;

    void reserve(std::size_t n)
    {
        if (n <= size_)
        {
            return;
        }
        release();
        void* p {nullptr};
        if (cudaMallocHost(&p, n * sizeof(std::uint64_t)) == cudaSuccess)
        {
            ptr_ = static_cast<std::uint64_t*>(p);
            pinned_ = true;
        }
        else
        {
            (void)cudaGetLastError();
            fallback_.resize(n);
            ptr_ = fallback_.data();
            pinned_ = false;
        }
        size_ = n;
    }

    std::uint64_t* data() noexcept
    {
        return ptr_;
    }

private:
    void release() noexcept
    {
        if (pinned_ && ptr_ != nullptr)
        {
            cudaFreeHost(ptr_);
        }
        fallback_.clear();
        fallback_.shrink_to_fit();
        ptr_ = nullptr;
        size_ = 0;
        pinned_ = false;
    }

    std::uint64_t* ptr_ {nullptr};
    std::size_t size_ {0};
    bool pinned_ {false};
    std::vector<std::uint64_t> fallback_;
};

// Tables the kernels need, passed by value so no device globals are required. Kernels copy
// them into shared memory first: constant memory serializes lane-divergent indexing.
struct alignas(4) device_tables
{
    std::uint8_t unset_bit[8];
    std::uint8_t bit_values[64];
    std::uint8_t keep_low[37];
    std::uint8_t keep_high[37];
    std::uint8_t class_of_residue[30];   // bit_of_residue for residues coprime to 30
    std::uint8_t inverse_mod_30[30];     // multiplicative inverse for residues coprime to 30
    std::uint8_t init_factor[30];        // wheel30 init: distance to the next coprime multiplier
    std::uint8_t init_phase[30];         // wheel30 init: phase of that multiplier
    std::uint8_t padding[2];
    std::uint32_t wheel[64];             // wheel30 element by state: unset | factor << 8 | correct << 16 | next << 24
};

__device__ __forceinline__ std::uint8_t wheel_unset(std::uint32_t w) { return static_cast<std::uint8_t>(w); }
__device__ __forceinline__ std::uint32_t wheel_factor(std::uint32_t w) { return (w >> 8) & 0xffu; }
__device__ __forceinline__ std::uint32_t wheel_correct(std::uint32_t w) { return (w >> 16) & 0xffu; }
__device__ __forceinline__ unsigned wheel_next(std::uint32_t w) { return w >> 24; }

static_assert(sizeof(device_tables) % 4 == 0, "device_tables must be copyable in 32-bit words");

// Block-cooperative copy of the kernel's table argument into shared memory.
__device__ __forceinline__ void load_tables(const device_tables& from, device_tables& to)
{
    const std::uint32_t* src {reinterpret_cast<const std::uint32_t*>(&from)};
    std::uint32_t* dst {reinterpret_cast<std::uint32_t*>(&to)};
    for (unsigned i {threadIdx.x}; i < sizeof(device_tables) / 4; i += blockDim.x)
    {
        dst[i] = src[i];
    }
    __syncthreads();
}

constexpr device_tables make_device_tables() noexcept
{
    device_tables t {};
    for (unsigned k {0}; k < 8; ++k)
    {
        t.unset_bit[k] = unset_bit[k];
    }
    for (unsigned b {0}; b < 64; ++b)
    {
        t.bit_values[b] = bit_values[b];
        t.wheel[b] = static_cast<std::uint32_t>(wheel30.wheel[b].unset_bit)
                   | (static_cast<std::uint32_t>(wheel30.wheel[b].next_multiple_factor) << 8)
                   | (static_cast<std::uint32_t>(wheel30.wheel[b].correct) << 16)
                   | (static_cast<std::uint32_t>(wheel30.wheel[b].next) << 24);
    }
    for (unsigned d {0}; d < 37; ++d)
    {
        t.keep_low[d] = keep_low[d];
        t.keep_high[d] = keep_high[d];
    }
    for (unsigned r {0}; r < 30; ++r)
    {
        t.class_of_residue[r] = bit_of_residue(r);
        t.init_factor[r] = wheel30.init[r].next_multiple_factor;
        t.init_phase[r] = wheel30.init[r].wheel_index;
        t.inverse_mod_30[r] = 0;
        for (unsigned x {1}; x < 30; ++x)
        {
            if ((r * x) % 30 == 1)
            {
                t.inverse_mod_30[r] = static_cast<std::uint8_t>(x);
            }
        }
    }
    return t;
}

// Built at compile time; copied into a kernel argument per launch.
inline constexpr device_tables device_tables_value = make_device_tables();

// Per sieving prime: p, chunk_low mod p (advanced per chunk) and the Barrett constant for p,
// so that segment_low mod p needs no 64-bit division.
struct device_prime
{
    std::uint32_t p;
    std::uint32_t residue;     // chunk_low mod p
    std::uint64_t barrett;     // floor(2^64 / p)
};

// r = t mod p for t < 2^63 using the precomputed floor(2^64 / p).
__device__ __forceinline__ std::uint32_t barrett_mod(std::uint64_t t, std::uint32_t p, std::uint64_t m)
{
    std::uint64_t q {__umul64hi(t, m)};
    std::uint64_t r {t - q * p};
    if (r >= p)
    {
        r -= p;
    }
    return static_cast<std::uint32_t>(r);
}

// Integer offset from segment_low and wheel state of the first multiple of p that is at least
// max(segment_low, p * p) and coprime to 30. segment_low_mod_p is segment_low modulo p.
__device__ __forceinline__ void first_multiple_offset(const device_tables& t, std::uint32_t p, std::uint64_t segment_low,
                                                      std::uint32_t segment_low_mod_p, std::uint64_t& d, unsigned& state)
{
    // d0: offset (in integers) of the first multiple >= segment_low
    d = segment_low_mod_p == 0 ? 0u : p - segment_low_mod_p;
    const unsigned p_mod_30 {p % 30u};
    unsigned q_mod_30 {};
    const std::uint64_t square {static_cast<std::uint64_t>(p) * p};
    if (segment_low + d < square)
    {
        // never cross off the prime itself or multiples a smaller prime already removed
        d = square - segment_low;
        q_mod_30 = p_mod_30;
    }
    else
    {
        // q0 mod 30 from d0 mod 30 (segment_low is a multiple of 30) and the inverse of p mod 30
        q_mod_30 = static_cast<unsigned>((d % 30u) * t.inverse_mod_30[p_mod_30] % 30u);
    }
    d += static_cast<std::uint64_t>(t.init_factor[q_mod_30]) * p;
    state = t.class_of_residue[p_mod_30] * 8u + t.init_phase[q_mod_30];
    // offsets below 7 belong to the previous segment's last byte: step once more
    while (d < 7)
    {
        const std::uint32_t w {t.wheel[state]};
        d += static_cast<std::uint64_t>(wheel_factor(w)) * p;
        state = wheel_next(w);
    }
}

// Byte index (relative to the segment) and wheel state of the first multiple in a segment.
// Multiples beyond the segment yield an index >= segment_bytes.
__device__ __forceinline__ void first_multiple_device(const device_tables& t, std::uint32_t p, std::uint64_t segment_low,
                                                      std::uint32_t segment_low_mod_p, std::uint32_t& byte_index, std::uint32_t& state)
{
    std::uint64_t d {};
    unsigned s {};
    first_multiple_offset(t, p, segment_low, segment_low_mod_p, d, s);
    const std::uint64_t byte {(d - 7) / 30};
    byte_index = byte > segment_bytes ? segment_bytes : static_cast<std::uint32_t>(byte);
    state = s;
}

__device__ __forceinline__ void clear_bit(std::uint32_t* words, std::uint32_t byte_index, std::uint8_t mask)
{
    const std::uint32_t word_mask {~(static_cast<std::uint32_t>(static_cast<std::uint8_t>(~mask)) << (8u * (byte_index & 3u)))};
    atomicAnd(&words[byte_index >> 2], word_mask);
}

// One block sieves one segment of the chunk. primes[0, n_warp) are handled per warp, primes
// [n_warp, n_thread) per thread. The finished segment is written to chunk_words.
__global__ void __launch_bounds__(block_threads)
sieve_segments_kernel(device_tables table_arg, const device_prime* primes, std::uint32_t n_warp, std::uint32_t n_thread,
                      const std::uint8_t* presieve0, std::uint32_t period0, const std::uint8_t* presieve1, std::uint32_t period1,
                      const std::uint8_t* presieve2, std::uint32_t period2,
                      std::uint32_t* chunk_words, std::uint64_t chunk_low, std::uint64_t start, std::uint64_t stop, std::uint32_t n_segments)
{
    __shared__ std::uint32_t words[segment_bytes / 4];
    __shared__ device_tables tables;
    const std::uint32_t segment {blockIdx.x};
    if (segment >= n_segments)
    {
        return;
    }
    load_tables(table_arg, tables);
    const std::uint64_t segment_low {chunk_low + static_cast<std::uint64_t>(segment) * segment_span};
    std::uint8_t* bytes {reinterpret_cast<std::uint8_t*>(words)};

    // bytes actually covering [segment_low + 7, stop]
    std::uint32_t used_bytes {segment_bytes};
    if (stop - segment_low < static_cast<std::uint64_t>(segment_span) + 1)
    {
        const std::uint64_t remaining {stop - segment_low};
        used_bytes = remaining < 7 ? 0u : static_cast<std::uint32_t>((remaining - 7) / 30 + 1);
    }

    // phase 0: pre-sieve 7..37 from the periodic tables
    {
        const std::uint32_t byte_low {static_cast<std::uint32_t>((segment_low / 30) % period0)};
        const std::uint32_t byte_low1 {static_cast<std::uint32_t>((segment_low / 30) % period1)};
        const std::uint32_t byte_low2 {static_cast<std::uint32_t>((segment_low / 30) % period2)};
        for (std::uint32_t i {threadIdx.x}; i < segment_bytes; i += block_threads)
        {
            std::uint32_t p0 {byte_low + i};
            if (p0 >= period0) { p0 -= period0; if (p0 >= period0) { p0 %= period0; } }
            std::uint32_t p1 {byte_low1 + i};
            if (p1 >= period1) { p1 -= period1; if (p1 >= period1) { p1 %= period1; } }
            std::uint32_t p2 {byte_low2 + i};
            if (p2 >= period2) { p2 -= period2; if (p2 >= period2) { p2 %= period2; } }
            bytes[i] = presieve0[p0] & presieve1[p1] & presieve2[p2];
        }
    }
    // the first segments of the number line contain the pre-sieved primes themselves
    if (threadIdx.x == 0 && segment_low <= presieve_limit)
    {
        const std::uint8_t small[9] = {7, 11, 13, 17, 19, 23, 29, 31, 37};
        for (unsigned k {0}; k < 9; ++k)
        {
            if (small[k] >= segment_low + 7)
            {
                bytes[(small[k] - segment_low - 7) / 30] |= static_cast<std::uint8_t>(1u << tables.class_of_residue[small[k] % 30]);
            }
        }
    }
    __syncthreads();

    // phase 1: warp per prime for the small primes
    {
        const unsigned lane {threadIdx.x & 31u};
        const unsigned warp {threadIdx.x >> 5};
        const unsigned warps {block_threads / 32u};
        for (std::uint32_t idx {warp}; idx < n_warp; idx += warps)
        {
            const device_prime pr {primes[idx]};
            const std::uint32_t span_mod {barrett_mod(segment_span, pr.p, pr.barrett)};
            const std::uint32_t low_mod {barrett_mod(static_cast<std::uint64_t>(pr.residue) + static_cast<std::uint64_t>(segment) * span_mod, pr.p, pr.barrett)};
            std::uint32_t i {};
            std::uint32_t state {};
            first_multiple_device(tables, pr.p, segment_low, low_mod, i, state);
            // one wheel turn is exactly p bytes; the eight offsets within a turn are fixed
            std::uint32_t offsets[8];
            std::uint8_t masks[8];
            const std::uint32_t p30 {pr.p / 30u};
            std::uint32_t o {0};
            for (unsigned k {0}; k < 8; ++k)
            {
                const std::uint32_t w {tables.wheel[state]};
                offsets[k] = o;
                masks[k] = wheel_unset(w);
                o += p30 * wheel_factor(w) + wheel_correct(w);
                state = wheel_next(w);
            }
            for (std::uint32_t turn_start {i + lane * pr.p}; turn_start < segment_bytes; turn_start += 32u * pr.p)
            {
                for (unsigned k {0}; k < 8; ++k)
                {
                    const std::uint32_t b {turn_start + offsets[k]};
                    if (b < segment_bytes)
                    {
                        clear_bit(words, b, masks[k]);
                    }
                }
            }
        }
    }

    // phase 2: thread per prime for the medium primes
    for (std::uint32_t idx {n_warp + threadIdx.x}; idx < n_thread; idx += block_threads)
    {
        const device_prime pr {primes[idx]};
        const std::uint32_t span_mod {barrett_mod(segment_span, pr.p, pr.barrett)};
        const std::uint32_t low_mod {barrett_mod(static_cast<std::uint64_t>(pr.residue) + static_cast<std::uint64_t>(segment) * span_mod, pr.p, pr.barrett)};
        std::uint32_t i {};
        std::uint32_t state {};
        first_multiple_device(tables, pr.p, segment_low, low_mod, i, state);
        const std::uint32_t p30 {pr.p / 30u};
        while (i < segment_bytes)
        {
            const std::uint32_t w {tables.wheel[state]};
            clear_bit(words, i, wheel_unset(w));
            i += p30 * wheel_factor(w) + wheel_correct(w);
            state = wheel_next(w);
        }
    }
    __syncthreads();

    // phase 3: range masks, tail zeroing, store to global memory
    if (threadIdx.x == 0)
    {
        if (segment_low < start && start - segment_low < 37)
        {
            bytes[0] &= tables.keep_low[start - segment_low];
        }
        if (used_bytes != 0 && used_bytes < segment_bytes)
        {
            const std::uint64_t e {stop - segment_low - 30u * static_cast<std::uint64_t>(used_bytes - 1)};
            bytes[used_bytes - 1] &= tables.keep_high[e > 36 ? 36 : static_cast<unsigned>(e)];
        }
        else if (used_bytes == segment_bytes && stop - segment_low < static_cast<std::uint64_t>(segment_span) + 1)
        {
            const std::uint64_t e {stop - segment_low - 30u * static_cast<std::uint64_t>(segment_bytes - 1)};
            bytes[segment_bytes - 1] &= tables.keep_high[e > 36 ? 36 : static_cast<unsigned>(e)];
        }
    }
    __syncthreads();
    for (std::uint32_t i {threadIdx.x}; i < segment_bytes; i += block_threads)
    {
        if (i >= used_bytes)
        {
            bytes[i] = 0;
        }
    }
    __syncthreads();
    std::uint32_t* out {chunk_words + static_cast<std::size_t>(segment) * (segment_bytes / 4)};
    for (std::uint32_t w {threadIdx.x}; w < segment_bytes / 4; w += block_threads)
    {
        out[w] = words[w];
    }
}

// Primes above the segment span: one thread per prime over the whole chunk bitset.
__global__ void __launch_bounds__(block_threads)
sieve_big_kernel(device_tables table_arg, const device_prime* primes, std::uint32_t first, std::uint32_t count,
                 std::uint32_t* chunk_words, std::uint64_t chunk_low, std::uint64_t stop, std::uint64_t chunk_bytes_used)
{
    __shared__ device_tables tables;
    load_tables(table_arg, tables);
    const std::uint32_t idx {first + blockIdx.x * block_threads + threadIdx.x};
    if (idx >= count)
    {
        return;
    }
    const device_prime pr {primes[idx]};
    std::uint64_t d {};
    unsigned state {};
    first_multiple_offset(tables, pr.p, chunk_low, pr.residue, d, state);
    const std::uint64_t limit {stop - chunk_low};
    while (d <= limit)
    {
        const std::uint64_t byte {(d - 7) / 30};
        if (byte >= chunk_bytes_used)
        {
            break;
        }
        const std::uint32_t w {tables.wheel[state]};
        const std::uint32_t word_mask {~(static_cast<std::uint32_t>(static_cast<std::uint8_t>(~wheel_unset(w))) << (8u * static_cast<unsigned>(byte & 3u)))};
        atomicAnd(&chunk_words[byte >> 2], word_mask);
        d += static_cast<std::uint64_t>(wheel_factor(w)) * pr.p;
        state = wheel_next(w);
    }
}

// residue = (residue + chunk_mod) mod p after a chunk advanced by chunk_numbers.
__global__ void advance_residues_kernel(device_prime* primes, std::uint32_t count)
{
    const std::uint32_t idx {blockIdx.x * block_threads + threadIdx.x};
    if (idx < count)
    {
        const std::uint32_t p {primes[idx].p};
        std::uint32_t r {primes[idx].residue + static_cast<std::uint32_t>(chunk_numbers % p)};
        if (r >= p)
        {
            r -= p;
        }
        primes[idx].residue = r;
    }
}

// Per-block popcounts over 256 * 4 words each (block_totals) and a global total.
__global__ void __launch_bounds__(block_threads)
count_kernel(const std::uint32_t* words, std::uint32_t n_words, std::uint32_t* block_totals, unsigned long long* total)
{
    __shared__ std::uint32_t partial[block_threads / 32];
    std::uint32_t c {0};
    const std::uint32_t base {blockIdx.x * block_threads * 4u};
    for (unsigned k {0}; k < 4; ++k)
    {
        const std::uint32_t w {base + k * block_threads + threadIdx.x};
        if (w < n_words)
        {
            c += static_cast<std::uint32_t>(__popc(words[w]));
        }
    }
    for (unsigned offset {16}; offset > 0; offset >>= 1)
    {
        c += __shfl_down_sync(0xffffffffu, c, offset);
    }
    if ((threadIdx.x & 31u) == 0)
    {
        partial[threadIdx.x >> 5] = c;
    }
    __syncthreads();
    if (threadIdx.x == 0)
    {
        std::uint32_t sum {0};
        for (unsigned i {0}; i < block_threads / 32; ++i)
        {
            sum += partial[i];
        }
        if (block_totals != nullptr)
        {
            block_totals[blockIdx.x] = sum;
        }
        if (total != nullptr)
        {
            atomicAdd(total, static_cast<unsigned long long>(sum));
        }
    }
}

// Exclusive scan of the block totals in place (single block), writing the grand total.
__global__ void scan_totals_kernel(std::uint32_t* totals, std::uint32_t n, unsigned long long* grand_total)
{
    __shared__ unsigned long long carry;
    __shared__ std::uint32_t values[1024];
    if (threadIdx.x == 0)
    {
        carry = 0;
    }
    __syncthreads();
    for (std::uint32_t base {0}; base < n; base += 1024)
    {
        const std::uint32_t i {base + threadIdx.x};
        const std::uint32_t v {i < n ? totals[i] : 0u};
        values[threadIdx.x] = v;
        __syncthreads();
        // Hillis-Steele inclusive scan in shared memory
        for (unsigned offset {1}; offset < 1024; offset <<= 1)
        {
            const std::uint32_t add {threadIdx.x >= offset ? values[threadIdx.x - offset] : 0u};
            __syncthreads();
            values[threadIdx.x] += add;
            __syncthreads();
        }
        const unsigned long long exclusive {carry + values[threadIdx.x] - v};
        if (i < n)
        {
            totals[i] = static_cast<std::uint32_t>(exclusive);
        }
        __syncthreads();
        if (threadIdx.x == 1023)
        {
            carry += values[1023];
        }
        __syncthreads();
    }
    if (threadIdx.x == 0 && grand_total != nullptr)
    {
        *grand_total = carry;
    }
}

// Writes the integers of the set bits, in order, using the scanned block offsets.
__global__ void __launch_bounds__(block_threads)
compact_kernel(device_tables table_arg, const std::uint32_t* words, std::uint32_t n_words, const std::uint32_t* block_offsets,
               std::uint64_t chunk_low, std::uint64_t* out)
{
    __shared__ std::uint32_t thread_counts[block_threads];
    __shared__ device_tables tables;
    load_tables(table_arg, tables);
    const std::uint32_t base {blockIdx.x * block_threads * 4u};
    std::uint32_t c {0};
    // each thread owns four consecutive words so that the output stays sorted
    for (unsigned k {0}; k < 4; ++k)
    {
        const std::uint32_t w {base + 4u * threadIdx.x + k};
        if (w < n_words)
        {
            c += static_cast<std::uint32_t>(__popc(words[w]));
        }
    }
    thread_counts[threadIdx.x] = c;
    __syncthreads();
    for (unsigned offset {1}; offset < block_threads; offset <<= 1)
    {
        const std::uint32_t add {threadIdx.x >= offset ? thread_counts[threadIdx.x - offset] : 0u};
        __syncthreads();
        thread_counts[threadIdx.x] += add;
        __syncthreads();
    }
    std::uint64_t pos {block_offsets[blockIdx.x] + (thread_counts[threadIdx.x] - c)};
    for (unsigned k {0}; k < 4; ++k)
    {
        const std::uint32_t w {base + 4u * threadIdx.x + k};
        if (w < n_words)
        {
            std::uint32_t bits {words[w]};
            const std::uint64_t low {chunk_low + static_cast<std::uint64_t>(w) * 120u};
            while (bits != 0)
            {
                const int b {__ffs(bits) - 1};
                out[pos++] = low + tables.bit_values[b];
                bits &= bits - 1;
            }
        }
    }
}

// Host-side state for one call.
class cuda_sieve
{
public:
    cuda_sieve(std::uint64_t start, std::uint64_t stop) : start_(start), stop_(stop), tables_(device_tables_value)
    {
        // the chunk grid begins at a multiple of 30 covering start
        chunk_low_ = 30u * ((start - 7) / 30);

        // sieving primes: 41..163 from the small table, 167.. from the CPU engine
        const std::uint64_t root {boost::math::tools::isqrt(stop)};
        std::vector<std::uint32_t> primes;
        for (std::size_t i {0}; i < 71 && small_primes_to_353[i] <= presieve_max_prime; ++i)
        {
            if (small_primes_to_353[i] > presieve_limit && small_primes_to_353[i] <= root)
            {
                primes.push_back(small_primes_to_353[i]);
            }
        }
#ifdef BOOST_MATH_PRIME_SIEVE_HAS_THREADS
        const std::vector<std::uint32_t> rest {sieving_primes_upto_parallel(root, prime_sieve_options {})};
#else
        const std::vector<std::uint32_t> rest {sieving_primes_upto(root)};
#endif
        primes.insert(primes.end(), rest.begin(), rest.end());

        std::vector<device_prime> records(primes.size());
        n_warp_ = 0;
        n_thread_ = 0;
        for (std::size_t i {0}; i < primes.size(); ++i)
        {
            const std::uint32_t p {primes[i]};
            records[i].p = p;
            records[i].residue = static_cast<std::uint32_t>(chunk_low_ % p);
            records[i].barrett = (std::numeric_limits<std::uint64_t>::max)() / p;
            if (p <= warp_prime_limit)
            {
                n_warp_ = static_cast<std::uint32_t>(i + 1);
            }
            if (p <= segment_span)
            {
                n_thread_ = static_cast<std::uint32_t>(i + 1);
            }
        }
        n_primes_ = static_cast<std::uint32_t>(records.size());
        if (!records.empty())
        {
            primes_.allocate(records.size());
            primes_.upload(records.data(), records.size());
        }

        // pre-sieve tables for 7..37
        const std::uint8_t groups[3][3] = {{7, 11, 13}, {17, 19, 23}, {29, 31, 37}};
        for (unsigned g {0}; g < 3; ++g)
        {
            std::size_t period {1};
            for (const std::uint8_t p : groups[g])
            {
                period *= p;
            }
            std::vector<std::uint8_t> pattern(period, 0xff);
            for (const std::uint8_t p : groups[g])
            {
                presieve_mark(pattern.data(), period, p);
            }
            presieve_period_[g] = static_cast<std::uint32_t>(period);
            presieve_[g].allocate(period);
            presieve_[g].upload(pattern.data(), period);
        }

        words_.allocate(chunk_bytes / 4);
        const std::size_t count_blocks {(chunk_bytes / 4 + block_threads * 4 - 1) / (block_threads * 4)};
        block_totals_.allocate(count_blocks);
        totals_.allocate(1);
    }

    // Runs f(chunk_low, chunk_stop) for every chunk with the chunk's bitset sieved into words_.
    template <class F>
    void for_each_chunk(F&& f)
    {
        std::uint64_t chunk_low {chunk_low_};
        while (true)
        {
            const std::uint64_t remaining {stop_ - chunk_low};
            if (remaining < 7)
            {
                break;
            }
            const std::uint64_t needed_bytes {(remaining - 7) / 30 + 1};
            const bool last {needed_bytes <= chunk_bytes};
            const std::uint64_t used_bytes {last ? needed_bytes : chunk_bytes};
            // A full chunk covers every candidate of its bytes, including offset 31 of the last
            // byte (chunk_low + chunk_numbers + 1), which the next chunk cannot represent.
            const std::uint64_t chunk_stop {last ? stop_ : chunk_low + chunk_numbers + 1};
            const std::uint32_t n_segments {static_cast<std::uint32_t>((used_bytes + segment_bytes - 1) / segment_bytes)};
            const std::size_t n_words {static_cast<std::size_t>(n_segments) * (segment_bytes / 4)};

            sieve_segments_kernel<<<n_segments, block_threads>>>(tables_, primes_.get(), n_warp_, n_thread_,
                presieve_[0].get(), presieve_period_[0], presieve_[1].get(), presieve_period_[1], presieve_[2].get(), presieve_period_[2],
                words_.get(), chunk_low, start_, chunk_stop, n_segments);
            cuda_check(cudaGetLastError(), "sieve_segments_kernel launch");
            if (n_primes_ > n_thread_)
            {
                const std::uint32_t big {n_primes_ - n_thread_};
                sieve_big_kernel<<<(big + block_threads - 1) / block_threads, block_threads>>>(tables_, primes_.get(), n_thread_, n_primes_,
                    words_.get(), chunk_low, chunk_stop, used_bytes);
                cuda_check(cudaGetLastError(), "sieve_big_kernel launch");
            }
            f(chunk_low, chunk_stop, static_cast<std::uint32_t>(n_words));
            if (last)
            {
                break;
            }
            advance_residues_kernel<<<(n_primes_ + block_threads - 1) / block_threads, block_threads>>>(primes_.get(), n_primes_);
            cuda_check(cudaGetLastError(), "advance_residues_kernel launch");
            chunk_low += chunk_numbers;
        }
    }

    std::uint64_t count()
    {
        unsigned long long total {0};
        for_each_chunk([&](std::uint64_t, std::uint64_t, std::uint32_t n_words)
        {
            const unsigned long long zero {0};
            cuda_check(cudaMemcpy(totals_.get(), &zero, sizeof(zero), cudaMemcpyHostToDevice), "cudaMemcpy total");
            const std::uint32_t blocks {(n_words + block_threads * 4 - 1) / (block_threads * 4)};
            count_kernel<<<blocks, block_threads>>>(words_.get(), n_words, nullptr, totals_.get());
            cuda_check(cudaGetLastError(), "count_kernel launch");
            unsigned long long chunk_total {0};
            cuda_check(cudaMemcpy(&chunk_total, totals_.get(), sizeof(chunk_total), cudaMemcpyDeviceToHost), "cudaMemcpy total");
            total += chunk_total;
        });
        return total;
    }

    template <class Consumer>
    void generate(Consumer& consume)
    {
        pinned_buffer host;
        for_each_chunk([&](std::uint64_t chunk_low, std::uint64_t, std::uint32_t n_words)
        {
            const std::uint32_t blocks {(n_words + block_threads * 4 - 1) / (block_threads * 4)};
            count_kernel<<<blocks, block_threads>>>(words_.get(), n_words, block_totals_.get(), nullptr);
            cuda_check(cudaGetLastError(), "count_kernel launch");
            scan_totals_kernel<<<1, 1024>>>(block_totals_.get(), blocks, totals_.get());
            cuda_check(cudaGetLastError(), "scan_totals_kernel launch");
            unsigned long long chunk_total {0};
            cuda_check(cudaMemcpy(&chunk_total, totals_.get(), sizeof(chunk_total), cudaMemcpyDeviceToHost), "cudaMemcpy total");
            if (chunk_total == 0)
            {
                return;
            }
            if (output_.size() < chunk_total)
            {
                output_.allocate(static_cast<std::size_t>(chunk_total));
            }
            compact_kernel<<<blocks, block_threads>>>(tables_, words_.get(), n_words, block_totals_.get(), chunk_low, output_.get());
            cuda_check(cudaGetLastError(), "compact_kernel launch");
            host.reserve(static_cast<std::size_t>(chunk_total));
            output_.download(host.data(), static_cast<std::size_t>(chunk_total));
            consume(host.data(), static_cast<std::size_t>(chunk_total));
        });
    }

private:
    std::uint64_t start_;
    std::uint64_t stop_;
    std::uint64_t chunk_low_ {0};
    device_tables tables_;
    device_buffer<device_prime> primes_;
    std::uint32_t n_primes_ {0};
    std::uint32_t n_warp_ {0};
    std::uint32_t n_thread_ {0};
    device_buffer<std::uint8_t> presieve_[3];
    std::uint32_t presieve_period_[3] {};
    device_buffer<std::uint32_t> words_;
    device_buffer<std::uint32_t> block_totals_;
    device_buffer<unsigned long long> totals_;
    device_buffer<std::uint64_t> output_;
};

} // namespace cuda_detail

// Counts the primes in [start, stop] (start >= 7) on the device.
inline std::uint64_t cuda_count(std::uint64_t start, std::uint64_t stop, const prime_sieve_options&)
{
    if (stop < start)
    {
        return 0;
    }
    cuda_detail::cuda_sieve sieve {start, stop};
    return sieve.count();
}

// Passes the primes in [start, stop] (start >= 7) in ascending batches to consume(const std::uint64_t*, std::size_t).
template <class Consumer>
void cuda_range(std::uint64_t start, std::uint64_t stop, const prime_sieve_options&, Consumer& consume)
{
    if (stop < start)
    {
        return;
    }
    cuda_detail::cuda_sieve sieve {start, stop};
    sieve.generate(consume);
}

} // namespace boost::math::detail::prime_sieve

#endif // __CUDACC__ && BOOST_MATH_ENABLE_CUDA && !BOOST_MATH_HAS_NVRTC
#endif // BOOST_MATH_SF_DETAIL_PRIME_SIEVE_CUDA_HPP
