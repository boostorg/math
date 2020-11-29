// Copyright 2020 Matt Borland and Jonathan Sorenson
//
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_SPECIAL_FUNCTIONS_DETAIL_PRIME_WHEEL_HPP
#define BOOST_MATH_SPECIAL_FUNCTIONS_DETAIL_PRIME_WHEEL_HPP

#include <array>
#include <memory>
#include <iomanip>
#include <iostream>
#include <cstdint>

namespace boost::math::detail::prime_sieve
{
template<typename Integer>
class Wheel
{    
private:
    struct Wheelrec
    {
        std::int_fast32_t rp;
        std::int_fast32_t dist;
        std::int_fast32_t pos;
        std::int_fast32_t inv;
    };

    std::unique_ptr<Wheelrec[]> W_;
    Integer M_;
    Integer k_;
    Integer phi_;

    static constexpr std::array<Integer, 8> P_ {2, 3, 5, 7, 11, 13, 17, 19};

    void build(Integer korsize);

public:
    Wheel() : W_{nullptr}, M_{0}, k_{0}, phi_{0} {};
    explicit Wheel(Integer korsize) { build(korsize); }
    explicit Wheel(const Wheel &x) { build(x.K()); }

    constexpr bool operator!() const noexcept { return W_ == nullptr; }
    constexpr const Wheelrec& operator[](const Integer i) const noexcept { return W_[i % M_]; }
    const Wheel& operator=(const Wheel &x)
    {
        if(this != &x)
        {
            build(x.K());
        }
        return *this;
    }
    
    constexpr Integer Size() const noexcept { return M_; }
    constexpr Integer K() const noexcept { return k_; }
    constexpr Integer Phi() const noexcept { return phi_; }

    constexpr Integer Next(const Integer i) const noexcept { return i + W_[i % M_].dist; }
    constexpr Integer MakeRP(const Integer i) const noexcept
    {
        if(W_[i % M_].rp)
        {
            return i;
        }
        return Next(i);
    } 
    constexpr Integer Prev(const Integer i) const noexcept { return i - W_[(M_ - (i % M_)) % M_].dist; }
    constexpr Integer Pos(const Integer i) const noexcept { return phi_ * (i / M_) + W_[i % M_].pos; }
    constexpr Integer Inv(const Integer i) const noexcept { return M_ * (i / phi_) + W_[i % phi_].inv; }

    void Print();
};

template<typename Integer>
void Wheel<Integer>::build(Integer korsize)
{
    // Calculate k_ and M_
    if(korsize >= 10)
    {
        --korsize;
        for(k_ = 0; korsize > 0; ++k_)
        {
            korsize /= P_[k_];
        }
    }
    else
    {
        k_ = korsize;
    }

    Integer i {0};
    Integer dist {0};
    Integer pos {1};

    for(M_ = 1; i < k_; ++i)
    {
        M_ *= P_[i];
    }
    
    W_ = std::make_unique<Wheelrec[]>(M_);

    // Compute the RP field
    for(i = 0; i < M_; ++i)
    {
        W_[i].rp = 1;
    }

    for(i = 0; i < k_; ++i)
    {
        for(Integer j {0}; j < M_; j += P_[i])
        {
            W_[j].rp = 0;
        }
    }

    // Compute the dist field
    W_[M_- 1].dist = 2;
    for(i = M_ - 2; i >= 0; --i)
    {
        W_[i].dist = ++dist;
        if(W_[i].rp)
        {
            dist = 0;
        }
    }

    // Copute pos and inv fields
    for(i = 0; i < M_; ++i)
    {
        W_[i].inv = 0;
        if(W_[i].rp)
        {
            W_[pos].inv = i;
            W_[i].pos = pos++;
        }
        else
        {
            W_[i].pos = 0;
        }
        
    }

    W_[0].inv = -1;
    phi_ = W_[M_- 1].pos;
}

template<typename Integer>
void Wheel<Integer>::Print()
{
    std::int_fast32_t i {};
    std::cout << "Wheel size = " << this->Size()
                << "\nk = " << this->K()
                << "\nphi(M) = " << this->Phi() << std::endl;

    // Verify size
    for(i = 0; i < this->Size(); ++i)
    {
        std::cout << std::setw(4) << i << ',';
        if(i % 25 == 24)
        {
            std::cout << std::endl;
        }
    }

    std::cout << "\n\nRP Field\n";
    for(i = 0; i < this->Size(); ++i)
    {
        std::cout << std::setw(3) << W_[i].rp << ',';
        if(i % 25 == 24)
        {
            std::cout << std::endl;
        }
    }

    std::cout << "\n\nDist Field\n";
    for(i = 0; i < this->Size(); ++i)
    {
        std::cout << std::setw(3) << W_[i].dist << ',';
        if(i % 25 == 24)
        {
            std::cout << std::endl;
        }
    }

    std::cout << "\n\nPos Field\n";
    for(i = 0; i < this->Size(); ++i)
    {
        std::cout << std::setw(3) << W_[i].pos << ',';
        if(i % 25 == 24)
        {
            std::cout << std::endl;
        }
    }

    std::cout << "\n\nInv Field\n";
    for(i = 0; i < this->Size(); ++i)
    {
        std::cout << std::setw(4) << W_[i].inv << ',';
        if(i % 25 == 24)
        {
            std::cout << std::endl;
        }
    }
    std::cout << std::endl;
}

// Pre-computed MOD 30 wheel
template <typename Integer>
class MOD30Wheel final
{
    static constexpr std::uint_fast8_t M_ {30};
    static constexpr std::uint_fast8_t k_ {3};
    static constexpr std::uint_fast8_t phi_ {8};

    static constexpr std::array<std::uint_fast8_t, M_> dist_
    {
        1,  6,  5,  4,  3,  2,  1,  4,  3,  2,  1,  2,  1,  4,  3,  2,  1,  2,  1,  4,  3,  2,  1,  6,  5, 4,  
        3,  2,  1,  2
    };
    static constexpr std::array<std::uint_fast8_t, phi_> spokes_
    {
        2, 6, 4, 2, 4, 2, 4, 6
    };

    std::size_t current_index_;
    Integer current_index_num_;

public:
    constexpr MOD30Wheel() = default;
    ~MOD30Wheel() = default;

    constexpr auto Size() const noexcept { return M_; }
    constexpr auto K() const noexcept { return k_; }
    constexpr auto Phi() const noexcept { return phi_; }

    constexpr Integer Next(const Integer i) const noexcept { return i + dist_[static_cast<std::size_t>(i % M_)]; }

    // Avoid using modulus if scanning quickly through the wheel. Uses the stored values to generate the next number
    inline Integer Next() noexcept 
    { 
        ++current_index_;
        if(current_index_ == phi_)
        {
            current_index_ = 0;
        }
        return current_index_num_ += spokes_[current_index_];
    }

    inline auto SetCurrentIndex(const Integer i) noexcept 
    {
        current_index_num_ = Next(i - 1);
        
        const std::size_t temp {static_cast<std::size_t>(current_index_num_ % M_)};
        switch (temp)
        {
        case 0: case 1:
            current_index_ = 0;
            break;
        case 2: case 3: case 4: case 5: case 6: case 7:
            current_index_ = 1;
            break;
        case 8: case 9: case 10: case 11:
            current_index_ = 2;
            break;
        case 12: case 13:
            current_index_ = 3;
            break;
        case 14: case 15: case 16: case 17:
            current_index_ = 4;
            break;
        case 18: case 19:
            current_index_ = 5;
            break;
        case 20: case 21: case 22: case 23:
            current_index_ = 6;
            break;
        case 24: case 25: case 26: case 27: case 28: case 29:
            current_index_ = 7;
            break;
        default:
            break;
        }

        return current_index_num_;
    }

    // Magic number is the number of possible primes in each turn of the mod 30 wheel
    // https://en.wikipedia.org/wiki/Wheel_factorization
    constexpr auto PrimeRatio() const noexcept { return static_cast<double>(phi_) / M_; }
};

// Pre-computed MOD 210 wheel
template <typename Integer>
class MOD210Wheel final
{
private:
    static constexpr auto M_ {210};
    static constexpr auto k_ {4};
    static constexpr auto phi_ {28};
    
    static constexpr std::array<std::int_fast8_t, M_> dist_
    {
        1, 10,  9,  8,  7,  6,  5,  4,  3,  2,  1,  2,  1,  4,  3,  2,  1,  2,  1,  4,  3,  2,  1,  6,  5,
        4,  3,  2,  1,  2,  1,  6,  5,  4,  3,  2,  1,  4,  3,  2,  1,  2,  1,  4,  3,  2,  1,  6,  5,  4,
        3,  2,  1,  6,  5,  4,  3,  2,  1,  2,  1,  6,  5,  4,  3,  2,  1,  4,  3,  2,  1,  2,  1,  6,  5,
        4,  3,  2,  1,  4,  3,  2,  1,  6,  5,  4,  3,  2,  1,  8,  7,  6,  5,  4,  3,  2,  1,  4,  3,  2,
        1,  2,  1,  4,  3,  2,  1,  2,  1,  4,  3,  2,  1,  8,  7,  6,  5,  4,  3,  2,  1,  6,  5,  4,  3,
        2,  1,  4,  3,  2,  1,  6,  5,  4,  3,  2,  1,  2,  1,  4,  3,  2,  1,  6,  5,  4,  3,  2,  1,  2,
        1,  6,  5,  4,  3,  2,  1,  6,  5,  4,  3,  2,  1,  4,  3,  2,  1,  2,  1,  4,  3,  2,  1,  6,  5,
        4,  3,  2,  1,  2,  1,  6,  5,  4,  3,  2,  1,  4,  3,  2,  1,  2,  1,  4,  3,  2,  1,  2,  1, 10,
        9,  8,  7,  6,  5,  4,  3,  2,  1,  2
    };
    
public:
    constexpr MOD210Wheel() = default;    
    ~MOD210Wheel() = default;

    constexpr auto Size() const noexcept { return M_; }
    constexpr auto K() const noexcept { return k_; }
    constexpr auto Phi() const noexcept { return phi_; }

    constexpr auto Next(const Integer i) const noexcept { return i + dist_[static_cast<std::size_t>(i % M_)]; }

    // Magic number is the number of possible primes in each turn of the mod 30 wheel
    // https://en.wikipedia.org/wiki/Wheel_factorization
    constexpr auto PrimeRatio() const noexcept { return 48.0 / 210.0; }
};
}

#endif // BOOST_MATH_SPECIAL_FUNCTIONS_DETAIL_PRIME_WHEEL_HPP
