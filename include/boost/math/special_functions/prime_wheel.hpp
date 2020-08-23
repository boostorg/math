// Copyright 2020 Matt Borland and Jonathan Sorenson
//
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_SPECIAL_FUNCTIONS_PRIME_WHEEL_HPP
#define BOOST_MATH_SPECIAL_FUNCTIONS_PRIME_WHEEL_HPP

#include <array>
#include <memory>
#include <iomanip>
#include <iostream>
#include <cstdint>

namespace boost::math::detail
{
template<class Integer>
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
        if(W[i % M].rp)
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

template<class Integer>
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

template<class Integer>
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

// Pre-computed MOD 210 wheel
template <class Integer>
class MOD210Wheel final
{
private:
    static constexpr auto M_ {210};
    static constexpr auto k_ {4};
    static constexpr auto phi_ {28};

    static constexpr std::array<std::int_fast8_t, M_> rp_
    {
        0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  1,  0,  0,  0,  1,  0,  1,  0,  0,  0,  1,  0,
        0,  0,  0,  0,  1,  0,  1,  0,  0,  0,  0,  0,  1,  0,  0,  0,  1,  0,  1,  0,  0,  0,  1,  0,  0,
        0,  0,  0,  1,  0,  0,  0,  0,  0,  1,  0,  1,  0,  0,  0,  0,  0,  1,  0,  0,  0,  1,  0,  1,  0,
        0,  0,  0,  0,  1,  0,  0,  0,  1,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,
        0,  1,  0,  1,  0,  0,  0,  1,  0,  1,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,
        0,  0,  1,  0,  0,  0,  1,  0,  0,  0,  0,  0,  1,  0,  1,  0,  0,  0,  1,  0,  0,  0,  0,  0,  1,
        0,  1,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  1,  0,  0,  0,  1,  0,  1,  0,  0,  0,  1,  0,
        0,  0,  0,  0,  1,  0,  1,  0,  0,  0,  0,  0,  1,  0,  0,  0,  1,  0,  1,  0,  0,  0,  1,  0,  1,
        0,  0,  0,  0,  0,  0,  0,  0,  0,  1
    };

    static constexpr std::array<std::int_fast16_t, M_> inv_
    {
       -1,   1,  11,  13,  17,  19,  23,  29,  31,  37,  41,  43,  47,  53,  59,  61,  67,  71,  73,  79,  83,  89,  97, 101, 103,
      107, 109, 113, 121, 127, 131, 137, 139, 143, 149, 151, 157, 163, 167, 169, 173, 179, 181, 187, 191, 193, 197, 199, 209,   0,
        0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
        0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
        0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
        0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
        0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
        0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
        0,   0,   0,   0,   0,   0,   0,   0,   0,   0
    };

    static constexpr std::array<std::int_fast8_t, M_> pos_
    {
        0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  0,  3,  0,  0,  0,  4,  0,  5,  0,  0,  0,  6,  0,
        0,  0,  0,  0,  7,  0,  8,  0,  0,  0,  0,  0,  9,  0,  0,  0, 10,  0, 11,  0,  0,  0, 12,  0,  0,
        0,  0,  0, 13,  0,  0,  0,  0,  0, 14,  0, 15,  0,  0,  0,  0,  0, 16,  0,  0,  0, 17,  0, 18,  0,
        0,  0,  0,  0, 19,  0,  0,  0, 20,  0,  0,  0,  0,  0, 21,  0,  0,  0,  0,  0,  0,  0, 22,  0,  0,
        0, 23,  0, 24,  0,  0,  0, 25,  0, 26,  0,  0,  0, 27,  0,  0,  0,  0,  0,  0,  0, 28,  0,  0,  0,
        0,  0, 29,  0,  0,  0, 30,  0,  0,  0,  0,  0, 31,  0, 32,  0,  0,  0, 33,  0,  0,  0,  0,  0, 34,
        0, 35,  0,  0,  0,  0,  0, 36,  0,  0,  0,  0,  0, 37,  0,  0,  0, 38,  0, 39,  0,  0,  0, 40,  0,
        0,  0,  0,  0, 41,  0, 42,  0,  0,  0,  0,  0, 43,  0,  0,  0, 44,  0, 45,  0,  0,  0, 46,  0, 47,
        0,  0,  0,  0,  0,  0,  0,  0,  0, 48
    };

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

    constexpr auto Next(const Integer i) const noexcept { return i + dist_[i % M_]; }
    constexpr auto MakeRP(const Integer i) const noexcept
    {
        if(rp_[i % M])
        {
            return i;
        }
        return Next(i);
    } 
    constexpr auto Prev(const Integer i) const noexcept { return i - dist_[(M_ - (i % M_)) % M_]; }
    constexpr auto Pos(const Integer i) const noexcept { return phi_ * (i / M_) + W_[i % M_].pos; }
    constexpr auto Inv(const Integer i) const noexcept { return M_ * (i / phi_) + W_[i % phi_].inv; }
};
}

#endif //BOOST_MATH_SPECIAL_FUNCTIONS_PRIME_WHEEL_HPP
