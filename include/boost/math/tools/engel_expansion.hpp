//  (C) Copyright Nick Thompson 2020.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_TOOLS_ENGEL_EXPANSION_HPP
#define BOOST_MATH_TOOLS_ENGEL_EXPANSION_HPP

#include <vector>
#include <ostream>
#include <iomanip>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace boost::math::tools {

template<typename Real, typename Z = int64_t>
class engel_expansion {
public:
    engel_expansion(Real x) : x_{x}
    {
        using std::floor;
        using std::abs;
        using std::sqrt;
        using std::isfinite;
        if (!isfinite(x))
        {
            throw std::domain_error("Cannot convert non-finites into a Lüroth representation.");
        }
        d_.reserve(50);
        Real dn = floor(x);
        d_.push_back(static_cast<Z>(dn));
        if (dn == x)
        {
           d_.shrink_to_fit();
           return;
        }
        x = x - dn;
        Real computed = dn;
        Real prod = 1;
        // Let the error bound grow by 1 ULP/iteration.
        // I haven't done the error analysis to show that this is an expected rate of error growth,
        // but if you don't do this, you can easily get into an infinite loop.
        Real i = 1;
        Real scale = std::numeric_limits<Real>::epsilon()*abs(x_)/2;
        while (abs(x_ - computed) > (i++)*scale)
        {
           Real recip = 1/x;
           Real dn = floor(recip);
           if (recip == dn) {
              d_.push_back(static_cast<Z>(dn - 1));
              break;
           }
           d_.push_back(static_cast<Z>(dn));
           Real tmp = 1/(dn+1);
           computed += prod*tmp;
           prod *= tmp/dn;
           x = dn*(dn+1)*(x - tmp);
        }

        for (size_t i = 1; i < d_.size(); ++i)
        {
            // Sanity check:
            if (d_[i] < d_[i-1])
            {
                throw std::domain_error("The digits of an Engel expansion must form a non-decreasing sequence.");
            }
        }
        d_.shrink_to_fit();
    }
    
    
    const std::vector<Z>& digits() const {
      return d_;
    }

    template<typename T, typename Z2>
    friend std::ostream& operator<<(std::ostream& out, engel_expansion<T, Z2>& eng);

private:
    const Real x_;
    std::vector<Z> d_;
};


template<typename Real, typename Z2>
std::ostream& operator<<(std::ostream& out, engel_expansion<Real, Z2>& engel)
{
    constexpr const int p = std::numeric_limits<Real>::max_digits10;
    if constexpr (p == 2147483647)
    {
        out << std::setprecision(engel.x_.backend().precision());
    }
    else
    {
        out << std::setprecision(p);
    }

    out << "{";
    for (size_t i = 0; i < engel.d_.size() - 1; ++i)
    {
        out << engel.d_[i] << ", ";
    }
    out << engel.d_.back();
    out << "}";
    return out;
}


}
#endif
