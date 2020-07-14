//  (C) Copyright Nick Thompson 2020.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_TOOLS_LUROTH_EXPANSION_HPP
#define BOOST_MATH_TOOLS_LUROTH_EXPANSION_HPP

#include <vector>
#include <ostream>
#include <iomanip>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace boost::math::tools {

template<typename Real, typename Z = int64_t>
class luroth_expansion {
public:
    luroth_expansion(Real x) : x_{x}
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
        if (dn == x) {
           d_.shrink_to_fit();
           return;
        }
        // This attempts to follow the notation of:
        // "Khinchine's constant for Luroth Representation", by Sophia Kalpazidou.
        x = x - dn;
        Real computed = dn;
        Real prod = 1;
        // Why scale by 2eps? Well, I tried it and this is a choice that rarely prints an incorrect digit.
        while (abs(x_ - computed) > 2*std::numeric_limits<Real>::epsilon()*abs(x_))
        {
           Real dn = floor(1/x);
           d_.push_back(static_cast<Z>(dn));
           Real tmp = 1/(dn+1);
           computed += prod*tmp;
           prod *= tmp/dn;
           x = dn*(dn+1)*(x - tmp);
        }
    }
    
    
    const std::vector<Z>& digits() const {
      return d_;
    }
    
    template<typename T, typename Z2>
    friend std::ostream& operator<<(std::ostream& out, luroth_expansion<T, Z2>& scf);

private:
    const Real x_;
    std::vector<Z> d_;
};


template<typename Real, typename Z2>
std::ostream& operator<<(std::ostream& out, luroth_expansion<Real, Z2>& luroth)
{
   constexpr const int p = std::numeric_limits<Real>::max_digits10;
   if constexpr (p == 2147483647)
   {
      out << std::setprecision(luroth.x_.backend().precision());
   }
   else
   {
      out << std::setprecision(p);
   }
   
   out << "((" << luroth.d_.front();
   if (luroth.d_.size() > 1)
   {
      out << "; ";
      for (size_t i = 1; i < luroth.d_.size() -1; ++i)
      {
         out << luroth.d_[i] << ", ";
      }
      out << luroth.d_.back();
   }
   out << "))";
   return out;
}


}
#endif
