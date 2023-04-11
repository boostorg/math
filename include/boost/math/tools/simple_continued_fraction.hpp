//  (C) Copyright Nick Thompson 2020.
//  (C) Copyright Matt Borland 2023.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_TOOLS_SIMPLE_CONTINUED_FRACTION_HPP
#define BOOST_MATH_TOOLS_SIMPLE_CONTINUED_FRACTION_HPP

#include <array>
#include <vector>
#include <ostream>
#include <iomanip>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <sstream>
#include <utility>
#include <cstdint>
#include <cassert>

#include <boost/math/tools/is_standalone.hpp>
#ifndef BOOST_MATH_STANDALONE
#include <boost/config.hpp>
#ifdef BOOST_NO_CXX17_IF_CONSTEXPR
#error "The header <boost/math/simple_continued_fraction.hpp> can only be used in C++17 and later."
#endif
#endif

#ifndef BOOST_MATH_STANDALONE
#include <boost/core/demangle.hpp>
#endif

namespace boost::math::tools {

template<typename Real, typename Z = int64_t>
class simple_continued_fraction {
public:
    typedef Z int_type;

    simple_continued_fraction(std::vector<Z> data) : b_{std::move(data)} {
        const size_t size_ = b_.size();
        if (size_ == 0) {
            throw std::length_error("Array of coefficients is empty.");
        }

        for (size_t i = 1; i < size_; ++i) {
            if (b_[i] <= 0) {
                std::ostringstream oss;
                oss << "Found a negative partial denominator: b[" << i << "] = " << b_[i] << ".";
                throw std::domain_error(oss.str());
            }
        }

        canonicalize();
    }

    simple_continued_fraction(Real x) : b_{} {
        using std::floor;
        using std::abs;
        using std::sqrt;
        using std::isfinite;
        if (!isfinite(x)) {
            throw std::domain_error("Cannot convert non-finites into continued fractions.");
        }

        if constexpr (std_precision == 2147483647) {
           precision_ = x.backend().precision();
        }

        b_.reserve(50);
        Real bj = floor(x);
        b_.push_back(static_cast<Z>(bj));
        if (bj == x) {
           b_.shrink_to_fit();
           return;
        }

        const Real orig_x = x;
        x = 1/(x-bj);
        Real f = bj;
        if (bj == 0) {
           f = 16*(std::numeric_limits<Real>::min)();
        }
        Real C = f;
        Real D = 0;
        // the "1 + i" lets the error bound grow slowly with the number of convergents.
        // I have not worked out the error propagation of the Modified Lentz's method to see if it does indeed grow at this rate.
        // Numerical Recipes claims that no one has worked out the error analysis of the modified Lentz's method.
        const Real eps_abs_orig_x = std::numeric_limits<Real>::epsilon()*abs(orig_x);
        for (int i = 0; abs(f - orig_x) >= (1 + i)*eps_abs_orig_x; ++i) {
          bj = floor(x);
          b_.push_back(static_cast<Z>(bj));
          x = 1/(x-bj);
          D += bj;
          if (D == 0) {
             D = 16*(std::numeric_limits<Real>::min)();
          }
          C = bj + 1/C;
          if (C==0) {
             C = 16*(std::numeric_limits<Real>::min)();
          }
          D = 1/D;
          f *= (C*D);
       }
       canonicalize();

       const size_t size_ = b_.size();
       for (size_t i = 1; i < size_; ++i) {
         if (b_[i] <= 0) {
            std::ostringstream oss;
            oss << "Found a negative partial denominator: b[" << i << "] = " << b_[i] << "."
                #ifndef BOOST_MATH_STANDALONE
                << " This means the integer type '" << boost::core::demangle(typeid(Z).name())
                #else
                << " This means the integer type '" << typeid(Z).name()
                #endif
                << "' has overflowed and you need to use a wider type,"
                << " or there is a bug.";
            throw std::overflow_error(oss.str());
         }
       }
    }

    Real khinchin_geometric_mean() const {
        const size_t size_ = b_.size();
        if (size_ == 1) {
         return std::numeric_limits<Real>::quiet_NaN();
        }
         using std::log;
         using std::exp;
         // Precompute the most probable logarithms. See the Gauss-Kuzmin distribution for details.
         // Example: b_i = 1 has probability -log_2(3/4) ~ .415:
         // A random partial denominator has ~80% chance of being in this table:
         const std::array<Real, 7> logs{std::numeric_limits<Real>::quiet_NaN(), static_cast<Real>(0), log(static_cast<Real>(2)), log(static_cast<Real>(3)), log(static_cast<Real>(4)), log(static_cast<Real>(5)), log(static_cast<Real>(6))};
         Real log_prod = 0;
         for (size_t i = 1; i < size_; ++i) {
            if (b_[i] < static_cast<Z>(logs.size())) {
               log_prod += logs[b_[i]];
            }
            else
            {
               log_prod += log(static_cast<Real>(b_[i]));
            }
         }
         log_prod /= (size_-1);
         return exp(log_prod);
    }

    Real khinchin_harmonic_mean() const {
        const size_t size_ = b_.size();
        if (size_ == 1) {
          return std::numeric_limits<Real>::quiet_NaN();
        }
        Real n = size_ - 1;
        Real denom = 0;
        for (size_t i = 1; i < size_; ++i) {
            denom += 1/static_cast<Real>(b_[i]);
        }
        return n/denom;
    }

    // Note that this also includes the integer part (i.e. all the coefficients)
    const std::vector<Z>& partial_denominators() const {
      return b_;
    }

    inline std::vector<Z>&& get_data() noexcept {
        return std::move(b_);
    }

    template<typename T, typename Z2>
    friend std::ostream& operator<<(std::ostream& out, simple_continued_fraction<T, Z2>& scf);
private:
    static constexpr int std_precision = std::numeric_limits<Real>::max_digits10;

    void canonicalize() {
        // Deal with non-uniqueness of continued fractions: [a0; a1, ..., an, 1] = a0; a1, ..., an + 1].
        // The shorter representation is considered the canonical representation,
        // so if we compute a non-canonical representation, change it to canonical:
        if (b_.size() > 2 && b_.back() == 1) {
            b_.pop_back();
            b_.back() += 1;
        }
        b_.shrink_to_fit();
    }

    std::vector<Z> b_;

    int precision_{std_precision};
};


template<typename Real, typename Z2>
std::ostream& operator<<(std::ostream& out, simple_continued_fraction<Real, Z2>& scf) {
   out << std::setprecision(scf.precision_);
   out << "[" << scf.b_.front();
   if (scf.b_.size() > 1)
   {
      out << "; ";
      for (size_t i = 1; i < scf.b_.size() -1; ++i)
      {
         out << scf.b_[i] << ", ";
      }
      out << scf.b_.back();
   }
   out << "]";
   return out;
}

template<typename Real, typename Z = std::int64_t>
inline auto simple_continued_fraction_coefficients(Real x)
{
    auto temp = simple_continued_fraction<Real, Z>(x);
    return temp.get_data();
}

// Can be used with `boost::rational` from <boost/rational.hpp>
template <typename Rational, typename Real, typename Z = std::int64_t>
inline Rational to_rational(const simple_continued_fraction<Real, Z>& scf)
{
    using int_t = typename Rational::int_type;

    auto& coefs = scf.partial_denominators();
    const size_t size_ = coefs.size();
    assert(size_ >= 1);
    if (size_ == 1) return static_cast<int_t>(coefs[0]);

    // p0 = a0, p1 = a1.a0 + 1, pn = an.pn-1 + pn-2 for 2 <= n
    // q0 = 1,  q1 = a1,        qn = an.qn-1 + qn-2 for 2 <= n

    int_t p0 = coefs[0];
    int_t p1 = p0*coefs[1] + 1;
    int_t q0 = 1;
    int_t q1 = coefs[1];
    for (size_t i = 2; i < size_; ++i) {
        const Z cn = coefs[i];
        const int_t pn = cn*p1 + p0;
        const int_t qn = cn*q1 + q0;
        p0 = std::exchange(p1, pn);
        q0 = std::exchange(q1, qn);
    }

    return {p1, q1};
}

}
#endif
