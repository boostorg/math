// Copyright Nick Thompson, 2020
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_INTERPOLATORS_DETAIL_CUBIC_HERMITE_DETAIL_HPP
#define BOOST_MATH_INTERPOLATORS_DETAIL_CUBIC_HERMITE_DETAIL_HPP
#include <stdexcept>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <sstream>
#include <limits>

#include <boost/math/policies/error_handling.hpp>

namespace boost {
namespace math {
namespace interpolators {
namespace detail {

template<class RandomAccessContainer, class Policy = policies::policy<>>
class cubic_hermite_detail {
public:
    using Real = typename RandomAccessContainer::value_type;
    using Size = typename RandomAccessContainer::size_type;

    cubic_hermite_detail(RandomAccessContainer && x, RandomAccessContainer && y, RandomAccessContainer dydx)
     : x_{std::move(x)}, y_{std::move(y)}, dydx_{std::move(dydx)}
    {
        static const char* function = "boost::math::interpolators::detail::cubic_hermite_detail::cubic_hermite_detail";
        using std::abs;
        using std::isnan;
        if (x_.size() != y_.size())
        {
            error_msg_ = "There must be the same number of ordinates as abscissas.";
            boost::math::policies::raise_domain_error(function, error_msg_.c_str(), false, Policy());
            valid_ = false;
            return;
        }
        if (x_.size() != dydx_.size())
        {
            error_msg_ = "There must be the same number of ordinates as derivative values.";
            boost::math::policies::raise_domain_error(function, error_msg_.c_str(), false, Policy());
            valid_ = false;
            return;
        }
        if (x_.size() < 2)
        {
            error_msg_ = "Must be at least two data points.";
            boost::math::policies::raise_domain_error(function, error_msg_.c_str(), false, Policy());
            valid_ = false;
            return;
        }
        Real x0 = x_[0];
        for (size_t i = 1; i < x_.size(); ++i)
        {
            Real x1 = x_[i];
            if (x1 <= x0)
            {
                std::ostringstream oss;
                oss.precision(std::numeric_limits<Real>::digits10+3);
                oss << "Abscissas must be listed in strictly increasing order x0 < x1 < ... < x_{n-1}, ";
                oss << "but at x[" << i - 1 << "] = " << x0 << ", and x[" << i << "] = " << x1 << ".\n";
                error_msg_ = oss.str();
                valid_ = false;
                return;
            }
            x0 = x1;
        }
        valid_ = true;
    }

    void push_back(Real x, Real y, Real dydx)
    {
        if ( ! valid_) return;  // We do not realize if the object is constructed correctly here.

        static const char* function = "boost::math::interpolators::detail::cubic_hermite_detail::push_back";
        using std::abs;
        using std::isnan;
        if (x <= x_.back())
        {
            return boost::math::policies::raise_domain_error<bool>(function, "Calling push_back must preserve the monotonicity of the x's", false, Policy());
        }
        x_.push_back(x);
        y_.push_back(y);
        dydx_.push_back(dydx);
    }

    Real operator()(Real x) const
    {
        if ( ! valid_) return std::numeric_limits<Real>::quiet_NaN();

        static const char* function = "boost::math::interpolators::detail::cubic_hermite_detail::operator()";
        if  (x < x_[0] || x > x_.back())
        {
            std::ostringstream oss;
            oss.precision(std::numeric_limits<Real>::digits10+3);
            oss << "Requested abscissa x = " << x << ", which is outside of allowed range ["
                << x_[0] << ", " << x_.back() << "]";
            return boost::math::policies::raise_domain_error<Real>(function, oss.str().c_str(), x, Policy());
        }
        // We need t := (x-x_k)/(x_{k+1}-x_k) \in [0,1) for this to work.
        // Sadly this neccessitates this loathesome check, otherwise we get t = 1 at x = xf.
        if (x == x_.back())
        {
            return y_.back();
        }

        auto it = std::upper_bound(x_.begin(), x_.end(), x);
        auto i = std::distance(x_.begin(), it) -1;
        Real x0 = *(it-1);
        Real x1 = *it;
        Real y0 = y_[i];
        Real y1 = y_[i+1];
        Real s0 = dydx_[i];
        Real s1 = dydx_[i+1];
        Real dx = (x1-x0);
        Real t = (x-x0)/dx;

        // See the section 'Representations' in the page
        // https://en.wikipedia.org/wiki/Cubic_Hermite_spline
        Real y = (1-t)*(1-t)*(y0*(1+2*t) + s0*(x-x0))
              + t*t*(y1*(3-2*t) + dx*s1*(t-1));
        return y;
    }

    Real prime(Real x) const
    {
        if ( ! valid_) return std::numeric_limits<Real>::quiet_NaN();

        static const char* function = "boost::math::interpolators::detail::cubic_hermite_detail::prime";
        if  (x < x_[0] || x > x_.back())
        {
            std::ostringstream oss;
            oss.precision(std::numeric_limits<Real>::digits10+3);
            oss << "Requested abscissa x = " << x << ", which is outside of allowed range ["
                << x_[0] << ", " << x_.back() << "]";
            return boost::math::policies::raise_domain_error<Real>(function, oss.str().c_str(), x, Policy());
        }
        if (x == x_.back())
        {
            return dydx_.back();
        }
        auto it = std::upper_bound(x_.begin(), x_.end(), x);
        auto i = std::distance(x_.begin(), it) -1;
        Real x0 = *(it-1);
        Real x1 = *it;
        Real y0 = y_[i];
        Real y1 = y_[i+1];
        Real s0 = dydx_[i];
        Real s1 = dydx_[i+1];
        Real dx = (x1-x0);

        Real d1 = (y1 - y0 - s0*dx)/(dx*dx);
        Real d2 = (s1 - s0)/(2*dx);
        Real c2 = 3*d1 - 2*d2;
        Real c3 = 2*(d2 - d1)/dx;
        return s0 + 2*c2*(x-x0) + 3*c3*(x-x0)*(x-x0); 
    }

    friend std::ostream& operator<<(std::ostream & os, const cubic_hermite_detail & m)
    {
        if ( ! m.valid_) return os;

        os << "(x,y,y') = {";
        for (size_t i = 0; i < m.x_.size() - 1; ++i)
        {
            os << "(" << m.x_[i] << ", " << m.y_[i] << ", " << m.dydx_[i] << "),  ";
        }
        auto n = m.x_.size()-1;
        os << "(" << m.x_[n] << ", " << m.y_[n] << ", " << m.dydx_[n] << ")}";
        return os;
    }

    Size size() const
    {
        if ( ! valid_) return 0;
        return x_.size();
    }

    int64_t bytes() const
    {
        if ( ! valid_) return 0;
        return 3*x_.size()*sizeof(Real) + 3*sizeof(x_);
    }

    std::pair<Real, Real> domain() const
    {
        if ( ! valid_) return {0, 0};
        return {x_.front(), x_.back()};
    }

    bool valid() const {
        return valid_;
    }

    std::string const& error_msg() const {
        return error_msg_;
    }

    RandomAccessContainer x_;
    RandomAccessContainer y_;
    RandomAccessContainer dydx_;
private:
    bool valid_ = false;
    std::string error_msg_;
};

template<class RandomAccessContainer, class Policy = policies::policy<>>
class cardinal_cubic_hermite_detail {
public:
    using Real = typename RandomAccessContainer::value_type;
    using Size = typename RandomAccessContainer::size_type;

    cardinal_cubic_hermite_detail(RandomAccessContainer && y, RandomAccessContainer dydx, Real x0, Real dx)
    : y_{std::move(y)}, dy_{std::move(dydx)}, x0_{x0}, inv_dx_{1/dx}
    {
        static const char* function = "boost::math::interpolators::detail::cardinal_cubic_hermite_detail::cardinal_cubic_hermite_detail";
        using std::abs;
        using std::isnan;
        if (y_.size() != dy_.size())
        {
            error_msg_ = "There must be the same number of derivatives as ordinates.";
            boost::math::policies::raise_domain_error(function, error_msg_.c_str(), false, Policy());
            valid_ = false;
            return;
        }
        if (y_.size() < 2)
        {
            error_msg_ = "Must be at least two data points.";
            boost::math::policies::raise_domain_error(function, error_msg_.c_str(), false, Policy());
            valid_ = false;
            return;
        }
        if (dx <= 0)
        {
            error_msg_ = "dx > 0 is required.";
            boost::math::policies::raise_domain_error(function, error_msg_.c_str(), false, Policy());
            valid_ = false;
            return;
        }

        for (auto & dy : dy_)
        {
            dy *= dx;
        }
        valid_ = true;
    }

    // Why not implement push_back? It's awkward: If the buffer is circular, x0_ += dx_.
    // If the buffer is not circular, x0_ is unchanged.
    // We need a concept for circular_buffer!

    inline Real operator()(Real x) const
    {
        static const char* function = "boost::math::interpolators::detail::cardinal_cubic_hermite_detail::operator()";
        if ( ! valid_) return std::numeric_limits<Real>::quiet_NaN();

        const Real xf = x0_ + (y_.size()-1)/inv_dx_;
        if  (x < x0_ || x > xf)
        {
            std::ostringstream oss;
            oss.precision(std::numeric_limits<Real>::digits10+3);
            oss << "Requested abscissa x = " << x << ", which is outside of allowed range ["
                << x0_ << ", " << xf << "]";
            return boost::math::policies::raise_domain_error<Real>(function, oss.str().c_str(), x, Policy());
        }
        if (x == xf)
        {
            return y_.back();
        }
        return this->unchecked_evaluation(x);
    }

    inline Real unchecked_evaluation(Real x) const
    {
        if ( ! valid_) return std::numeric_limits<Real>::quiet_NaN();
        using std::floor;
        Real s = (x-x0_)*inv_dx_;
        Real ii = floor(s);
        auto i = static_cast<decltype(y_.size())>(ii);
        Real t = s - ii;
        Real y0 = y_[i];
        Real y1 = y_[i+1];
        Real dy0 = dy_[i];
        Real dy1 = dy_[i+1];

        Real r = 1-t;
        return r*r*(y0*(1+2*t) + dy0*t)
              + t*t*(y1*(3-2*t) - dy1*r);
    }

    inline Real prime(Real x) const
    {
        if ( ! valid_) return std::numeric_limits<Real>::quiet_NaN();
        static const char* function = "boost::math::interpolators::detail::cardinal_cubic_hermite_detail::prime";
        const Real xf = x0_ + (y_.size()-1)/inv_dx_;
        if  (x < x0_ || x > xf)
        {
            std::ostringstream oss;
            oss.precision(std::numeric_limits<Real>::digits10+3);
            oss << "Requested abscissa x = " << x << ", which is outside of allowed range ["
                << x0_ << ", " << xf << "]";
            return boost::math::policies::raise_domain_error<Real>(function, oss.str().c_str(), x, Policy());
        }
        if (x == xf)
        {
            return dy_.back()*inv_dx_;
        }
        return this->unchecked_prime(x);
    }

    inline Real unchecked_prime(Real x) const
    {
        if ( ! valid_) return std::numeric_limits<Real>::quiet_NaN();

        using std::floor;
        Real s = (x-x0_)*inv_dx_;
        Real ii = floor(s);
        auto i = static_cast<decltype(y_.size())>(ii);
        Real t = s - ii;
        Real y0 = y_[i];
        Real y1 = y_[i+1];
        Real dy0 = dy_[i];
        Real dy1 = dy_[i+1];

        Real dy = 6*t*(1-t)*(y1 - y0)  + (3*t*t - 4*t+1)*dy0 + t*(3*t-2)*dy1;
        return dy*inv_dx_;
    }

    Size size() const
    {
        if ( ! valid_) return 0;
        return y_.size();
    }

    int64_t bytes() const
    {
        if ( ! valid_) return 0;
        return 2*y_.size()*sizeof(Real) + 2*sizeof(y_) + 2*sizeof(Real);
    }

    std::pair<Real, Real> domain() const
    {
        if ( ! valid_) return {0, 0};
        Real xf = x0_ + (y_.size()-1)/inv_dx_;
        return {x0_, xf};
    }

    bool valid() const {
        return valid_;
    }

    std::string const& error_msg() const {
        return error_msg_;
    }
private:
    RandomAccessContainer y_;
    RandomAccessContainer dy_;
    Real x0_;
    Real inv_dx_;

    bool valid_ = false;
    std::string error_msg_;
};


template<class RandomAccessContainer, class Policy = policies::policy<>>
class cardinal_cubic_hermite_detail_aos {
public:
    using Point = typename RandomAccessContainer::value_type;
    using Real = typename Point::value_type;
    using Size = typename RandomAccessContainer::size_type;

    cardinal_cubic_hermite_detail_aos(RandomAccessContainer && dat, Real x0, Real dx)
        : dat_{std::move(dat)}, x0_{x0}, inv_dx_{1/dx}
    {
        static const char* function = "boost::math::interpolators::detail::cardinal_cubic_hermite_detail_aos::cardinal_cubic_hermite_detail_aos";
        if (dat_.size() < 2)
        {
            error_msg_ = "Must be at least two data points.";
            boost::math::policies::raise_domain_error(function, error_msg_.c_str(), false, Policy());
            valid_ = false;
            return;
        }
        if (dat_[0].size() != 2)
        {
            error_msg_ = "Each datum must contain (y, y'), and nothing else.";
            boost::math::policies::raise_domain_error(function, error_msg_.c_str(), false, Policy());
            valid_ = false;
            return;
        }
        if (dx <= 0)
        {
            error_msg_ = "dx > 0 is required.";
            boost::math::policies::raise_domain_error(function, error_msg_.c_str(), false, Policy());
            valid_ = false;
            return;
        }

        for (auto & d : dat_)
        {
            d[1] *= dx;
        }
        valid_ = true;
    }

    inline Real operator()(Real x) const
    {
        if ( ! valid_) return std::numeric_limits<Real>::quiet_NaN();

        static const char* function = "boost::math::interpolators::detail::cardinal_cubic_hermite_detail_aos::operator()";
        const Real xf = x0_ + (dat_.size()-1)/inv_dx_;
        if  (x < x0_ || x > xf)
        {
            std::ostringstream oss;
            oss.precision(std::numeric_limits<Real>::digits10+3);
            oss << "Requested abscissa x = " << x << ", which is outside of allowed range ["
                << x0_ << ", " << xf << "]";
            return boost::math::policies::raise_domain_error<Real>(function, oss.str().c_str(), x, Policy());
        }
        if (x == xf)
        {
            return dat_.back()[0];
        }
        return this->unchecked_evaluation(x);
    }

    inline Real unchecked_evaluation(Real x) const
    {
        if ( ! valid_) return std::numeric_limits<Real>::quiet_NaN();

        using std::floor;
        Real s = (x-x0_)*inv_dx_;
        Real ii = floor(s);
        auto i = static_cast<decltype(dat_.size())>(ii);

        Real t = s - ii;
        // If we had infinite precision, this would never happen.
        // But we don't have infinite precision.
        if (t == 0)
        {
            return dat_[i][0];
        }
        Real y0 = dat_[i][0];
        Real y1 = dat_[i+1][0];
        Real dy0 = dat_[i][1];
        Real dy1 = dat_[i+1][1];

        Real r = 1-t;
        return r*r*(y0*(1+2*t) + dy0*t)
              + t*t*(y1*(3-2*t) - dy1*r);
    }

    inline Real prime(Real x) const
    {
        if ( ! valid_) return std::numeric_limits<Real>::quiet_NaN();

        static const char* function = "boost::math::interpolators::detail::cardinal_cubic_hermite_detail_aos::prime";
        const Real xf = x0_ + (dat_.size()-1)/inv_dx_;
        if  (x < x0_ || x > xf)
        {
            std::ostringstream oss;
            oss.precision(std::numeric_limits<Real>::digits10+3);
            oss << "Requested abscissa x = " << x << ", which is outside of allowed range ["
                << x0_ << ", " << xf << "]";
            return boost::math::policies::raise_domain_error<Real>(function, oss.str().c_str(), x, Policy());
        }
        if (x == xf)
        {
            return dat_.back()[1]*inv_dx_;
        }
        return this->unchecked_prime(x);
    }

    inline Real unchecked_prime(Real x) const
    {
        if ( ! valid_) return std::numeric_limits<Real>::quiet_NaN();

        using std::floor;
        Real s = (x-x0_)*inv_dx_;
        Real ii = floor(s);
        auto i = static_cast<decltype(dat_.size())>(ii);
        Real t = s - ii;
        if (t == 0)
        {
            return dat_[i][1]*inv_dx_;
        }
        Real y0 = dat_[i][0];
        Real dy0 = dat_[i][1];
        Real y1 = dat_[i+1][0];
        Real dy1 = dat_[i+1][1];

        Real dy = 6*t*(1-t)*(y1 - y0)  + (3*t*t - 4*t+1)*dy0 + t*(3*t-2)*dy1;
        return dy*inv_dx_;
    }

    Size size() const
    {
        if ( ! valid_) return 0;
        return dat_.size();
    }

    int64_t bytes() const
    {
        if ( ! valid_) return 0;
        return dat_.size()*dat_[0].size()*sizeof(Real) + sizeof(dat_) + 2*sizeof(Real);
    }

    std::pair<Real, Real> domain() const
    {
        if ( ! valid_) return {0, 0};
        Real xf = x0_ + (dat_.size()-1)/inv_dx_;
        return {x0_, xf};
    }

    bool valid() const {
        return valid_;
    }

    std::string const& error_msg() const {
        return error_msg_;
    }

private:
    RandomAccessContainer dat_;
    Real x0_;
    Real inv_dx_;

    bool valid_ = false;
    std::string error_msg_;
};

}
}
}
}
#endif

