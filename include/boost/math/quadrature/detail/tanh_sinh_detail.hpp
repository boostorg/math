// Copyright Nick Thompson, 2017
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_QUADRATURE_DETAIL_TANH_SINH_DETAIL_HPP
#define BOOST_MATH_QUADRATURE_DETAIL_TANH_SINH_DETAIL_HPP

#include <cmath>
#include <vector>
#include <typeinfo>
#include <boost/math/constants/constants.hpp>
#include <boost/math/quadrature/detail/tanh_sinh_constants.hpp>
#include <boost/math/special_functions/next.hpp>


namespace boost{ namespace math{ namespace quadrature { namespace detail{


// Returns the tanh-sinh quadrature of a function f over the open interval (-1, 1)

template<class Real, class Policy>
class tanh_sinh_detail 
   : public tanh_sinh_detail_constants<Real, 
       !std::numeric_limits<Real>::is_specialized || (std::numeric_limits<Real>::radix != 2) ? 0 :
       (std::numeric_limits<Real>::digits <= 24) ? 24 :
       (std::numeric_limits<Real>::digits <= 64) ? 64 : 
#ifdef BOOST_HAS_FLOAT128
       (std::numeric_limits<Real>::digits <= 113) ? 113 : 
#endif
       (std::numeric_limits<Real>::digits <= 334) ? 334 :
       0>
{
   typedef tanh_sinh_detail_constants<Real, !std::numeric_limits<Real>::is_specialized || (std::numeric_limits<Real>::radix != 2) ? 0 :
      std::numeric_limits<Real>::digits <= std::numeric_limits<float>::digits ? 24 : 
      (std::numeric_limits<Real>::digits <= 64) ? 64 : 
#ifdef BOOST_HAS_FLOAT128
      (std::numeric_limits<Real>::digits <= 113) ? 113 :
#endif
      (std::numeric_limits<Real>::digits <= 334) ? 334 :
      0> base_type;
public:
    tanh_sinh_detail(Real tol, size_t max_refinements);

    template<class F>
    Real integrate(const F f, Real* error, Real* L1, const char* function, Real left_min_complement = tools::epsilon<Real>(), Real right_min_complement = tools::epsilon<Real>()) const;

private:
    Real m_tol;
};

template<class Real, class Policy>
tanh_sinh_detail<Real, Policy>::tanh_sinh_detail(Real tol, std::size_t max_refinements) : tanh_sinh_detail_constants(max_refinements), m_tol(tol) {}

template<class Real, class Policy>
template<class F>
Real tanh_sinh_detail<Real, Policy>::integrate(const F f, Real* error, Real* L1, const char* function, Real left_min_complement, Real right_min_complement) const
{
    using std::abs;
    using std::floor;
    using std::tanh;
    using std::sinh;
    using std::sqrt;
    using boost::math::constants::half;
    using boost::math::constants::half_pi;
    //
    // We maintain 4 integer values:
    // max_left_position is the logical index of the abscissa value closest to the
    // left endpoint of the range that we can call f(x_i) on without rounding error
    // inside f(x_i) causing evaluation at the endpoint.
    // max_left_index is the actual position in the current row that has a logical index
    // no higher than max_left_position.  Remember that since we only store odd numbered
    // indexes in each row, this may actually be one position to the left of max_left_position
    // in the case that is even.  Then, if we only evaluate f(-x_i) for abscissa values 
    // i <= max_left_index we will never evaluate f(-x_i) at the left endpoint.
    // max_right_position and max_right_index are defained similarly for the right boundary
    // and are used to guard evaluation of f(x_i).
    //
    // max_left_position and max_right_position start off as the last element in row zero:
    // 
    std::size_t max_left_position(m_abscissas[0].size() - 1);
    std::size_t max_left_index, max_right_position(max_left_position), max_right_index;
    //
    // Decrement max_left_position and max_right_position until the complement 
    // of the abscissa value is greater than the smallest permitted (as specified
    // by the function caller):
    //
    while (max_left_position && fabs(m_abscissas[0][max_left_position]) < left_min_complement)
       --max_left_position;
    while (max_right_position && fabs(m_abscissas[0][max_right_position]) < right_min_complement)
       --max_right_position;
    //
    // Assumption: left_min_complement/right_min_complement are sufficiently small that we only
    // ever decrement through the stored values that are complements (the negative ones), and 
    // never ever hit the true abscissa values (positive stored values).
    //
    BOOST_ASSERT(m_abscissas[0][max_left_position] < 0);
    BOOST_ASSERT(m_abscissas[0][max_right_position] < 0);


    Real h = m_t_max / m_inital_row_length;
    Real I0 = half_pi<Real>()*f(0, 1);
    Real L1_I0 = abs(I0);
    for(size_t i = 1; i < m_abscissas[0].size(); ++i)
    {
        if ((i > max_right_position) && (i > max_left_position))
            break;
        Real x = m_abscissas[0][i];
        Real xc = x;
        Real w = m_weights[0][i];
        if ((boost::math::signbit)(x))
        {
           // We have stored x - 1:
           x = 1 + xc;
        }
        else
           xc = x - 1;
        Real yp, ym;
        yp = i <= max_right_position ? f(x, -xc) : 0;
        ym = i <= max_left_position ? f(-x, xc) : 0;
        I0 += (yp + ym)*w;
        L1_I0 += (abs(yp) + abs(ym))*w;
    }
    //
    // We have:
    // k = current row.
    // I0 = last integral value.
    // I1 = current integral value.
    // L1_I0 and L1_I1 are the absolute integral values.
    //
    size_t k = 1;
    Real I1 = I0;
    Real L1_I1 = L1_I0;
    Real err = 0;

    while (k < 4 || (k < m_weights.size() && k < m_max_refinements) )
    {
        I0 = I1;
        L1_I0 = L1_I1;

        I1 = half<Real>()*I0;
        L1_I1 = half<Real>()*L1_I0;
        h *= half<Real>();
        Real sum = 0;
        Real absum = 0;
        auto const& abscissa_row = this->get_abscissa_row(k);
        auto const& weight_row = this->get_weight_row(k);
        std::size_t first_complement_index = this->get_first_complement_index(k);
        //
        // At the start of each new row we need to update the max left/right indexes
        // at which we can evaluate f(x_i).  The new logical position is simply twice
        // the old value.  The new max index is one position to the left of the new
        // logical value (remmember each row contains only odd numbered positions).
        // Then we have to make a single check, to see if one position to the right
        // is also in bounds (this is the new abscissa value in this row which is 
        // known to be in between a value known to be in bounds, and one known to be
        // not in bounds).
        // Thus, we filter which abscissa values generate a call to f(x_i), with a single
        // floating point comparison per loop.  Everything else is integer logic.
        //
        max_left_index = max_left_position - 1;
        max_left_position *= 2;
        max_right_index = max_right_position - 1;
        max_right_position *= 2;
        if ((abscissa_row.size() > max_left_index + 1) && (fabs(abscissa_row[max_left_index + 1]) > left_min_complement))
        {
           ++max_left_position;
           ++max_left_index;
        }
        if ((abscissa_row.size() > max_right_index + 1) && (fabs(abscissa_row[max_right_index + 1]) > right_min_complement))
        {
           ++max_right_position;
           ++max_right_index;
        }

        for(size_t j = 0; j < weight_row.size(); ++j)
        {
            // If both left and right abscissa values are out of bounds at this step
            // we can just stop this loop right now:
            if ((j > max_left_index) && (j > max_right_index))
                break;
            Real x = abscissa_row[j];
            Real xc = x;
            Real w = weight_row[j];
            if (j >= first_complement_index)
            {
               // We have stored x - 1:
               BOOST_ASSERT(x < 0);
               x = 1 + xc;
            }
            else
            {
               BOOST_ASSERT(x >= 0);
               xc = x - 1;
            }

            Real yp = j > max_right_index ? 0 : f(x, -xc);
            Real ym = j > max_left_index ? 0 : f(-x, xc);
            Real term = (yp + ym)*w;
            sum += term;

            // A question arises as to how accurately we actually need to estimate the L1 integral.
            // For simple integrands, computing the L1 norm makes the integration 20% slower,
            // but for more complicated integrands, this calculation is not noticeable.
            Real abterm = (abs(yp) + abs(ym))*w;
            absum += abterm;
        }

        I1 += sum*h;
        L1_I1 += absum*h;
        ++k;
        err = abs(I0 - I1);
        // std::cout << "Estimate:        " << I1 << " Error estimate at level " << k  << " = " << err << std::endl;

        if (!(boost::math::isfinite)(I1))
        {
            return policies::raise_evaluation_error(function, "The tanh_sinh quadrature evaluated your function at a singular point at got %1%. Please narrow the bounds of integration or check your function for singularities.", I1, Policy());
        }
        //
        // Termination condition:
        // No more levels are considered once the error is less than the specified tolerance.
        // Note however, that we always go down at least 4 levels, otherwise we risk missing
        // features of interest in f() - imagine for example a function which flatlines, except
        // for a very small "spike". An example would be the incomplete beta integral with large
        // parameters.  We could keep hunting until we find something, but that would handicap
        // integrals which really are zero.... so a compromise then!
        //
        if ((k > 4) && (err <= m_tol*L1_I1))
        {
            break;
        }

    }
    if (error)
    {
        *error = err;
    }

    if (L1)
    {
        *L1 = L1_I1;
    }

    return I1;
}

}}}}  // namespaces
#endif
