//  (C) Copyright Nick Thompson 2018.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_SPECIAL_KUMMER_M_HPP
#define BOOST_MATH_SPECIAL_KUMMER_M_HPP
#include <boost/math/tools/promotion.hpp>

namespace boost { namespace math {

template<class T1, class T2, class T3>
inline typename tools::promote_args<T1, T2, T3>::type kummer_m(T1 const & a, T2 const & b, T3 const & x)
{
    using std::floor;

    if (floor(b) == b && b <= 0)
    {
        throw std::logic_error("The Kummer-M function M(a,b,x) is not defined when b is a negative integer.\n");
    }
    if(abs(x) >= 1)
    {
        throw std::logic_error("not yet implemented\n");
    }
    typedef typename boost::math::tools::promote_args<T1, T2, T3>::type Real;


    Real term = 1;
    Real result = term;
    Real an = a;
    Real bn = b;
    size_t n = 1;
    int consecutive_small_terms = 0;
    while (consecutive_small_terms < 3)
    {
        term *= an*x/(bn*n);
        an += 1;
        bn += 1;
        ++n;
        result += term;
        if (abs(term) <= std::numeric_limits<Real>::epsilon()*abs(result))
        {
            ++consecutive_small_terms;
        }
        else
        {
          consecutive_small_terms = 0;
        }
    }
    return result;
}


}}

#endif
