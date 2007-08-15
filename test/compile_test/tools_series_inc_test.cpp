//  Copyright John Maddock 2006.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Basic sanity check that header <boost/math/tools/series.hpp>
// #includes all the files that it needs to.
//
#include <boost/math/tools/series.hpp>

struct Functor
{
   typedef double result_type;
   double operator()();
};
#define U double

template Functor::result_type boost::math::tools::sum_series<Functor>(Functor& func, int bits);
template Functor::result_type boost::math::tools::sum_series<Functor>(Functor& func, int bits, boost::uintmax_t& max_terms);
template Functor::result_type boost::math::tools::sum_series<Functor, U>(Functor& func, int bits, U init_value);
template Functor::result_type boost::math::tools::sum_series<Functor, U>(Functor& func, int bits, boost::uintmax_t& max_terms, U init_value);
template Functor::result_type boost::math::tools::kahan_sum_series<Functor>(Functor& func, int bits);
template Functor::result_type boost::math::tools::kahan_sum_series<Functor>(Functor& func, int bits, boost::uintmax_t& max_terms);


