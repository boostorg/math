//  Copyright John Maddock 2006.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Basic sanity check that header <boost/math/special_functions/gamma.hpp>
// #includes all the files that it needs to.
//
#include <boost/math/special_functions/gamma.hpp>

template float boost::math::tgamma<float>(float);
template double boost::math::tgamma<double>(double);
template long double boost::math::tgamma<long double>(long double);

template float boost::math::lgamma<float>(float);
template double boost::math::lgamma<double>(double);
template long double boost::math::lgamma<long double>(long double);

template float boost::math::gamma_p<float>(float, float);
template double boost::math::gamma_p<double>(double, double);
template long double boost::math::gamma_p<long double>(long double, long double);

template float boost::math::gamma_q<float>(float, float);
template double boost::math::gamma_q<double>(double, double);
template long double boost::math::gamma_q<long double>(long double, long double);

template float boost::math::gamma_p_inv<float>(float, float);
template double boost::math::gamma_p_inv<double>(double, double);
template long double boost::math::gamma_p_inv<long double>(long double, long double);

template float boost::math::gamma_q_inv<float>(float, float);
template double boost::math::gamma_q_inv<double>(double, double);
template long double boost::math::gamma_q_inv<long double>(long double, long double);

template float boost::math::gamma_p_inva<float>(float, float);
template double boost::math::gamma_p_inva<double>(double, double);
template long double boost::math::gamma_p_inva<long double>(long double, long double);

template float boost::math::gamma_q_inva<float>(float, float);
template double boost::math::gamma_q_inva<double>(double, double);
template long double boost::math::gamma_q_inva<long double>(long double, long double);

template float boost::math::gamma_p_derivative<float>(float, float);
template double boost::math::gamma_p_derivative<double>(double, double);
template long double boost::math::gamma_p_derivative<long double>(long double, long double);

template float boost::math::tgamma_ratio<float>(float, float);
template double boost::math::tgamma_ratio<double>(double, double);
template long double boost::math::tgamma_ratio<long double>(long double, long double);

template float boost::math::tgamma_delta_ratio<float>(float, float);
template double boost::math::tgamma_delta_ratio<double>(double, double);
template long double boost::math::tgamma_delta_ratio<long double>(long double, long double);

