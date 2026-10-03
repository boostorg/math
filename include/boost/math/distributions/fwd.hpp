// fwd.hpp Forward declarations of Boost.Math distributions.

// Copyright Paul A. Bristow 2007, 2010, 2012, 2014.
// Copyright John Maddock 2007.

// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_DISTRIBUTIONS_FWD_HPP
#define BOOST_MATH_DISTRIBUTIONS_FWD_HPP

#include <boost/math/tools/config.hpp>

BOOST_MATH_NAMESPACE_BEGIN

BOOST_MATH_EXPORT template <class RealType, class Policy>
class arcsine_distribution;

BOOST_MATH_EXPORT template <class RealType, class Policy>
class bernoulli_distribution;

BOOST_MATH_EXPORT template <class RealType, class Policy>
class beta_distribution;

BOOST_MATH_EXPORT template <class RealType, class Policy>
class binomial_distribution;

BOOST_MATH_EXPORT template <class RealType, class Policy>
class cauchy_distribution;

BOOST_MATH_EXPORT template <class RealType, class Policy>
class chi_squared_distribution;

BOOST_MATH_EXPORT template <class RealType, class Policy>
class exponential_distribution;

BOOST_MATH_EXPORT template <class RealType, class Policy>
class extreme_value_distribution;

BOOST_MATH_EXPORT template <class RealType, class Policy>
class fisher_f_distribution;

BOOST_MATH_EXPORT template <class RealType, class Policy>
class gamma_distribution;

BOOST_MATH_EXPORT template <class RealType, class Policy>
class geometric_distribution;

BOOST_MATH_EXPORT template <class RealType, class Policy>
class hyperexponential_distribution;

BOOST_MATH_EXPORT template <class RealType, class Policy>
class hypergeometric_distribution;

BOOST_MATH_EXPORT template <class RealType, class Policy>
class inverse_chi_squared_distribution;

BOOST_MATH_EXPORT template <class RealType, class Policy>
class inverse_gamma_distribution;

BOOST_MATH_EXPORT template <class RealType, class Policy>
class inverse_gaussian_distribution;

BOOST_MATH_EXPORT template <class RealType, class Policy>
class kolmogorov_smirnov_distribution;

BOOST_MATH_EXPORT template <class RealType, class Policy>
class landau_distribution;

BOOST_MATH_EXPORT template <class RealType, class Policy>
class mapairy_distribution;

BOOST_MATH_EXPORT template <class RealType, class Policy>
class holtsmark_distribution;

BOOST_MATH_EXPORT template <class RealType, class Policy>
class saspoint5_distribution;

BOOST_MATH_EXPORT template <class RealType, class Policy>
class laplace_distribution;

BOOST_MATH_EXPORT template <class RealType, class Policy>
class logistic_distribution;

BOOST_MATH_EXPORT template <class RealType, class Policy>
class lognormal_distribution;

BOOST_MATH_EXPORT template <class RealType, class Policy>
class negative_binomial_distribution;

BOOST_MATH_EXPORT template <class RealType, class Policy>
class non_central_beta_distribution;

BOOST_MATH_EXPORT template <class RealType, class Policy>
class non_central_chi_squared_distribution;

BOOST_MATH_EXPORT template <class RealType, class Policy>
class non_central_f_distribution;

BOOST_MATH_EXPORT template <class RealType, class Policy>
class non_central_t_distribution;

BOOST_MATH_EXPORT template <class RealType, class Policy>
class normal_distribution;

BOOST_MATH_EXPORT template <class RealType, class Policy>
class pareto_distribution;

BOOST_MATH_EXPORT template <class RealType, class Policy>
class poisson_distribution;

BOOST_MATH_EXPORT template <class RealType, class Policy>
class rayleigh_distribution;

BOOST_MATH_EXPORT template <class RealType, class Policy>
class skew_normal_distribution;

BOOST_MATH_EXPORT template <class RealType, class Policy>
class students_t_distribution;

BOOST_MATH_EXPORT template <class RealType, class Policy>
class triangular_distribution;

BOOST_MATH_EXPORT template <class RealType, class Policy>
class uniform_distribution;

BOOST_MATH_EXPORT template <class RealType, class Policy>
class weibull_distribution;

BOOST_MATH_EXPORT template <class RealType, class Policy>
class von_mises_distribution;

BOOST_MATH_NAMESPACE_END

#define BOOST_MATH_DECLARE_DISTRIBUTIONS(Type, Policy)\
   typedef BOOST_MATH_NAMESPACE::arcsine_distribution<Type, Policy> arcsine;\
   typedef BOOST_MATH_NAMESPACE::bernoulli_distribution<Type, Policy> bernoulli;\
   typedef BOOST_MATH_NAMESPACE::beta_distribution<Type, Policy> beta;\
   typedef BOOST_MATH_NAMESPACE::binomial_distribution<Type, Policy> binomial;\
   typedef BOOST_MATH_NAMESPACE::cauchy_distribution<Type, Policy> cauchy;\
   typedef BOOST_MATH_NAMESPACE::chi_squared_distribution<Type, Policy> chi_squared;\
   typedef BOOST_MATH_NAMESPACE::exponential_distribution<Type, Policy> exponential;\
   typedef BOOST_MATH_NAMESPACE::extreme_value_distribution<Type, Policy> extreme_value;\
   typedef BOOST_MATH_NAMESPACE::fisher_f_distribution<Type, Policy> fisher_f;\
   typedef BOOST_MATH_NAMESPACE::gamma_distribution<Type, Policy> gamma;\
   typedef BOOST_MATH_NAMESPACE::geometric_distribution<Type, Policy> geometric;\
   typedef BOOST_MATH_NAMESPACE::hypergeometric_distribution<Type, Policy> hypergeometric;\
   typedef BOOST_MATH_NAMESPACE::kolmogorov_smirnov_distribution<Type, Policy> kolmogorov_smirnov;\
   typedef BOOST_MATH_NAMESPACE::inverse_chi_squared_distribution<Type, Policy> inverse_chi_squared;\
   typedef BOOST_MATH_NAMESPACE::inverse_gaussian_distribution<Type, Policy> inverse_gaussian;\
   typedef BOOST_MATH_NAMESPACE::inverse_gamma_distribution<Type, Policy> inverse_gamma;\
   typedef BOOST_MATH_NAMESPACE::landau_distribution<Type, Policy> landau;\
   typedef BOOST_MATH_NAMESPACE::mapairy_distribution<Type, Policy> mapairy;\
   typedef BOOST_MATH_NAMESPACE::holtsmark_distribution<Type, Policy> holtsmark;\
   typedef BOOST_MATH_NAMESPACE::saspoint5_distribution<Type, Policy> saspoint5;\
   typedef BOOST_MATH_NAMESPACE::laplace_distribution<Type, Policy> laplace;\
   typedef BOOST_MATH_NAMESPACE::logistic_distribution<Type, Policy> logistic;\
   typedef BOOST_MATH_NAMESPACE::lognormal_distribution<Type, Policy> lognormal;\
   typedef BOOST_MATH_NAMESPACE::negative_binomial_distribution<Type, Policy> negative_binomial;\
   typedef BOOST_MATH_NAMESPACE::non_central_beta_distribution<Type, Policy> non_central_beta;\
   typedef BOOST_MATH_NAMESPACE::non_central_chi_squared_distribution<Type, Policy> non_central_chi_squared;\
   typedef BOOST_MATH_NAMESPACE::non_central_f_distribution<Type, Policy> non_central_f;\
   typedef BOOST_MATH_NAMESPACE::non_central_t_distribution<Type, Policy> non_central_t;\
   typedef BOOST_MATH_NAMESPACE::normal_distribution<Type, Policy> normal;\
   typedef BOOST_MATH_NAMESPACE::pareto_distribution<Type, Policy> pareto;\
   typedef BOOST_MATH_NAMESPACE::poisson_distribution<Type, Policy> poisson;\
   typedef BOOST_MATH_NAMESPACE::rayleigh_distribution<Type, Policy> rayleigh;\
   typedef BOOST_MATH_NAMESPACE::skew_normal_distribution<Type, Policy> skew_normal;\
   typedef BOOST_MATH_NAMESPACE::students_t_distribution<Type, Policy> students_t;\
   typedef BOOST_MATH_NAMESPACE::triangular_distribution<Type, Policy> triangular;\
   typedef BOOST_MATH_NAMESPACE::uniform_distribution<Type, Policy> uniform;\
   typedef BOOST_MATH_NAMESPACE::weibull_distribution<Type, Policy> weibull; \
   typedef BOOST_MATH_NAMESPACE::von_mises_distribution<Type, Policy> von_mises;

#endif // BOOST_MATH_DISTRIBUTIONS_FWD_HPP
