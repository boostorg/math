
// Copyright John Maddock 2007.

// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <boost/math/policy/policy.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/test/included/test_exec_monitor.hpp> // for test_main

template <class P1, class P2>
bool check_same(const P1&, const P2&)
{
   if(!boost::is_same<P1, P2>::value)
   {
      std::cout << "P1 = " << typeid(P1).name() << std::endl;
      std::cout << "P2 = " << typeid(P2).name() << std::endl;
   }
   return boost::is_same<P1, P2>::value;
}


int test_main(int, char* [])
{
   using namespace boost::math::policies;
   using namespace boost;
   BOOST_CHECK(is_domain_error<domain_error<ignore_error> >::value);
   BOOST_CHECK(0 == is_domain_error<pole_error<ignore_error> >::value);
   BOOST_CHECK(is_pole_error<pole_error<ignore_error> >::value);
   BOOST_CHECK(0 == is_pole_error<domain_error<ignore_error> >::value);
   BOOST_CHECK(is_digits10<digits10<ignore_error> >::value);
   BOOST_CHECK(0 == is_digits10<digits2<ignore_error> >::value);

   BOOST_CHECK((is_same<policy<>::domain_error_type, domain_error<BOOST_MATH_DOMAIN_ERROR_POLICY> >::value));
   BOOST_CHECK((is_same<policy<>::evaluation_error_type, evaluation_error<BOOST_MATH_EVALUATION_ERROR_POLICY> >::value));
	
   BOOST_CHECK((is_same<policy<domain_error<ignore_error> >::domain_error_type, domain_error<ignore_error> >::value));
   BOOST_CHECK((is_same<policy<domain_error<ignore_error> >::pole_error_type, pole_error<BOOST_MATH_POLE_ERROR_POLICY> >::value));
   BOOST_CHECK((is_same<policy<domain_error<ignore_error> >::overflow_error_type, overflow_error<BOOST_MATH_OVERFLOW_ERROR_POLICY> >::value));
   BOOST_CHECK((is_same<policy<domain_error<ignore_error> >::underflow_error_type, underflow_error<BOOST_MATH_UNDERFLOW_ERROR_POLICY> >::value));
   BOOST_CHECK((is_same<policy<domain_error<ignore_error> >::denorm_error_type, denorm_error<BOOST_MATH_DENORM_ERROR_POLICY> >::value));
   BOOST_CHECK((is_same<policy<domain_error<ignore_error> >::evaluation_error_type, evaluation_error<BOOST_MATH_EVALUATION_ERROR_POLICY> >::value));
   BOOST_CHECK((is_same<policy<domain_error<ignore_error> >::precision_type, policy<>::precision_type>::value));
   BOOST_CHECK((is_same<policy<domain_error<ignore_error> >::promote_float_type, policy<>::promote_float_type>::value));
   BOOST_CHECK((is_same<policy<domain_error<ignore_error> >::promote_double_type, policy<>::promote_double_type>::value));
   BOOST_CHECK((is_same<policy<domain_error<ignore_error> >::discrete_quantile_type, policy<>::discrete_quantile_type>::value));
   
   BOOST_CHECK((is_same<policy<pole_error<user_error> >::domain_error_type, policy<>::domain_error_type >::value));
   BOOST_CHECK((is_same<policy<pole_error<user_error> >::pole_error_type, pole_error<user_error> >::value));
   BOOST_CHECK((is_same<policy<pole_error<user_error> >::overflow_error_type, overflow_error<BOOST_MATH_OVERFLOW_ERROR_POLICY> >::value));
   BOOST_CHECK((is_same<policy<pole_error<user_error> >::underflow_error_type, underflow_error<BOOST_MATH_UNDERFLOW_ERROR_POLICY> >::value));
   BOOST_CHECK((is_same<policy<pole_error<user_error> >::denorm_error_type, denorm_error<BOOST_MATH_DENORM_ERROR_POLICY> >::value));
   BOOST_CHECK((is_same<policy<pole_error<user_error> >::evaluation_error_type, evaluation_error<BOOST_MATH_EVALUATION_ERROR_POLICY> >::value));
   BOOST_CHECK((is_same<policy<pole_error<user_error> >::precision_type, policy<>::precision_type>::value));
   BOOST_CHECK((is_same<policy<pole_error<user_error> >::promote_float_type, policy<>::promote_float_type>::value));
   BOOST_CHECK((is_same<policy<pole_error<user_error> >::promote_double_type, policy<>::promote_double_type>::value));
   BOOST_CHECK((is_same<policy<pole_error<user_error> >::discrete_quantile_type, policy<>::discrete_quantile_type>::value));
   
   BOOST_CHECK((is_same<policy<overflow_error<errno_on_error> >::domain_error_type, policy<>::domain_error_type >::value));
   BOOST_CHECK((is_same<policy<overflow_error<errno_on_error> >::pole_error_type, policy<>::pole_error_type >::value));
   BOOST_CHECK((is_same<policy<overflow_error<errno_on_error> >::overflow_error_type, overflow_error<errno_on_error> >::value));
   BOOST_CHECK((is_same<policy<overflow_error<errno_on_error> >::underflow_error_type, underflow_error<BOOST_MATH_UNDERFLOW_ERROR_POLICY> >::value));
   BOOST_CHECK((is_same<policy<overflow_error<errno_on_error> >::denorm_error_type, denorm_error<BOOST_MATH_DENORM_ERROR_POLICY> >::value));
   BOOST_CHECK((is_same<policy<overflow_error<errno_on_error> >::evaluation_error_type, evaluation_error<BOOST_MATH_EVALUATION_ERROR_POLICY> >::value));
   BOOST_CHECK((is_same<policy<overflow_error<errno_on_error> >::precision_type, policy<>::precision_type>::value));
   BOOST_CHECK((is_same<policy<overflow_error<errno_on_error> >::promote_float_type, policy<>::promote_float_type>::value));
   BOOST_CHECK((is_same<policy<overflow_error<errno_on_error> >::promote_double_type, policy<>::promote_double_type>::value));
   BOOST_CHECK((is_same<policy<overflow_error<errno_on_error> >::discrete_quantile_type, policy<>::discrete_quantile_type>::value));
   
   BOOST_CHECK((is_same<policy<underflow_error<errno_on_error> >::domain_error_type, policy<>::domain_error_type >::value));
   BOOST_CHECK((is_same<policy<underflow_error<errno_on_error> >::pole_error_type, policy<>::pole_error_type >::value));
   BOOST_CHECK((is_same<policy<underflow_error<errno_on_error> >::overflow_error_type, policy<>::overflow_error_type >::value));
   BOOST_CHECK((is_same<policy<underflow_error<errno_on_error> >::underflow_error_type, underflow_error<errno_on_error> >::value));
   BOOST_CHECK((is_same<policy<underflow_error<errno_on_error> >::denorm_error_type, denorm_error<BOOST_MATH_DENORM_ERROR_POLICY> >::value));
   BOOST_CHECK((is_same<policy<underflow_error<errno_on_error> >::evaluation_error_type, evaluation_error<BOOST_MATH_EVALUATION_ERROR_POLICY> >::value));
   BOOST_CHECK((is_same<policy<underflow_error<errno_on_error> >::precision_type, policy<>::precision_type>::value));
   BOOST_CHECK((is_same<policy<underflow_error<errno_on_error> >::promote_float_type, policy<>::promote_float_type>::value));
   BOOST_CHECK((is_same<policy<underflow_error<errno_on_error> >::promote_double_type, policy<>::promote_double_type>::value));
   BOOST_CHECK((is_same<policy<underflow_error<errno_on_error> >::discrete_quantile_type, policy<>::discrete_quantile_type>::value));
   
   BOOST_CHECK((is_same<policy<denorm_error<errno_on_error> >::domain_error_type, policy<>::domain_error_type >::value));
   BOOST_CHECK((is_same<policy<denorm_error<errno_on_error> >::pole_error_type, policy<>::pole_error_type >::value));
   BOOST_CHECK((is_same<policy<denorm_error<errno_on_error> >::overflow_error_type, policy<>::overflow_error_type >::value));
   BOOST_CHECK((is_same<policy<denorm_error<errno_on_error> >::underflow_error_type, policy<>::underflow_error_type >::value));
   BOOST_CHECK((is_same<policy<denorm_error<errno_on_error> >::denorm_error_type, denorm_error<errno_on_error> >::value));
   BOOST_CHECK((is_same<policy<denorm_error<errno_on_error> >::evaluation_error_type, evaluation_error<BOOST_MATH_EVALUATION_ERROR_POLICY> >::value));
   BOOST_CHECK((is_same<policy<denorm_error<errno_on_error> >::precision_type, policy<>::precision_type>::value));
   BOOST_CHECK((is_same<policy<denorm_error<errno_on_error> >::promote_float_type, policy<>::promote_float_type>::value));
   BOOST_CHECK((is_same<policy<denorm_error<errno_on_error> >::promote_double_type, policy<>::promote_double_type>::value));
   BOOST_CHECK((is_same<policy<denorm_error<errno_on_error> >::discrete_quantile_type, policy<>::discrete_quantile_type>::value));
   
   BOOST_CHECK((is_same<policy<evaluation_error<errno_on_error> >::domain_error_type, policy<>::domain_error_type >::value));
   BOOST_CHECK((is_same<policy<evaluation_error<errno_on_error> >::pole_error_type, policy<>::pole_error_type >::value));
   BOOST_CHECK((is_same<policy<evaluation_error<errno_on_error> >::overflow_error_type, policy<>::overflow_error_type >::value));
   BOOST_CHECK((is_same<policy<evaluation_error<errno_on_error> >::underflow_error_type, policy<>::underflow_error_type >::value));
   BOOST_CHECK((is_same<policy<evaluation_error<errno_on_error> >::denorm_error_type, policy<>::denorm_error_type >::value));
   BOOST_CHECK((is_same<policy<evaluation_error<errno_on_error> >::evaluation_error_type, evaluation_error<errno_on_error> >::value));
   BOOST_CHECK((is_same<policy<evaluation_error<errno_on_error> >::precision_type, policy<>::precision_type>::value));
   BOOST_CHECK((is_same<policy<evaluation_error<errno_on_error> >::promote_float_type, policy<>::promote_float_type>::value));
   BOOST_CHECK((is_same<policy<evaluation_error<errno_on_error> >::promote_double_type, policy<>::promote_double_type>::value));
   BOOST_CHECK((is_same<policy<evaluation_error<errno_on_error> >::discrete_quantile_type, policy<>::discrete_quantile_type>::value));
   
   BOOST_CHECK((is_same<policy<digits2<20> >::domain_error_type, policy<>::domain_error_type >::value));
   BOOST_CHECK((is_same<policy<digits2<20> >::pole_error_type, policy<>::pole_error_type >::value));
   BOOST_CHECK((is_same<policy<digits2<20> >::overflow_error_type, policy<>::overflow_error_type >::value));
   BOOST_CHECK((is_same<policy<digits2<20> >::underflow_error_type, policy<>::underflow_error_type >::value));
   BOOST_CHECK((is_same<policy<digits2<20> >::denorm_error_type, policy<>::denorm_error_type >::value));
   BOOST_CHECK((is_same<policy<digits2<20> >::evaluation_error_type, policy<>::evaluation_error_type >::value));
   BOOST_CHECK((is_same<policy<digits2<20> >::precision_type, digits2<20> >::value));
   BOOST_CHECK((is_same<policy<digits2<20> >::promote_float_type, policy<>::promote_float_type>::value));
   BOOST_CHECK((is_same<policy<digits2<20> >::promote_double_type, policy<>::promote_double_type>::value));
   BOOST_CHECK((is_same<policy<digits2<20> >::discrete_quantile_type, policy<>::discrete_quantile_type>::value));
   
   BOOST_CHECK((is_same<policy<promote_float<false> >::domain_error_type, policy<>::domain_error_type >::value));
   BOOST_CHECK((is_same<policy<promote_float<false> >::pole_error_type, policy<>::pole_error_type >::value));
   BOOST_CHECK((is_same<policy<promote_float<false> >::overflow_error_type, policy<>::overflow_error_type >::value));
   BOOST_CHECK((is_same<policy<promote_float<false> >::underflow_error_type, policy<>::underflow_error_type >::value));
   BOOST_CHECK((is_same<policy<promote_float<false> >::denorm_error_type, policy<>::denorm_error_type >::value));
   BOOST_CHECK((is_same<policy<promote_float<false> >::evaluation_error_type, policy<>::evaluation_error_type >::value));
   BOOST_CHECK((is_same<policy<promote_float<false> >::precision_type, policy<>::precision_type >::value));
   BOOST_CHECK((is_same<policy<promote_float<false> >::promote_float_type, promote_float<false> >::value));
   BOOST_CHECK((is_same<policy<promote_float<false> >::promote_double_type, policy<>::promote_double_type>::value));
   BOOST_CHECK((is_same<policy<promote_float<false> >::discrete_quantile_type, policy<>::discrete_quantile_type>::value));
   
   BOOST_CHECK((is_same<policy<promote_double<false> >::domain_error_type, policy<>::domain_error_type >::value));
   BOOST_CHECK((is_same<policy<promote_double<false> >::pole_error_type, policy<>::pole_error_type >::value));
   BOOST_CHECK((is_same<policy<promote_double<false> >::overflow_error_type, policy<>::overflow_error_type >::value));
   BOOST_CHECK((is_same<policy<promote_double<false> >::underflow_error_type, policy<>::underflow_error_type >::value));
   BOOST_CHECK((is_same<policy<promote_double<false> >::denorm_error_type, policy<>::denorm_error_type >::value));
   BOOST_CHECK((is_same<policy<promote_double<false> >::evaluation_error_type, policy<>::evaluation_error_type >::value));
   BOOST_CHECK((is_same<policy<promote_double<false> >::precision_type, policy<>::precision_type >::value));
   BOOST_CHECK((is_same<policy<promote_double<false> >::promote_float_type,  policy<>::promote_float_type>::value));
   BOOST_CHECK((is_same<policy<promote_double<false> >::promote_double_type, promote_double<false> >::value));
   BOOST_CHECK((is_same<policy<promote_double<false> >::discrete_quantile_type, policy<>::discrete_quantile_type>::value));
   
   BOOST_CHECK((is_same<policy<discrete_quantile<integer_above> >::domain_error_type, policy<>::domain_error_type >::value));
   BOOST_CHECK((is_same<policy<discrete_quantile<integer_above> >::pole_error_type, policy<>::pole_error_type >::value));
   BOOST_CHECK((is_same<policy<discrete_quantile<integer_above> >::overflow_error_type, policy<>::overflow_error_type >::value));
   BOOST_CHECK((is_same<policy<discrete_quantile<integer_above> >::underflow_error_type, policy<>::underflow_error_type >::value));
   BOOST_CHECK((is_same<policy<discrete_quantile<integer_above> >::denorm_error_type, policy<>::denorm_error_type >::value));
   BOOST_CHECK((is_same<policy<discrete_quantile<integer_above> >::evaluation_error_type, policy<>::evaluation_error_type >::value));
   BOOST_CHECK((is_same<policy<discrete_quantile<integer_above> >::precision_type, policy<>::precision_type >::value));
   BOOST_CHECK((is_same<policy<discrete_quantile<integer_above> >::promote_float_type,  policy<>::promote_float_type>::value));
   BOOST_CHECK((is_same<policy<discrete_quantile<integer_above> >::promote_double_type, policy<>::promote_double_type>::value));
   BOOST_CHECK((is_same<policy<discrete_quantile<integer_above> >::discrete_quantile_type, discrete_quantile<integer_above> >::value));

   //
   // Now try again with 2 policies:
   //
   BOOST_CHECK((is_same<policy<domain_error<ignore_error>, overflow_error<ignore_error> >::domain_error_type, domain_error<ignore_error> >::value));
   BOOST_CHECK((is_same<policy<domain_error<ignore_error>, overflow_error<ignore_error> >::pole_error_type, pole_error<BOOST_MATH_POLE_ERROR_POLICY> >::value));
   BOOST_CHECK((is_same<policy<domain_error<ignore_error>, overflow_error<ignore_error> >::overflow_error_type, overflow_error<ignore_error> >::value));
   BOOST_CHECK((is_same<policy<domain_error<ignore_error>, overflow_error<ignore_error> >::underflow_error_type, underflow_error<BOOST_MATH_UNDERFLOW_ERROR_POLICY> >::value));
   BOOST_CHECK((is_same<policy<domain_error<ignore_error>, overflow_error<ignore_error> >::denorm_error_type, denorm_error<BOOST_MATH_DENORM_ERROR_POLICY> >::value));
   BOOST_CHECK((is_same<policy<domain_error<ignore_error>, overflow_error<ignore_error> >::evaluation_error_type, evaluation_error<BOOST_MATH_EVALUATION_ERROR_POLICY> >::value));
   BOOST_CHECK((is_same<policy<domain_error<ignore_error>, overflow_error<ignore_error> >::precision_type, policy<>::precision_type>::value));
   BOOST_CHECK((is_same<policy<domain_error<ignore_error>, overflow_error<ignore_error> >::promote_float_type, policy<>::promote_float_type>::value));
   BOOST_CHECK((is_same<policy<domain_error<ignore_error>, overflow_error<ignore_error> >::promote_double_type, policy<>::promote_double_type>::value));
   BOOST_CHECK((is_same<policy<domain_error<ignore_error>, overflow_error<ignore_error> >::discrete_quantile_type, policy<>::discrete_quantile_type>::value));
   
   BOOST_CHECK((is_same<policy<denorm_error<throw_on_error>, domain_error<ignore_error> >::domain_error_type, domain_error<ignore_error> >::value));
   BOOST_CHECK((is_same<policy<denorm_error<throw_on_error>, domain_error<ignore_error> >::pole_error_type, pole_error<BOOST_MATH_POLE_ERROR_POLICY> >::value));
   BOOST_CHECK((is_same<policy<denorm_error<throw_on_error>, domain_error<ignore_error> >::overflow_error_type, overflow_error<BOOST_MATH_OVERFLOW_ERROR_POLICY> >::value));
   BOOST_CHECK((is_same<policy<denorm_error<throw_on_error>, domain_error<ignore_error> >::underflow_error_type, underflow_error<BOOST_MATH_UNDERFLOW_ERROR_POLICY> >::value));
   BOOST_CHECK((is_same<policy<denorm_error<throw_on_error>, domain_error<ignore_error> >::denorm_error_type, denorm_error<throw_on_error> >::value));
   BOOST_CHECK((is_same<policy<denorm_error<throw_on_error>, domain_error<ignore_error> >::evaluation_error_type, evaluation_error<BOOST_MATH_EVALUATION_ERROR_POLICY> >::value));
   BOOST_CHECK((is_same<policy<denorm_error<throw_on_error>, domain_error<ignore_error> >::precision_type, policy<>::precision_type>::value));
   BOOST_CHECK((is_same<policy<denorm_error<throw_on_error>, domain_error<ignore_error> >::promote_float_type, policy<>::promote_float_type>::value));
   BOOST_CHECK((is_same<policy<denorm_error<throw_on_error>, domain_error<ignore_error> >::promote_double_type, policy<>::promote_double_type>::value));
   BOOST_CHECK((is_same<policy<denorm_error<throw_on_error>, domain_error<ignore_error> >::discrete_quantile_type, policy<>::discrete_quantile_type>::value));
   
   BOOST_CHECK((is_same<policy<digits2<20>, promote_float<false>, discrete_quantile<integer_below>, denorm_error<throw_on_error>, domain_error<ignore_error> >::domain_error_type, domain_error<ignore_error> >::value));
   BOOST_CHECK((is_same<policy<digits2<20>, promote_float<false>, discrete_quantile<integer_below>, denorm_error<throw_on_error>, domain_error<ignore_error> >::pole_error_type, pole_error<BOOST_MATH_POLE_ERROR_POLICY> >::value));
   BOOST_CHECK((is_same<policy<digits2<20>, promote_float<false>, discrete_quantile<integer_below>, denorm_error<throw_on_error>, domain_error<ignore_error> >::overflow_error_type, overflow_error<BOOST_MATH_OVERFLOW_ERROR_POLICY> >::value));
   BOOST_CHECK((is_same<policy<digits2<20>, promote_float<false>, discrete_quantile<integer_below>, denorm_error<throw_on_error>, domain_error<ignore_error> >::underflow_error_type, underflow_error<BOOST_MATH_UNDERFLOW_ERROR_POLICY> >::value));
   BOOST_CHECK((is_same<policy<digits2<20>, promote_float<false>, discrete_quantile<integer_below>, denorm_error<throw_on_error>, domain_error<ignore_error> >::denorm_error_type, denorm_error<throw_on_error> >::value));
   BOOST_CHECK((is_same<policy<digits2<20>, promote_float<false>, discrete_quantile<integer_below>, denorm_error<throw_on_error>, domain_error<ignore_error> >::evaluation_error_type, evaluation_error<BOOST_MATH_EVALUATION_ERROR_POLICY> >::value));
   BOOST_CHECK((is_same<policy<digits2<20>, promote_float<false>, discrete_quantile<integer_below>, denorm_error<throw_on_error>, domain_error<ignore_error> >::precision_type, digits2<20> >::value));
   BOOST_CHECK((is_same<policy<digits2<20>, promote_float<false>, discrete_quantile<integer_below>, denorm_error<throw_on_error>, domain_error<ignore_error> >::promote_float_type, promote_float<false> >::value));
   BOOST_CHECK((is_same<policy<digits2<20>, promote_float<false>, discrete_quantile<integer_below>, denorm_error<throw_on_error>, domain_error<ignore_error> >::promote_double_type, policy<>::promote_double_type>::value));
   BOOST_CHECK((is_same<policy<digits2<20>, promote_float<false>, discrete_quantile<integer_below>, denorm_error<throw_on_error>, domain_error<ignore_error> >::discrete_quantile_type, discrete_quantile<integer_below> >::value));
   
   BOOST_CHECK((is_same<normalise<policy<digits2<20>, promote_float<false>, discrete_quantile<integer_below>, denorm_error<throw_on_error>, domain_error<ignore_error> > >::type::domain_error_type, domain_error<ignore_error> >::value));
   BOOST_CHECK((is_same<normalise<policy<digits2<20>, promote_float<false>, discrete_quantile<integer_below>, denorm_error<throw_on_error>, domain_error<ignore_error> > >::type::pole_error_type, pole_error<BOOST_MATH_POLE_ERROR_POLICY> >::value));
   BOOST_CHECK((is_same<normalise<policy<digits2<20>, promote_float<false>, discrete_quantile<integer_below>, denorm_error<throw_on_error>, domain_error<ignore_error> > >::type::overflow_error_type, overflow_error<BOOST_MATH_OVERFLOW_ERROR_POLICY> >::value));
   BOOST_CHECK((is_same<normalise<policy<digits2<20>, promote_float<false>, discrete_quantile<integer_below>, denorm_error<throw_on_error>, domain_error<ignore_error> > >::type::underflow_error_type, underflow_error<BOOST_MATH_UNDERFLOW_ERROR_POLICY> >::value));
   BOOST_CHECK((is_same<normalise<policy<digits2<20>, promote_float<false>, discrete_quantile<integer_below>, denorm_error<throw_on_error>, domain_error<ignore_error> > >::type::denorm_error_type, denorm_error<throw_on_error> >::value));
   BOOST_CHECK((is_same<normalise<policy<digits2<20>, promote_float<false>, discrete_quantile<integer_below>, denorm_error<throw_on_error>, domain_error<ignore_error> > >::type::evaluation_error_type, evaluation_error<BOOST_MATH_EVALUATION_ERROR_POLICY> >::value));
   BOOST_CHECK((is_same<normalise<policy<digits2<20>, promote_float<false>, discrete_quantile<integer_below>, denorm_error<throw_on_error>, domain_error<ignore_error> > >::type::precision_type, digits2<20> >::value));
   BOOST_CHECK((is_same<normalise<policy<digits2<20>, promote_float<false>, discrete_quantile<integer_below>, denorm_error<throw_on_error>, domain_error<ignore_error> > >::type::promote_float_type, promote_float<false> >::value));
   BOOST_CHECK((is_same<normalise<policy<digits2<20>, promote_float<false>, discrete_quantile<integer_below>, denorm_error<throw_on_error>, domain_error<ignore_error> > >::type::promote_double_type, policy<>::promote_double_type>::value));
   BOOST_CHECK((is_same<normalise<policy<digits2<20>, promote_float<false>, discrete_quantile<integer_below>, denorm_error<throw_on_error>, domain_error<ignore_error> > >::type::discrete_quantile_type, discrete_quantile<integer_below> >::value));
   
   BOOST_CHECK((is_same<normalise<policy<>, digits2<20>, promote_float<false>, discrete_quantile<integer_below>, denorm_error<throw_on_error>, domain_error<ignore_error> >::type::domain_error_type, domain_error<ignore_error> >::value));
   BOOST_CHECK((is_same<normalise<policy<>, digits2<20>, promote_float<false>, discrete_quantile<integer_below>, denorm_error<throw_on_error>, domain_error<ignore_error> >::type::pole_error_type, pole_error<BOOST_MATH_POLE_ERROR_POLICY> >::value));
   BOOST_CHECK((is_same<normalise<policy<>, digits2<20>, promote_float<false>, discrete_quantile<integer_below>, denorm_error<throw_on_error>, domain_error<ignore_error>  >::type::overflow_error_type, overflow_error<BOOST_MATH_OVERFLOW_ERROR_POLICY> >::value));
   BOOST_CHECK((is_same<normalise<policy<>, digits2<20>, promote_float<false>, discrete_quantile<integer_below>, denorm_error<throw_on_error>, domain_error<ignore_error>  >::type::underflow_error_type, underflow_error<BOOST_MATH_UNDERFLOW_ERROR_POLICY> >::value));
   BOOST_CHECK((is_same<normalise<policy<>, digits2<20>, promote_float<false>, discrete_quantile<integer_below>, denorm_error<throw_on_error>, domain_error<ignore_error>  >::type::denorm_error_type, denorm_error<throw_on_error> >::value));
   BOOST_CHECK((is_same<normalise<policy<>, digits2<20>, promote_float<false>, discrete_quantile<integer_below>, denorm_error<throw_on_error>, domain_error<ignore_error>  >::type::evaluation_error_type, evaluation_error<BOOST_MATH_EVALUATION_ERROR_POLICY> >::value));
   BOOST_CHECK((is_same<normalise<policy<>, digits2<20>, promote_float<false>, discrete_quantile<integer_below>, denorm_error<throw_on_error>, domain_error<ignore_error>  >::type::precision_type, digits2<20> >::value));
   BOOST_CHECK((is_same<normalise<policy<>, digits2<20>, promote_float<false>, discrete_quantile<integer_below>, denorm_error<throw_on_error>, domain_error<ignore_error>  >::type::promote_float_type, promote_float<false> >::value));
   BOOST_CHECK((is_same<normalise<policy<>, digits2<20>, promote_float<false>, discrete_quantile<integer_below>, denorm_error<throw_on_error>, domain_error<ignore_error>  >::type::promote_double_type, policy<>::promote_double_type>::value));
   BOOST_CHECK((is_same<normalise<policy<>, digits2<20>, promote_float<false>, discrete_quantile<integer_below>, denorm_error<throw_on_error>, domain_error<ignore_error>  >::type::discrete_quantile_type, discrete_quantile<integer_below> >::value));

   BOOST_CHECK(check_same(make_policy(), policy<>()));
   BOOST_CHECK(check_same(make_policy(denorm_error<ignore_error>()), normalise<policy<denorm_error<ignore_error> > >::type()));
   BOOST_CHECK(check_same(make_policy(digits2<20>()), normalise<policy<digits2<20> > >::type()));
   BOOST_CHECK(check_same(make_policy(promote_float<false>()), normalise<policy<promote_float<false> > >::type()));
   BOOST_CHECK(check_same(make_policy(domain_error<ignore_error>()), normalise<policy<domain_error<ignore_error> > >::type()));
   BOOST_CHECK(check_same(make_policy(pole_error<ignore_error>()), normalise<policy<pole_error<ignore_error> > >::type()));

   BOOST_CHECK(check_same(make_policy(domain_error<ignore_error>()), policy<domain_error<ignore_error> >()));
   BOOST_CHECK(check_same(make_policy(domain_error<ignore_error>(), pole_error<ignore_error>()), policy<domain_error<ignore_error>, pole_error<ignore_error> >()));
   BOOST_CHECK(check_same(make_policy(domain_error<ignore_error>(), pole_error<ignore_error>(), overflow_error<ignore_error>()), policy<domain_error<ignore_error>, pole_error<ignore_error>, overflow_error<ignore_error> >()));
   BOOST_CHECK(check_same(make_policy(domain_error<ignore_error>(), pole_error<ignore_error>(), overflow_error<ignore_error>(), underflow_error<throw_on_error>()), policy<domain_error<ignore_error>, pole_error<ignore_error>, overflow_error<ignore_error>, underflow_error<throw_on_error> >()));
   BOOST_CHECK(check_same(make_policy(domain_error<ignore_error>(), pole_error<ignore_error>(), overflow_error<ignore_error>(), underflow_error<throw_on_error>(), denorm_error<throw_on_error>()), policy<domain_error<ignore_error>, pole_error<ignore_error>, overflow_error<ignore_error>, underflow_error<throw_on_error>, denorm_error<throw_on_error> >()));
   BOOST_CHECK(check_same(make_policy(domain_error<ignore_error>(), pole_error<ignore_error>(), overflow_error<ignore_error>(), underflow_error<throw_on_error>(), denorm_error<throw_on_error>(), evaluation_error<ignore_error>()), policy<domain_error<ignore_error>, pole_error<ignore_error>, overflow_error<ignore_error>, underflow_error<throw_on_error>, denorm_error<throw_on_error>, evaluation_error<ignore_error> >()));
   BOOST_CHECK(check_same(make_policy(domain_error<ignore_error>(), pole_error<ignore_error>(), overflow_error<ignore_error>(), underflow_error<throw_on_error>(), denorm_error<throw_on_error>(), evaluation_error<ignore_error>(), digits2<10>()), policy<domain_error<ignore_error>, pole_error<ignore_error>, overflow_error<ignore_error>, underflow_error<throw_on_error>, denorm_error<throw_on_error>, evaluation_error<ignore_error>, digits2<10> >()));
   BOOST_CHECK(check_same(make_policy(domain_error<ignore_error>(), pole_error<ignore_error>(), overflow_error<ignore_error>(), underflow_error<throw_on_error>(), denorm_error<throw_on_error>(), evaluation_error<ignore_error>(), digits10<5>()), policy<domain_error<ignore_error>, pole_error<ignore_error>, overflow_error<ignore_error>, underflow_error<throw_on_error>, denorm_error<throw_on_error>, evaluation_error<ignore_error>, digits2<19> >()));
   BOOST_CHECK(check_same(make_policy(domain_error<ignore_error>(), pole_error<ignore_error>(), overflow_error<ignore_error>(), underflow_error<throw_on_error>(), denorm_error<throw_on_error>(), evaluation_error<ignore_error>(), digits2<10>(), promote_float<false>()), policy<domain_error<ignore_error>, pole_error<ignore_error>, overflow_error<ignore_error>, underflow_error<throw_on_error>, denorm_error<throw_on_error>, evaluation_error<ignore_error>, digits2<10>, promote_float<false> >()));
   BOOST_CHECK(check_same(make_policy(domain_error<ignore_error>(), pole_error<ignore_error>(), overflow_error<ignore_error>(), underflow_error<throw_on_error>(), denorm_error<throw_on_error>(), evaluation_error<ignore_error>(), digits2<10>(), promote_float<false>(), promote_double<false>()), policy<domain_error<ignore_error>, pole_error<ignore_error>, overflow_error<ignore_error>, underflow_error<throw_on_error>, denorm_error<throw_on_error>, evaluation_error<ignore_error>, digits2<10>, promote_float<false>, promote_double<false> >()));
   BOOST_CHECK(check_same(make_policy(domain_error<ignore_error>(), pole_error<ignore_error>(), overflow_error<ignore_error>(), underflow_error<throw_on_error>(), denorm_error<throw_on_error>(), evaluation_error<ignore_error>(), digits2<10>(), promote_float<false>(), promote_double<false>(), discrete_quantile<integer_below>()), policy<domain_error<ignore_error>, pole_error<ignore_error>, overflow_error<ignore_error>, underflow_error<throw_on_error>, denorm_error<throw_on_error>, evaluation_error<ignore_error>, digits2<10>, promote_float<false>, promote_double<false>, discrete_quantile<integer_below> >()));

   return 0;
} // int test_main(int, char* [])



