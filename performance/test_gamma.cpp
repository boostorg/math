
#include "performance_measure.hpp"

#include <boost/math/special_functions/gamma.hpp>
#include <boost/array.hpp>

#define T double
#include "../test/test_gamma_data.ipp"

template <std::size_t N>
double gamma_evaluate2(const boost::array<boost::array<T, 3>, N>& data)
{
   double result = 0;
   for(unsigned i = 0; i < N; ++i)
      result += boost::math::tgamma(data[i][0]);
   return result;
}

BOOST_MATH_PERFORMANCE_TEST(gamma_test, "gamma")
{
   double result = gamma_evaluate2(factorials);
   result += gamma_evaluate2(near_1);
   result += gamma_evaluate2(near_2);
   result += gamma_evaluate2(near_0);
   result += gamma_evaluate2(near_m10);
   result += gamma_evaluate2(near_m55);

   consume_result(result);
   set_call_count(
      (sizeof(factorials) 
      + sizeof(near_1) 
      + sizeof(near_2)
      + sizeof(near_0)
      + sizeof(near_m10)
      + sizeof(near_m55)) / sizeof(factorials[0]));
}

template <std::size_t N>
double lgamma_evaluate2(const boost::array<boost::array<T, 3>, N>& data)
{
   double result = 0;
   for(unsigned i = 0; i < N; ++i)
      result += boost::math::lgamma(data[i][0]);
   return result;
}

BOOST_MATH_PERFORMANCE_TEST(lgamma_test, "lgamma")
{
   double result = lgamma_evaluate2(factorials);
   result += lgamma_evaluate2(near_1);
   result += lgamma_evaluate2(near_2);
   result += lgamma_evaluate2(near_0);
   result += lgamma_evaluate2(near_m10);
   result += lgamma_evaluate2(near_m55);

   consume_result(result);
   set_call_count(
      (sizeof(factorials) 
      + sizeof(near_1) 
      + sizeof(near_2)
      + sizeof(near_0)
      + sizeof(near_m10)
      + sizeof(near_m55)) / sizeof(factorials[0]));
}
