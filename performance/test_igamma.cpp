
#include "performance_measure.hpp"

#include <boost/math/special_functions/gamma.hpp>
#include <boost/array.hpp>

#define T double
#include "../test/igamma_big_data.ipp"
#include "../test/igamma_int_data.ipp"
#include "../test/igamma_med_data.ipp"
#include "../test/igamma_small_data.ipp"

template <std::size_t N>
double igamma_evaluate2(const boost::array<boost::array<T, 6>, N>& data)
{
   double result = 0;
   for(unsigned i = 0; i < N; ++i)
      result += boost::math::gamma_p(data[i][0], data[i][1]);
   return result;
}

BOOST_MATH_PERFORMANCE_TEST(igamma_test, "igamma")
{
   double result = igamma_evaluate2(igamma_big_data);
   result += igamma_evaluate2(igamma_int_data);
   result += igamma_evaluate2(igamma_med_data);
   result += igamma_evaluate2(igamma_small_data);

   consume_result(result);
   set_call_count(
      (sizeof(igamma_big_data) 
      + sizeof(igamma_int_data) 
      + sizeof(igamma_med_data)
      + sizeof(igamma_small_data)) / sizeof(igamma_big_data[0]));
}

template <std::size_t N>
double igamma_inv_evaluate2(const boost::array<boost::array<T, 6>, N>& data)
{
   double result = 0;
   for(unsigned i = 0; i < N; ++i)
      result += boost::math::gamma_p_inv(data[i][0], data[i][5]);
   return result;
}

BOOST_MATH_PERFORMANCE_TEST(igamma_inv_test, "igamma_inv")
{
   double result = igamma_inv_evaluate2(igamma_big_data);
   result += igamma_inv_evaluate2(igamma_int_data);
   result += igamma_inv_evaluate2(igamma_med_data);
   result += igamma_inv_evaluate2(igamma_small_data);

   consume_result(result);
   set_call_count(
      (sizeof(igamma_big_data) 
      + sizeof(igamma_int_data) 
      + sizeof(igamma_med_data)
      + sizeof(igamma_small_data)) / sizeof(igamma_big_data[0]));
}
