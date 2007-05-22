
#include "performance_measure.hpp"

#include <boost/math/special_functions/gamma.hpp>
#include <boost/array.hpp>

#define T double
#include "../test/erf_data.ipp"
#include "../test/erf_large_data.ipp"
#include "../test/erf_small_data.ipp"

template <std::size_t N>
double erf_evaluate2(const boost::array<boost::array<T, 3>, N>& data)
{
   double result = 0;
   for(unsigned i = 0; i < N; ++i)
      result += boost::math::erf(data[i][0]);
   return result;
}

BOOST_MATH_PERFORMANCE_TEST(erf_test, "erf")
{
   double result = erf_evaluate2(erf_data);
   result += erf_evaluate2(erf_large_data);
   result += erf_evaluate2(erf_small_data);

   consume_result(result);
   set_call_count((sizeof(erf_data) + sizeof(erf_large_data) + sizeof(erf_small_data)) / sizeof(erf_data[0]));
}

template <std::size_t N>
double erf_inv_evaluate2(const boost::array<boost::array<T, 3>, N>& data)
{
   double result = 0;
   for(unsigned i = 0; i < N; ++i)
      result += boost::math::erf_inv(data[i][1]);
   return result;
}

BOOST_MATH_PERFORMANCE_TEST(erf_inv_test, "erf_inv")
{
   double result = erf_inv_evaluate2(erf_data);
   result += erf_inv_evaluate2(erf_large_data);
   result += erf_inv_evaluate2(erf_small_data);

   consume_result(result);
   set_call_count((sizeof(erf_data) + sizeof(erf_large_data) + sizeof(erf_small_data)) / sizeof(erf_data[0]));
}