
#include "performance_measure.hpp"

#include <boost/math/special_functions/beta.hpp>
#include <boost/array.hpp>

#define T double
#include "../test/ibeta_data.ipp"
#include "../test/ibeta_int_data.ipp"
#include "../test/ibeta_large_data.ipp"
#include "../test/ibeta_small_data.ipp"

template <std::size_t N>
double ibeta_evaluate2(const boost::array<boost::array<T, 7>, N>& data)
{
   double result = 0;
   for(unsigned i = 0; i < N; ++i)
      result += boost::math::ibeta(data[i][0], data[i][1], data[i][2]);
   return result;
}

BOOST_MATH_PERFORMANCE_TEST(ibeta_test, "ibeta")
{
   double result = ibeta_evaluate2(ibeta_data);
   result += ibeta_evaluate2(ibeta_int_data);
   result += ibeta_evaluate2(ibeta_large_data);
   result += ibeta_evaluate2(ibeta_small_data);

   consume_result(result);
   set_call_count(
      (sizeof(ibeta_data) 
      + sizeof(ibeta_int_data) 
      + sizeof(ibeta_large_data)
      + sizeof(ibeta_small_data)) / sizeof(ibeta_data[0]));
}

template <std::size_t N>
double ibeta_inv_evaluate2(const boost::array<boost::array<T, 7>, N>& data)
{
   double result = 0;
   for(unsigned i = 0; i < N; ++i)
      result += boost::math::ibeta_inv(data[i][0], data[i][1], data[i][5]);
   return result;
}

BOOST_MATH_PERFORMANCE_TEST(ibeta_inv_test, "ibeta_inv")
{
   double result = ibeta_inv_evaluate2(ibeta_data);
   result += ibeta_inv_evaluate2(ibeta_int_data);
   result += ibeta_inv_evaluate2(ibeta_large_data);
   result += ibeta_inv_evaluate2(ibeta_small_data);

   consume_result(result);
   set_call_count(
      (sizeof(ibeta_data) 
      + sizeof(ibeta_int_data) 
      + sizeof(ibeta_large_data)
      + sizeof(ibeta_small_data)) / sizeof(ibeta_data[0]));
}
