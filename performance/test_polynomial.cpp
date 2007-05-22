
#include "performance_measure.hpp"

#include <boost/math/tools/rational.hpp>
#include <boost/cstdint.hpp>

BOOST_MATH_PERFORMANCE_TEST(poly_test, "polynomial")
{
      typedef double T;
      static const T num[13] = {
         static_cast<T>(56906521.91347156388090791033559122686859L),
         static_cast<T>(103794043.1163445451906271053616070238554L),
         static_cast<T>(86363131.28813859145546927288977868422342L),
         static_cast<T>(43338889.32467613834773723740590533316085L),
         static_cast<T>(14605578.08768506808414169982791359218571L),
         static_cast<T>(3481712.15498064590882071018964774556468L),
         static_cast<T>(601859.6171681098786670226533699352302507L),
         static_cast<T>(75999.29304014542649875303443598909137092L),
         static_cast<T>(6955.999602515376140356310115515198987526L),
         static_cast<T>(449.9445569063168119446858607650988409623L),
         static_cast<T>(19.51992788247617482847860966235652136208L),
         static_cast<T>(0.5098416655656676188125178644804694509993L),
         static_cast<T>(0.006061842346248906525783753964555936883222L)
      };
      static const double denom[13] = {
         static_cast<boost::uint32_t>(0u),
         static_cast<boost::uint32_t>(39916800u),
         static_cast<boost::uint32_t>(120543840u),
         static_cast<boost::uint32_t>(150917976u),
         static_cast<boost::uint32_t>(105258076u),
         static_cast<boost::uint32_t>(45995730u),
         static_cast<boost::uint32_t>(13339535u),
         static_cast<boost::uint32_t>(2637558u),
         static_cast<boost::uint32_t>(357423u),
         static_cast<boost::uint32_t>(32670u),
         static_cast<boost::uint32_t>(1925u),
         static_cast<boost::uint32_t>(66u),
         static_cast<boost::uint32_t>(1u)
      };
   double result = 0;
   for(double i = 1; i < 100; i += 0.5)
      result += boost::math::tools::evaluate_polynomial(num, i) / boost::math::tools::evaluate_polynomial(denom, i);
   consume_result(result);
   set_call_count(200);
}

BOOST_MATH_PERFORMANCE_TEST(poly_mixed_test, "polynomial-mixed")
{
      typedef double T;
      static const T num[13] = {
         static_cast<T>(56906521.91347156388090791033559122686859L),
         static_cast<T>(103794043.1163445451906271053616070238554L),
         static_cast<T>(86363131.28813859145546927288977868422342L),
         static_cast<T>(43338889.32467613834773723740590533316085L),
         static_cast<T>(14605578.08768506808414169982791359218571L),
         static_cast<T>(3481712.15498064590882071018964774556468L),
         static_cast<T>(601859.6171681098786670226533699352302507L),
         static_cast<T>(75999.29304014542649875303443598909137092L),
         static_cast<T>(6955.999602515376140356310115515198987526L),
         static_cast<T>(449.9445569063168119446858607650988409623L),
         static_cast<T>(19.51992788247617482847860966235652136208L),
         static_cast<T>(0.5098416655656676188125178644804694509993L),
         static_cast<T>(0.006061842346248906525783753964555936883222L)
      };
      static const boost::uint32_t denom[13] = {
         static_cast<boost::uint32_t>(0u),
         static_cast<boost::uint32_t>(39916800u),
         static_cast<boost::uint32_t>(120543840u),
         static_cast<boost::uint32_t>(150917976u),
         static_cast<boost::uint32_t>(105258076u),
         static_cast<boost::uint32_t>(45995730u),
         static_cast<boost::uint32_t>(13339535u),
         static_cast<boost::uint32_t>(2637558u),
         static_cast<boost::uint32_t>(357423u),
         static_cast<boost::uint32_t>(32670u),
         static_cast<boost::uint32_t>(1925u),
         static_cast<boost::uint32_t>(66u),
         static_cast<boost::uint32_t>(1u)
      };
   double result = 0;
   for(double i = 1; i < 100; i += 0.5)
      result += boost::math::tools::evaluate_polynomial(num, i) / boost::math::tools::evaluate_polynomial(denom, i);
   consume_result(result);
   set_call_count(200);
}

BOOST_MATH_PERFORMANCE_TEST(rat_test, "rational")
{
      typedef double T;
      static const T num[13] = {
         static_cast<T>(56906521.91347156388090791033559122686859L),
         static_cast<T>(103794043.1163445451906271053616070238554L),
         static_cast<T>(86363131.28813859145546927288977868422342L),
         static_cast<T>(43338889.32467613834773723740590533316085L),
         static_cast<T>(14605578.08768506808414169982791359218571L),
         static_cast<T>(3481712.15498064590882071018964774556468L),
         static_cast<T>(601859.6171681098786670226533699352302507L),
         static_cast<T>(75999.29304014542649875303443598909137092L),
         static_cast<T>(6955.999602515376140356310115515198987526L),
         static_cast<T>(449.9445569063168119446858607650988409623L),
         static_cast<T>(19.51992788247617482847860966235652136208L),
         static_cast<T>(0.5098416655656676188125178644804694509993L),
         static_cast<T>(0.006061842346248906525783753964555936883222L)
      };
      static const double denom[13] = {
         static_cast<boost::uint32_t>(0u),
         static_cast<boost::uint32_t>(39916800u),
         static_cast<boost::uint32_t>(120543840u),
         static_cast<boost::uint32_t>(150917976u),
         static_cast<boost::uint32_t>(105258076u),
         static_cast<boost::uint32_t>(45995730u),
         static_cast<boost::uint32_t>(13339535u),
         static_cast<boost::uint32_t>(2637558u),
         static_cast<boost::uint32_t>(357423u),
         static_cast<boost::uint32_t>(32670u),
         static_cast<boost::uint32_t>(1925u),
         static_cast<boost::uint32_t>(66u),
         static_cast<boost::uint32_t>(1u)
      };
   double result = 0;
   for(double i = 1; i < 100; i += 0.5)
      result += boost::math::tools::evaluate_rational(num, denom, i);
   consume_result(result);
   set_call_count(100);
}

BOOST_MATH_PERFORMANCE_TEST(rat_mixed_test, "rational-mixed")
{
      typedef double T;
      static const T num[13] = {
         static_cast<T>(56906521.91347156388090791033559122686859L),
         static_cast<T>(103794043.1163445451906271053616070238554L),
         static_cast<T>(86363131.28813859145546927288977868422342L),
         static_cast<T>(43338889.32467613834773723740590533316085L),
         static_cast<T>(14605578.08768506808414169982791359218571L),
         static_cast<T>(3481712.15498064590882071018964774556468L),
         static_cast<T>(601859.6171681098786670226533699352302507L),
         static_cast<T>(75999.29304014542649875303443598909137092L),
         static_cast<T>(6955.999602515376140356310115515198987526L),
         static_cast<T>(449.9445569063168119446858607650988409623L),
         static_cast<T>(19.51992788247617482847860966235652136208L),
         static_cast<T>(0.5098416655656676188125178644804694509993L),
         static_cast<T>(0.006061842346248906525783753964555936883222L)
      };
      static const boost::uint32_t denom[13] = {
         static_cast<boost::uint32_t>(0u),
         static_cast<boost::uint32_t>(39916800u),
         static_cast<boost::uint32_t>(120543840u),
         static_cast<boost::uint32_t>(150917976u),
         static_cast<boost::uint32_t>(105258076u),
         static_cast<boost::uint32_t>(45995730u),
         static_cast<boost::uint32_t>(13339535u),
         static_cast<boost::uint32_t>(2637558u),
         static_cast<boost::uint32_t>(357423u),
         static_cast<boost::uint32_t>(32670u),
         static_cast<boost::uint32_t>(1925u),
         static_cast<boost::uint32_t>(66u),
         static_cast<boost::uint32_t>(1u)
      };
   double result = 0;
   for(double i = 1; i < 100; i += 0.5)
      result += boost::math::tools::evaluate_rational(num, denom, i);
   consume_result(result);
   set_call_count(100);
}

