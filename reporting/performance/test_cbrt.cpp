//  Copyright John Maddock 2015.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifdef _MSC_VER
#  pragma warning (disable : 4224)
#endif

#include <boost/math/special_functions/cbrt.hpp>
#include <boost/array.hpp>
#include "../../test/table_type.hpp"
#include "table_helper.hpp"
#include "performance.hpp"
#include <iostream>

int main()
{
   typedef double T;
#define SC_(x) static_cast<double>(x)
#  include "../../test/cbrt_data.ipp"

   add_data(cbrt_data);


   screen_data([](const std::vector<double>& v){  return boost::math::cbrt(v[1]);  }, [](const std::vector<double>& v){ return v[0];  });

#ifdef TEST_C99
   screen_data([](const std::vector<double>& v){  return ::cbrt(v[1]);  }, [](const std::vector<double>& v){ return v[0];  });
#endif

   double time = exec_timed_test([](const std::vector<double>& v){  return boost::math::cbrt(v[1]);  });
   report_execution_time_multi_compilation(time, "cbrt", "boost::math::cbrt", get_compiler_options_name());
   std::cout << time << std::endl;
#ifdef TEST_C99
   time = exec_timed_test([](const std::vector<double>& v){  return ::cbrt(v[1]);  });
   std::cout << time << std::endl;
   report_execution_time_multi_compilation(time, "cbrt", "::cbrt", get_compiler_options_name());
#endif


   return 0;
}

