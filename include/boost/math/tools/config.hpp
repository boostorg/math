
#ifndef BOOST_MATH_TOOLS_CONFIG_HPP
#define BOOST_MATH_TOOLS_CONFIG_HPP

#include <boost/math/tools/error_handling.hpp>

#define BOOST_MATH_MAX_ITER 1000000

namespace boost{ namespace math{ namespace tools{

inline void check_series_iterations(const char* function, boost::uintmax_t max_iter)
{
   if(max_iter >= BOOST_MATH_MAX_ITER)
      tools::logic_error<boost::uintmax_t>(
         function, 
         "Series evaluation exceeded %1% iterations, giving up now.", max_iter);
}

}}} // namespaces

#endif // BOOST_MATH_TOOLS_CONFIG_HPP
