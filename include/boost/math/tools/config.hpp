#ifndef BOOST_MATH_TOOLS_CONFIG_HPP
#define BOOST_MATH_TOOLS_CONFIG_HPP

#include <boost/math/tools/error_handling.hpp>
#include <boost/cstdint.hpp> // for boost::uintmax_t

#define BOOST_MATH_MAX_ITER 1000000

namespace boost{ namespace math{
namespace tools
{

inline void check_series_iterations(const char* function, boost::uintmax_t max_iter)
{
   if(max_iter >= BOOST_MATH_MAX_ITER)
      tools::logic_error<boost::uintmax_t>(
         function,
         "Series evaluation exceeded %1% iterations, giving up now.", max_iter);
}

} // namespace tools
}} // namespace boost namespace math

#ifdef __linux__

	#include <fenv.h>

	namespace boost{ namespace math{
	namespace detail
	{
	struct fpu_guard
	{
		fpu_guard()
		{
			fegetexceptflag(&m_flags, FE_ALL_EXCEPT);
			feclearexcept(FE_ALL_EXCEPT);
		}
		~fpu_guard()
		{
			fesetexceptflag(&m_flags, FE_ALL_EXCEPT);
		}
	private:
		fexcept_t m_flags;
	};

	} // namespace detail
	}} // namespaces

	#define BOOST_FPU_EXCEPTION_GUARD boost::math::detail::fpu_guard local_guard_object;
#else // All other platforms.
  #define BOOST_FPU_EXCEPTION_GUARD
#endif

#endif // BOOST_MATH_TOOLS_CONFIG_HPP
