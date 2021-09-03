#define BOOST_MATH_INSTRUMENT

#include <boost/cstdfloat.hpp>
#include <iostream>

#include <boost/math/special_functions/bessel.hpp>


int main()
{
   std::cout << std::setprecision(18) << boost::math::cyl_bessel_j(0, 11.791534423828125) << std::endl;
   std::cout << std::setprecision(35) << boost::math::cyl_bessel_j(0, 11.791534423828125L) << std::endl;
   std::cout << std::setprecision(35) << boost::math::cyl_bessel_j(0, 11.791534423828125Q) << std::endl;
   return 0;
}
