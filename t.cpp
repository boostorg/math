
#define BOOST_MATH_INSTRUMENT

#include <boost/math/special_functions/bessel.hpp>


int main()
{
   long double a = -1000;
   long double b = 700;
   long double expect = 6.51561979144735818903553852606382e-31L;
   long double result = boost::math::cyl_bessel_k(a, b);
   std::cout << std::setprecision(std::numeric_limits<long double>::max_digits10) << result << std::endl;
   long double err = fabs((expect - result) / expect);
   std::cout << err << std::endl;
   std::cout << err / std::numeric_limits<long double>::epsilon() << std::endl;
}
