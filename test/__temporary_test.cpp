//
// This test file exists to output diagnostic info for tests failing in the online matrix
// for perplexing reasons, it's contents are subject to constant change!!
//
#define BOOST_MATH_INSTRUMENT

#include <boost/math/special_functions/ellint_d.hpp>
#include <iostream>

int main()
{
   std::cout << std::setprecision(20) << boost::math::ellint_d(-1.0, 1.0) << std::endl;
}

