#include <boost/math/concepts/real_concept.hpp>
#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>

#include <utility>

#include <boost/array.hpp>
#include <boost/math/tools/polynomial.hpp>

using namespace boost::math;
using namespace std;

typedef tools::polynomial<double> PR;
typedef tools::polynomial<int> PZ;

BOOST_AUTO_TEST_CASE( test_main )
{
    boost::array<int, 4> const a = {3, -4, -6, 10};
    boost::array<int, 2> const b = {1, -2};
    PZ const x(a.data(), a.static_size), y(b.data(), b.static_size);
    
    std::pair<PZ, PZ> z = tools::quotient_remainder(x, y);
    
}