#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>

#include <utility>

#include <boost/array.hpp>
#include <boost/math/tools/polynomial.hpp>

using namespace boost::math::tools;
using namespace std;

template <typename T>
struct question
{
    polynomial<T> dividend;
    polynomial<T> divisor;
};

template <typename T>
struct answer
{
    answer(std::pair< polynomial<T>, polynomial<T> > const &x) :
        quotient(x.first), remainder(x.second) {}
        
    polynomial<T> quotient;
    polynomial<T> remainder;
};

typedef polynomial<int> PZ;

BOOST_AUTO_TEST_CASE( test_main )
{
    boost::array<int, 4> const d3 = {3, -4, -6, 10};
    boost::array<int, 2> const d1 = {1, -2};
    boost::array<int, 3> const d3_div_d1 = {3, 2, -2};
    boost::array<int, 1> const d3_rem_d1 = {6};
    PZ const x(d3.data(), d3.static_size), y(d1.data(), d1.static_size);
    PZ const zero;
    PZ const q(d3_div_d1.data(), d3_div_d1.static_size);
    PZ const r(d3_rem_d1.data(), d3_rem_d1.static_size);
    
    // boost::array<pair< question<int>, answer<int> >, 1> q_and_a;
    answer<int> result = quotient_remainder(x, y);
    BOOST_CHECK_EQUAL(result.quotient, q);
    BOOST_CHECK_EQUAL(result.remainder, r);
}
