#include "math_unit_test.hpp"
#include <boost/math/fft/algorithms.hpp>

using namespace boost::math::fft;

void test_is_prime()
{
  CHECK_EQUAL(false,detail::is_prime(1));
  CHECK_EQUAL(true,detail::is_prime(2));
  CHECK_EQUAL(true,detail::is_prime(3));
  CHECK_EQUAL(false,detail::is_prime(4));
  CHECK_EQUAL(true,detail::is_prime(5));
  CHECK_EQUAL(false,detail::is_prime(6));
  CHECK_EQUAL(true,detail::is_prime(7));
  CHECK_EQUAL(false,detail::is_prime(8));
  CHECK_EQUAL(false,detail::is_prime(9));
  CHECK_EQUAL(false,detail::is_prime(10));
}
void test_primitive_root()
{
  for(auto p : std::vector<long>{3,5,7,11,13,17,23,29})
  {
    const long r = detail::primitive_root(p);
    std::cerr << "p = " << p << ", root = " << r << '\n';
    const long phi = p-1;
    
    long r_p = r;
    for(int i=1;i<phi;++i)
    {
      CHECK_EQUAL(false,r_p == 1);
      r_p = (r_p * r) % p;
    }
      CHECK_EQUAL(true,r_p == 1);
  }
}
int main()
{
  test_is_prime(); 
  test_primitive_root();
  return boost::math::test::report_errors();
}

