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
  
  // factorial primes 
  CHECK_EQUAL(true,detail::is_prime(719));
  CHECK_EQUAL(true,detail::is_prime(5039));
  CHECK_EQUAL(true,detail::is_prime(39916801));
  CHECK_EQUAL(true,detail::is_prime(479001599));
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

void test_power2()
{
  CHECK_EQUAL(false,detail::is_power2(-4));   
  CHECK_EQUAL(false,detail::is_power2(-3));   
  CHECK_EQUAL(false,detail::is_power2(-2));   
  CHECK_EQUAL(false,detail::is_power2(-1));   
  CHECK_EQUAL(false,detail::is_power2(0));   
  CHECK_EQUAL(false,detail::is_power2(3));   
  CHECK_EQUAL(false,detail::is_power2(5));   
  CHECK_EQUAL(false,detail::is_power2(6));   
  CHECK_EQUAL(false,detail::is_power2(7));   
  
  for(int i=9;i<16;++i)
    CHECK_EQUAL(false,detail::is_power2(i));   
  
  CHECK_EQUAL(true,detail::is_power2(1));   
  CHECK_EQUAL(true,detail::is_power2(2));   
  CHECK_EQUAL(true,detail::is_power2(4));   
  CHECK_EQUAL(true,detail::is_power2(8));   
  CHECK_EQUAL(true,detail::is_power2(16));   
  
  for(int i=-20;i<=0;++i)
    CHECK_EQUAL(0,detail::lower_bound_power2(i));   
    
  CHECK_EQUAL(1,detail::lower_bound_power2(1));   
  CHECK_EQUAL(2,detail::lower_bound_power2(2));   
  CHECK_EQUAL(2,detail::lower_bound_power2(3));   
  CHECK_EQUAL(4,detail::lower_bound_power2(4));   
  CHECK_EQUAL(4,detail::lower_bound_power2(5));   
  CHECK_EQUAL(4,detail::lower_bound_power2(6));   
  CHECK_EQUAL(4,detail::lower_bound_power2(7));   
  for(int i=8;i<16;++i)
    CHECK_EQUAL(8,detail::lower_bound_power2(i));   
  for(int i=16;i<32;++i)
    CHECK_EQUAL(16,detail::lower_bound_power2(i));   
  
  for(int i=-20;i<=0;++i)
    CHECK_EQUAL(1,detail::upper_bound_power2(i));   
    
  CHECK_EQUAL(1,detail::upper_bound_power2(1));   
  CHECK_EQUAL(2,detail::upper_bound_power2(2));   
  CHECK_EQUAL(4,detail::upper_bound_power2(3));   
  CHECK_EQUAL(4,detail::upper_bound_power2(4));   
  CHECK_EQUAL(8,detail::upper_bound_power2(5));   
  CHECK_EQUAL(8,detail::upper_bound_power2(6));   
  CHECK_EQUAL(8,detail::upper_bound_power2(7));   
  CHECK_EQUAL(8,detail::upper_bound_power2(8));   
  for(int i=9;i<=16;++i)
    CHECK_EQUAL(16,detail::upper_bound_power2(i));   
  for(int i=17;i<=32;++i)
    CHECK_EQUAL(32,detail::upper_bound_power2(i));   
  CHECK_EQUAL(0,detail::upper_bound_power2(std::numeric_limits<int>::max()));   
  CHECK_EQUAL(0,detail::upper_bound_power2((1<<30)+1));   
  CHECK_EQUAL(1<<30,detail::upper_bound_power2(1<<30));   
  CHECK_EQUAL(1<<30,detail::upper_bound_power2((1<<30) - 1));   
}

int main()
{
  test_power2();
  test_is_prime(); 
  test_primitive_root();
  return boost::math::test::report_errors();
}

