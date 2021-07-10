///////////////////////////////////////////////////////////////////
//  Copyright Eduardo Quintana 2021
//  Copyright Janek Kozicki 2021
//  Copyright Christopher Kormanyos 2021
//  Distributed under the Boost Software License,
//  Version 1.0. (See accompanying file LICENSE_1_0.txt
//  or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <algorithm>
#include <limits>
#include <type_traits>

#ifndef BOOST_MATH_FFT_DISCRETE_HPP
  #define BOOST_MATH_FFT_DISCRETE_HPP
  
  namespace boost { namespace math {  namespace fft {
  namespace detail {
  
  template <typename integer>
  integer modulo(integer a, integer b)
  // it computes 'a mod b' in the mathematical sense, not 'a%b'
  // it can handle a<0 values
  // precondition: b>0
  {
    return ( a%b + b ) % b;
  }
  
  template <typename integer>
  struct euclid_return
  // helper class for the extended euclidean algorithm
  {
    integer gcd, x, y;
  };
  
  template <typename integer>
  euclid_return<integer> extended_euclid(integer a, integer b)
  // extended euclidean algorithm. It solves Bezout equation for ANY integer
  // pair (a,b):
  // g = a x  + b y
  {
    using std::swap;
    static_assert(std::numeric_limits<integer>::is_signed,
                  "Signed integers are needed");
    integer x = 1, y = 0, g = a;
    for (integer x1 = 0, y1 = 1, g1 = b; g1;)
    {
      integer q = g / g1;
      x = x - x1 * q;
      swap(x, x1);
      y = y - y1 * q;
      swap(y, y1);
      g = g - g1 * q;
      swap(g, g1);
    }
    return {g, x, y};
  }
  
  template <typename integer>
  integer gcd(integer a, integer b)
  // Greatest Commond Divisor
  {
    euclid_return<integer> ret = extended_euclid(a,b);
    return ret.gcd;
  }
  
  template <typename T,typename integer>
  T power(const T& x, integer n)
    // WARNING: 
    // -> for n<0 we cannot decide the outcome, because we don't know the inverse of x
    // -> for n=0 we cannot decide the outcome, because we don't know the
    // multiplication neutral element of T.
  // precondition: n>0
  {
    if(n<=0)
      return x;
      
    bool identity = true;
    T r{x};

    for (T aux{x}; n; n /= 2)
    {
      if (n % 2)
      {
        r = identity ? aux : r * aux;
        identity = false;
      }
      aux *= aux;
    }
    return r;
  }
  
  template<typename integer>
  constexpr 
  typename std::enable_if<std::is_integral<integer>::value,bool>::type 
  // enabled only for builtin integer types
  is_power2(integer x) 
  { 
    return (x>0) && (x == (x & -x));
  }
  
  template<typename integer>
  typename std::enable_if<std::is_integral<integer>::value,integer>::type 
  // enabled only for builtin integer types
  lower_bound_power2(integer x)
  // Returns the biggest power of two that is smaller than or equal to x.
  // returns 1 if no such power of two exists, ie. for x<=0
  {
    if(x<=0)
      return 0;
      
    for(integer lsb = x & -x; lsb!=x;)
    {
      x-=lsb;
      lsb = x&-x;
    }
    return x;
  }
  
  
  template<typename integer>
  typename std::enable_if<std::is_integral<integer>::value,integer>::type 
  // enabled only for builtin integer types
  upper_bound_power2(integer x)
  // Returns the smallest power of two that is greater or equal to x.
  // returns 0 if no such power of two exits, ie. when the value overflows
  // the integer
  {
    if(x<=1)
      return 1;
    
    const integer up_max = lower_bound_power2(
      std::numeric_limits<integer>::max() );
    if(up_max<x)
      return 0;
    
    integer up=1;
    for(;up<x;up<<=1);
    return up;
  }
  template <typename integer, typename Exp>
  bool power_greater(integer A, integer B, Exp p)
    // precondition: A,B,p > 0
    // returns true if A^p > B
  {
    if (A == 1)
      return 1 > B;
    if (p > std::numeric_limits<integer>::digits)
      return true;

    while (p--)
      B /= A;
    return 1 > B;
  }
  template <typename integer, typename exp_t>
  integer root(const integer N, const exp_t p)
    // binary search to find the biggest r such that r^p <= N,
    // ie. r = floor( pRoot(N)  )
  {
    if (p == 1)
      return N;
    // no overflow
    integer down = 1, up = N;
    while (up - down > 1)
    {
      integer mid = down + (up - down) / 2;
      if (power_greater(mid, N, p))
        up = mid;
      else
        down = mid;
    }
    return down;
  }
  template <typename integer>
  class mint
  /*
    Helper class for modular arithmetics,
    with modulo known at runtime.
    WARNING: the big flaw is that there are no checks for the homogeneity
    of the modulo, ie. you could a+b where a.mod != b.mod
  */
  {
    static bool overflow_by_multiplication(integer m)
    // checks if the give modular base 'm' can overflow the direct
    // 'integer' multiplication
    {
      integer sqrt_max = root(std::numeric_limits<integer>::max(), 2);
      return m - 1 > sqrt_max;
    }
  
    const integer mod{};
    integer x{};

   public:
    mint(){}
    mint(integer m,integer _x) :  mod{m}, x{_x}
    {
      x = modulo(x,mod);
    }
    mint(const mint& that) : mod{that.mod}, x{that.x} {}

    mint& operator=(const mint& that) { return x = that.x, *this; }

    explicit operator bool() const { return x == integer{0}; }
    operator integer() const { return x; }

    mint& operator++() { return x = (x + 1) % mod, *this; }
    mint& operator--() { return x = (mod + x - 1) % mod, *this;}
    mint operator++(int)
    { /*post*/
      mint temp(*this);
      this->operator++();
      return temp;
    }
    mint operator--(int)
    { /*post*/
      mint temp(*this);
      this->operator--();
      return temp;
    }
    mint& operator+=(const mint& that){return x = (x + that.x) % mod, *this;} 
    mint& operator-=(const mint& t){return x = (mod + x - t.x) % mod, *this;}

    mint& operator*=(const mint& t)
    {
      if (overflow_by_multiplication(mod))
      {
        /// multiplication based on sums
        integer res{0};
        for (integer a = x, b = t.x; b; b /= 2)
        {
          if (b % 2)
            res = (res + a) % mod;
          a = (a + a) % mod;
        }
        x = res;
        return *this;
      }
      // direct multiplication
      x = (x * t.x) % mod;
      return *this;
    }

    mint<integer> inverse() const
    {
        euclid_return<integer> ee = extended_euclid(x, mod);
        if (ee.gcd != 1)
            return mint<integer>(0,0); // no inverse found
        integer inv_x = ee.x;
        return mint(mod,inv_x);
        // return power(*this,mod-2);
    }
    bool operator==(const mint& that) const { return x == that.x; }
    bool operator<(const mint& that) const { return x < that.x; }
  };
  template <typename T>
  mint<T> operator*(const mint<T>& A, const mint<T>& B)
  {
    mint<T> C{A};
    return C *= B;
  }
  
  template<typename integer>
  integer power_mod(integer x, integer n, integer m)
  // Computes x^n mod m
  {
    if(x==0 || x==1)
        return x;
    if(n<=0)
      return 1;
    mint<integer> xm(m,x);
    // n = modulo(n,euler_phi(m));
    xm = power(xm,n);
    return integer(xm);
  }
  
  inline bool is_prime(int n)
  {
    // Naive O(sqrt(n)) prime decision
    // TODO: replace with something more sophisticated like Rabin-Miller's test
    
    if(n<2)
      return false;
    if(n==2 || n==3)
      return true;
    
    const int sqrt_n = root(n,2);
    for(int x=2;x<=sqrt_n;++x)
    {
      if(n%x == 0)
        return false;
    }
    return true;
  }
  
  template<class Iterator>
  int prime_factorization(int n,Iterator out, const bool unique=false)
  // Naive O(sqrt(n)) prime factorization.
  // returns a list of prime numbers that when multiplied they give n
  {
    int count{0};
    const int sqrt_n = root(n,2);
    for (int x = 2; x <= sqrt_n;++x)
    {
      if(n % x == 0)
      {
        *out = x;
        ++out;++count;
        n /= x;
      }
      if(unique)
      {
        while(n % x == 0) n /= x;
      }else
      {
        while(n % x == 0)
        { 
          *out = x;
          ++out;++count;
          n /= x;
        }
      }
    }
    if (n > 1)
    {
      *out = n;
      ++out;++count;
    }
    return count;
  }
  
  template<class Iterator>
  int prime_factorization_unique(int n, Iterator out)
  // returns the list of unique prime factors divisors of n
  {
    return prime_factorization(n,out,true);
  }

  inline int euler_phi(int n)
  // Euler Phi function
  {
    int r = n;
    std::array<int,32> F;
    const int count = prime_factorization_unique(n,F.begin());
    for (int i=0;i<count;++i)
        r -= r / F[i];
    return r;
  }
  
  inline int primitive_root(int n)
  // it finds a primitive root r of n,
  // ie. r^phi(n) = 1 mod n
  {
    const int phi = euler_phi(n);
    std::array<int,32> F;
    const int count = prime_factorization_unique(phi,F.begin());
    for(int i=1;i<n;++i)
    {
      int g = gcd(i,n);
      if(g!=1) continue;
      
      bool ok = true;
      for(int j=0;j<count;++j)
      {
        if(power_mod(i,phi/F[j],n)==1) // not a root
        {
          ok=false;
          break;
        }
      }
      
      if(ok)
        return i;
    }
    return 0; // no roots found
  }

  } // namespace detail
  } } } // namespace boost::math::fft
#endif
