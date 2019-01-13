//           Copyright Matthew Pulver 2018 - 2019.
// Distributed under the Boost Software License, Version 1.0.
//      (See accompanying file LICENSE_1_0.txt or copy at
//           https://www.boost.org/LICENSE_1_0.txt)

#include <boost/math/differentiation/autodiff.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <iostream>

template<typename T>
T f(const T& w, const T& x, const T& y, const T& z)
{
  using namespace std;
  return exp(w*sin(x*log(y)/z) + sqrt(w*z/(x*y))) + w*w/tan(z);
}

int main()
{
  using cpp_dec_float_100 = boost::multiprecision::cpp_dec_float_100;
  using namespace boost::math::differentiation;

  constexpr int Nw=3; // Max order of derivative to calculate for w
  constexpr int Nx=2; // Max order of derivative to calculate for x
  constexpr int Ny=4; // Max order of derivative to calculate for y
  constexpr int Nz=3; // Max order of derivative to calculate for z
  using var = autodiff::variable<cpp_dec_float_100,Nw,Nx,Ny,Nz>;
  const var w = autodiff::variable<cpp_dec_float_100,Nw>(11);
  const var x = autodiff::variable<cpp_dec_float_100,0,Nx>(12);
  const var y = autodiff::variable<cpp_dec_float_100,0,0,Ny>(13);
  const var z = autodiff::variable<cpp_dec_float_100,0,0,0,Nz>(14);
  const var v = f(w,x,y,z);
  // Calculated from Mathematica symbolic differentiation. See multiprecision.nb for script.
  const cpp_dec_float_100 answer("1976.31960074779771777988187529041872090812118921875499076582535951111845769110560421820940516423255314");
  std::cout << std::setprecision(std::numeric_limits<cpp_dec_float_100>::digits10)
    << "mathematica   : " << answer << '\n'
    << "autodiff      : " << v.derivative(Nw,Nx,Ny,Nz) << '\n'
    << "relative error: " << std::setprecision(3) << (v.derivative(Nw,Nx,Ny,Nz)/answer-1) << std::endl;
  return 0;
}
/*
Output:
mathematica   : 1976.319600747797717779881875290418720908121189218754990765825359511118457691105604218209405164232553
autodiff      : 1976.319600747797717779881875290418720908121189218754990765825359511118457691105604218209405164232566
relative error: 6.47e-99
**/
