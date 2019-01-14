//           Copyright Matthew Pulver 2018 - 2019.
// Distributed under the Boost Software License, Version 1.0.
//      (See accompanying file LICENSE_1_0.txt or copy at
//           https://www.boost.org/LICENSE_1_0.txt)

#include <boost/math/differentiation/autodiff.hpp>
#include <iostream>

using namespace boost::math::differentiation;

// Equations and function/variable names are from
// https://en.wikipedia.org/wiki/Greeks_(finance)#Formulas_for_European_option_Greeks

// Standard normal probability density function
template<typename X>
X phi(const X& x)
{
  return boost::math::constants::one_div_root_two_pi<double>()*exp(-0.5*x*x);
}

// Standard normal cumulative distribution function
template<typename X>
X Phi(const X& x)
{
  return 0.5*erfc(-boost::math::constants::one_div_root_two<double>()*x);
}

enum CP { call, put };

// Assume zero annual dividend yield (q=0).
template<typename Price,typename Sigma,typename Tau,typename Rate>
autodiff::promote<Price,Sigma,Tau,Rate>
    black_scholes_option_price(CP cp, double K, const Price& S, const Sigma& sigma, const Tau& tau, const Rate& r)
{
  using namespace std;
  const auto d1 = (log(S/K) + (r+sigma*sigma/2)*tau) / (sigma*sqrt(tau));
  const auto d2 = (log(S/K) + (r-sigma*sigma/2)*tau) / (sigma*sqrt(tau));
  if (cp == call)
    return S*Phi(d1) - exp(-r*tau)*K*Phi(d2);
  else
    return exp(-r*tau)*K*Phi(-d2) - S*Phi(-d1);
}

int main()
{
  const double K = 100.0; // Strike price.
  const autodiff::variable<double,3> S(105); // Stock price.
  const autodiff::variable<double,0,3> sigma(5); // Volatility.
  const autodiff::variable<double,0,0,1> tau(30.0/365); // Time to expiration in years. (30 days).
  const autodiff::variable<double,0,0,0,1> r(1.25/100); // Interest rate.
  const auto call_price = black_scholes_option_price(call, K, S, sigma, tau, r);
  const auto put_price  = black_scholes_option_price(put,  K, S, sigma, tau, r);

  // Compare automatically calculated greeks by autodiff with formulas for greeks.
  // https://en.wikipedia.org/wiki/Greeks_(finance)#Formulas_for_European_option_Greeks
  const double d1 = static_cast<double>((log(S/K) + (r+sigma*sigma/2)*tau) / (sigma*sqrt(tau)));
  const double d2 = static_cast<double>((log(S/K) + (r-sigma*sigma/2)*tau) / (sigma*sqrt(tau)));
  const double formula_call_delta = +Phi(+d1);
  const double formula_put_delta  = -Phi(-d1);
  const double formula_vega = static_cast<double>(S*phi(d1)*sqrt(tau));
  const double formula_call_theta = static_cast<double>(-S*phi(d1)*sigma/(2*sqrt(tau))-r*K*exp(-r*tau)*Phi(+d2));
  const double formula_put_theta  = static_cast<double>(-S*phi(d1)*sigma/(2*sqrt(tau))+r*K*exp(-r*tau)*Phi(-d2));
  const double formula_call_rho = static_cast<double>(+K*tau*exp(-r*tau)*Phi(+d2));
  const double formula_put_rho  = static_cast<double>(-K*tau*exp(-r*tau)*Phi(-d2));
  const double formula_gamma = static_cast<double>(phi(d1)/(S*sigma*sqrt(tau)));
  const double formula_vanna = static_cast<double>(-phi(d1)*d2/sigma);
  const double formula_charm = static_cast<double>(phi(d1)*(d2*sigma*sqrt(tau)-2*r*tau)/(2*tau*sigma*sqrt(tau)));
  const double formula_vomma = static_cast<double>(S*phi(d1)*sqrt(tau)*d1*d2/sigma);
  const double formula_veta = static_cast<double>(-S*phi(d1)*sqrt(tau)*(r*d1/(sigma*sqrt(tau))-(1+d1*d2)/(2*tau)));
  const double formula_speed = static_cast<double>(-phi(d1)*(d1/(sigma*sqrt(tau))+1)/(S*S*sigma*sqrt(tau)));
  const double formula_zomma = static_cast<double>(phi(d1)*(d1*d2-1)/(S*sigma*sigma*sqrt(tau)));
  const double formula_color =
    static_cast<double>(-phi(d1)/(2*S*tau*sigma*sqrt(tau))*(1+(2*r*tau-d2*sigma*sqrt(tau))*d1/(sigma*sqrt(tau))));
  const double formula_ultima = -formula_vega*static_cast<double>((d1*d2*(1-d1*d2)+d1*d1+d2*d2)/(sigma*sigma));

  std::cout << std::setprecision(std::numeric_limits<double>::digits10)
    << "autodiff black-scholes call price = " << call_price.derivative(0,0,0,0) << '\n'
    << "autodiff black-scholes put  price = " << put_price.derivative(0,0,0,0) << '\n'
    << "\n## First-order Greeks\n"
    << "autodiff call delta = " << call_price.derivative(1,0,0,0) << '\n'
    << " formula call delta = " << formula_call_delta << '\n'
    << "autodiff call vega  = " << call_price.derivative(0,1,0,0) << '\n'
    << " formula call vega  = " << formula_vega << '\n'
    << "autodiff call theta = " << -call_price.derivative(0,0,1,0) << '\n' // minus sign due to tau = T-time
    << " formula call theta = " << formula_call_theta << '\n'
    << "autodiff call rho   = " << call_price.derivative(0,0,0,1) << '\n'
    << " formula call rho   = " << formula_call_rho << '\n'
    << '\n'
    << "autodiff put delta = " << put_price.derivative(1,0,0,0) << '\n'
    << " formula put delta = " << formula_put_delta << '\n'
    << "autodiff put vega  = " << put_price.derivative(0,1,0,0) << '\n'
    << " formula put vega  = " << formula_vega << '\n'
    << "autodiff put theta = " << -put_price.derivative(0,0,1,0) << '\n'
    << " formula put theta = " << formula_put_theta << '\n'
    << "autodiff put rho   = " << put_price.derivative(0,0,0,1) << '\n'
    << " formula put rho   = " << formula_put_rho << '\n'
    << "\n## Second-order Greeks\n"
    << "autodiff call gamma = " << call_price.derivative(2,0,0,0) << '\n'
    << "autodiff put  gamma = " << put_price.derivative(2,0,0,0) << '\n'
    << "      formula gamma = " << formula_gamma << '\n'
    << "autodiff call vanna = " << call_price.derivative(1,1,0,0) << '\n'
    << "autodiff put  vanna = " << put_price.derivative(1,1,0,0) << '\n'
    << "      formula vanna = " << formula_vanna << '\n'
    << "autodiff call charm = " << -call_price.derivative(1,0,1,0) << '\n'
    << "autodiff put  charm = " << -put_price.derivative(1,0,1,0) << '\n'
    << "      formula charm = " << formula_charm << '\n'
    << "autodiff call vomma = " << call_price.derivative(0,2,0,0) << '\n'
    << "autodiff put  vomma = " << put_price.derivative(0,2,0,0) << '\n'
    << "      formula vomma = " << formula_vomma << '\n'
    << "autodiff call veta = " << call_price.derivative(0,1,1,0) << '\n'
    << "autodiff put  veta = " << put_price.derivative(0,1,1,0) << '\n'
    << "      formula veta = " << formula_veta << '\n'
    << "\n## Third-order Greeks\n"
    << "autodiff call speed = " << call_price.derivative(3,0,0,0) << '\n'
    << "autodiff put  speed = " << put_price.derivative(3,0,0,0) << '\n'
    << "      formula speed = " << formula_speed << '\n'
    << "autodiff call zomma = " << call_price.derivative(2,1,0,0) << '\n'
    << "autodiff put  zomma = " << put_price.derivative(2,1,0,0) << '\n'
    << "      formula zomma = " << formula_zomma << '\n'
    << "autodiff call color = " << call_price.derivative(2,0,1,0) << '\n'
    << "autodiff put  color = " << put_price.derivative(2,0,1,0) << '\n'
    << "      formula color = " << formula_color << '\n'
    << "autodiff call ultima = " << call_price.derivative(0,3,0,0) << '\n'
    << "autodiff put  ultima = " << put_price.derivative(0,3,0,0) << '\n'
    << "      formula ultima = " << formula_ultima << '\n'
    ;
  return 0;
}
/*
Output:
autodiff black-scholes call price = 56.5136030677739
autodiff black-scholes put  price = 51.4109161009333

## First-order Greeks
autodiff call delta = 0.773818444921273
 formula call delta = 0.773818444921274
autodiff call vega  = 9.05493427705736
 formula call vega  = 9.05493427705736
autodiff call theta = -275.73013426444
 formula call theta = -275.73013426444
autodiff call rho   = 2.03320550539396
 formula call rho   = 2.03320550539396

autodiff put delta = -0.226181555078726
 formula put delta = -0.226181555078726
autodiff put vega  = 9.05493427705736
 formula put vega  = 9.05493427705736
autodiff put theta = -274.481417851526
 formula put theta = -274.481417851526
autodiff put rho   = -6.17753255212599
 formula put rho   = -6.17753255212599

## Second-order Greeks
autodiff call gamma = 0.00199851912993254
autodiff put  gamma = 0.00199851912993254
      formula gamma = 0.00199851912993254
autodiff call vanna = 0.0410279463126531
autodiff put  vanna = 0.0410279463126531
      formula vanna = 0.0410279463126531
autodiff call charm = -1.2505564233679
autodiff put  charm = -1.2505564233679
      formula charm = -1.2505564233679
autodiff call vomma = -0.928114149313108
autodiff put  vomma = -0.928114149313108
      formula vomma = -0.928114149313107
autodiff call veta = 26.7947073115641
autodiff put  veta = 26.7947073115641
      formula veta = 26.7947073115641

## Third-order Greeks
autodiff call speed = -2.90117322380992e-05
autodiff put  speed = -2.90117322380992e-05
      formula speed = -2.90117322380992e-05
autodiff call zomma = -0.000604548369901419
autodiff put  zomma = -0.000604548369901419
      formula zomma = -0.000604548369901419
autodiff call color = -0.0184014426606065
autodiff put  color = -0.0184014426606065
      formula color = -0.0184014426606065
autodiff call ultima = -0.0922426864775683
autodiff put  ultima = -0.0922426864775683
      formula ultima = -0.0922426864775685
**/
