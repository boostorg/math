//           Copyright Maksym Zhelyenzyakov 2025-2026.
// Distributed under the Boost Software License, Version 1.0.
//      (See accompanying file LICENSE_1_0.txt or copy at
//           https://www.boost.org/LICENSE_1_0.txt)
#include <boost/math/differentiation/autodiff_reverse.hpp>

using namespace boost::math::differentiation::reverse_mode;
using namespace boost::math::constants;

template<typename Real>
Real phi(Real const& x)
{
    return one_div_root_two_pi<Real>() * exp(-0.5 * x * x);
}

template<typename Real>
Real Phi(Real const& x)
{
    return 0.5 * erfc(-one_div_root_two<Real>() * x);
}

enum class CP { call, put };

template<typename T>
T black_scholes_option_price(CP cp, double K, T const& S, T const& sigma, T const& tau, T const& r)
{
    using namespace std;
    auto const d1 = (log(S / K) + (r + sigma * sigma / 2) * tau) / (sigma * sqrt(tau));
    auto const d2 = (log(S / K) + (r - sigma * sigma / 2) * tau) / (sigma * sqrt(tau));
    switch (cp) {
    case CP::call:
        return S * Phi<T>(d1) - exp(-r * tau) * K * Phi<T>(d2);
    case CP::put:
        return exp(-r * tau) * K * Phi<T>(-d2) - S * Phi<T>(-d1);
    default:
        throw runtime_error("Invalid CP value.");
    }
}

int main()
{
    double const    K         = 100.0;
    double          S_val     = 105.0;
    double          sigma_val = 5.0;
    double          tau_val   = 30.0 / 365;
    double          r_val     = 1.25 / 100;
    rvar<double, 3> S         = make_rvar<double, 3>(S_val);
    rvar<double, 3> sigma     = make_rvar<double, 3>(sigma_val);
    rvar<double, 3> tau       = make_rvar<double, 3>(tau_val);
    rvar<double, 3> r         = make_rvar<double, 3>(r_val);

    rvar<double, 3> call_price
        = black_scholes_option_price<rvar<double, 3>>(CP::call, K, S, sigma, tau, r);
    rvar<double, 3> put_price
        = black_scholes_option_price<rvar<double, 3>>(CP::put, K, S, sigma, tau, r);

    double const d1 = ((log(S_val / K) + (r_val + sigma_val * sigma_val / 2) * tau_val)
                       / (sigma_val * sqrt(tau_val)));
    double const d2 = ((log(S_val / K) + (r_val - sigma_val * sigma_val / 2) * tau_val)
                       / (sigma_val * sqrt(tau_val)));
    double const formula_call_delta = +Phi(+d1);
    double const formula_put_delta  = -Phi(-d1);
    double const formula_vega       = (S_val * phi(d1) * sqrt(tau_val));
    double const formula_call_theta = (-S_val * phi(d1) * sigma_val / (2 * sqrt(tau_val))
                                       - r_val * K * exp(-r_val * tau_val) * Phi(+d2));

    double const formula_put_theta = (-S_val * phi(d1) * sigma_val / (2 * sqrt(tau_val))
                                      + r_val * K * exp(-r_val * tau_val) * Phi(-d2));
    double const formula_call_rho  = (+K * tau_val * exp(-r_val * tau_val) * Phi(+d2));
    double const formula_put_rho   = (-K * tau_val * exp(-r_val * tau_val) * Phi(-d2));
    double const formula_gamma     = (phi(d1) / (S_val * sigma_val * sqrt(tau_val)));
    double const formula_vanna     = (-phi(d1) * d2 / sigma_val);
    double const formula_charm  = (phi(d1) * (d2 * sigma_val * sqrt(tau_val) - 2 * r_val * tau_val)
                                  / (2 * tau_val * sigma_val * sqrt(tau_val)));
    double const formula_vomma  = (S_val * phi(d1) * sqrt(tau_val) * d1 * d2 / sigma_val);
    double const formula_veta   = (-S_val * phi(d1) * sqrt(tau_val)
                                 * (r_val * d1 / (sigma_val * sqrt(tau_val))
                                    - (1 + d1 * d2) / (2 * tau_val)));
    double const formula_speed  = (-phi(d1) * (d1 / (sigma_val * sqrt(tau_val)) + 1)
                                  / (S_val * S_val * sigma_val * sqrt(tau_val)));
    double const formula_zomma  = (phi(d1) * (d1 * d2 - 1)
                                  / (S_val * sigma_val * sigma_val * sqrt(tau_val)));
    double const formula_color  = (-phi(d1) / (2 * S_val * tau_val * sigma_val * sqrt(tau_val))
                                  * (1
                                     + (2 * r_val * tau_val - d2 * sigma_val * sqrt(tau_val)) * d1
                                           / (sigma_val * sqrt(tau_val))));
    double const formula_ultima = -formula_vega
                                  * ((d1 * d2 * (1 - d1 * d2) + d1 * d1 + d2 * d2)
                                     / (sigma_val * sigma_val));

    auto call_greeks = grad(call_price, &S, &sigma, &tau, &r);
    auto put_greeks  = grad(put_price, &S, &sigma, &tau, &r);

    auto call_greeks_2nd_order = grad_nd<2>(call_price, &S, &sigma, &tau, &r);
    auto put_greeks_2nd_order  = grad_nd<2>(put_price, &S, &sigma, &tau, &r);

    auto call_greeks_3rd_order = grad_nd<3>(call_price, &S, &sigma, &tau, &r);
    auto put_greeks_3rd_order  = grad_nd<3>(put_price, &S, &sigma, &tau, &r);

    std::cout << std::setprecision(std::numeric_limits<double>::digits10)
              << "autodiff black-scholes call price = " << call_price.item() << "\n"
              << "autodiff black-scholes put price = " << put_price.item() << "\n"
              << "\n ## First-order Greeks \n"
              << "autodiff call delta = " << call_greeks[0].item() << "\n"
              << "formula call delta = " << formula_call_delta << "\n"
              << "autodiff call vega  = " << call_greeks[1].item() << '\n'
              << " formula call vega  = " << formula_vega << '\n'
              << "autodiff call theta = " << -call_greeks[2].item() << '\n'
              << " formula call theta = " << formula_call_theta << '\n'
              << "autodiff call rho   = " << call_greeks[3].item() << 'n'
              << " formula call rho   = " << formula_call_rho << '\n'
              << '\n'
              << "autodiff put delta = " << put_greeks[0].item() << 'n'
              << " formula put delta = " << formula_put_delta << '\n'
              << "autodiff put vega  = " << put_greeks[1].item() << '\n'
              << " formula put vega  = " << formula_vega << '\n'
              << "autodiff put theta = " << -put_greeks[2].item() << '\n'
              << " formula put theta = " << formula_put_theta << '\n'
              << "autodiff put rho   = " << put_greeks[3].item() << '\n'
              << " formula put rho   = " << formula_put_rho << '\n'

              << "\n## Second-order Greeks\n"
              << "autodiff call gamma = " << call_greeks_2nd_order[0][0].item() << '\n'
              << "autodiff put  gamma = " << put_greeks_2nd_order[0][0].item() << '\n'
              << "      formula gamma = " << formula_gamma << '\n'
              << "autodiff call vanna = " << call_greeks_2nd_order[0][1].item() << '\n'
              << "autodiff put  vanna = " << put_greeks_2nd_order[0][1].item() << '\n'
              << "      formula vanna = " << formula_vanna << '\n'
              << "autodiff call charm = " << -call_greeks_2nd_order[0][2].item() << '\n'
              << "autodiff put  charm = " << -put_greeks_2nd_order[0][2].item() << '\n'
              << "      formula charm = " << formula_charm << '\n'
              << "autodiff call vomma = " << call_greeks_2nd_order[1][1].item() << '\n'
              << "autodiff put  vomma = " << put_greeks_2nd_order[1][1].item() << '\n'
              << "      formula vomma = " << formula_vomma << '\n'
              << "autodiff call veta = " << call_greeks_2nd_order[1][2].item() << '\n'
              << "autodiff put  veta = " << put_greeks_2nd_order[1][2].item() << '\n'
              << "      formula veta = " << formula_veta << '\n'

              << "\n## Third-order Greeks\n"
              << "autodiff call speed = " << call_greeks_3rd_order[0][0][0] << '\n'
              << "autodiff put  speed = " << put_greeks_3rd_order[0][0][0] << '\n'
              << "      formula speed = " << formula_speed << '\n'
              << "autodiff call zomma = " << call_greeks_3rd_order[0][0][1] << '\n'
              << "autodiff put  zomma = " << put_greeks_3rd_order[0][0][1] << '\n'
              << "      formula zomma = " << formula_zomma << '\n'
              << "autodiff call color = " << call_greeks_3rd_order[0][0][2] << '\n'
              << "autodiff put  color = " << put_greeks_3rd_order[0][0][2] << '\n'
              << "      formula color = " << formula_color << '\n'
              << "autodiff call ultima = " << call_greeks_3rd_order[1][1][1] << '\n'
              << "autodiff put  ultima = " << put_greeks_3rd_order[1][1][1] << '\n'
              << "      formula ultima = " << formula_ultima << '\n';
    return 0;
}
