#ifndef BOOST_MATH_DISTRIBUTIONS_ANDERSON_DARLING_HPP
#define BOOST_MATH_DISTRIBUTIONS_ANDERSON_DARLING_HPP

#include <cmath>
#include <algorithm>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/tools/univariate_statistics.hpp>
#include <boost/math/quadrature/exp_sinh.hpp>
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <boost/math/special_functions/erf.hpp>

namespace boost { namespace math {

template<class RandomAccessContainer>
auto anderson_darling_normality_test(RandomAccessContainer const & v,
                                     typename RandomAccessContainer::value_type mu = std::numeric_limits<typename RandomAccessContainer::value_type>::quiet_NaN(),
                                     typename RandomAccessContainer::value_type sd = std::numeric_limits<typename RandomAccessContainer::value_type>::quiet_NaN())
{
    using Real = typename RandomAccessContainer::value_type;
    using std::log;
    using std::pow;
    using std::sqrt;
    using boost::math::erfc;
    using boost::math::erf;

    if (std::isnan(mu)) {
        mu = boost::math::tools::mean(v);
    }
    if (std::isnan(sd)) {
        sd = sqrt(boost::math::tools::sample_variance(v));
    }

    auto normal = boost::math::normal_distribution(mu, sd);

    // This is where Knuth's literate programming could really come in handy!
    // I need some LaTeX. The idea is that before any observation, the ecdf is identically zero.
    // So we need to compute:
    // \int_{0}^{v_0} \frac{F(x)F'(x)}{1- F(x)} \, \mathrm{d}x, where F(x) := \frac{1}{2}[1+\erf(\frac{x-\mu}{\sigma \sqrt{2}})]
    // Astonishingly, there is an analytic evaluation to this integral, as you can validate with the following Mathematica command:
    // Integrate[(1/2 (1 + Erf[(x - mu)/Sqrt[2*sigma^2]])*Exp[-(x - mu)^2/(2*sigma^2)]*1/Sqrt[2*\[Pi]*sigma^2])/(1 - 1/2 (1 + Erf[(x - mu)/Sqrt[2*sigma^2]])),
    // {x, -Infinity, x0}, Assumptions -> {x0 \[Element] Reals && mu \[Element] Reals && sigma > 0}]

    Real z0 = (v[0] - mu)/sd;
    Real erfcz0 = erfc(z0/sqrt(Real(2)));
    Real left_tail = -1 + erfcz0/2 + log(2/erfcz0);

    // For the right tail, the ecdf is identically 1.
    // Hence we need the integral:
    // \int_{v_{n-1}}^{\infty} \frac{(1-F(x))F'(x)}{F(x)} \, \mathrm{d}x
    // This also has an analytic evaluation! It can be found via the following Mathematica command:
    // Integrate[(E^(-(z^2/2)) *(1 - 1/2 (1 + Erf[z/Sqrt[2]])))/(Sqrt[2 \[Pi]] (1/2 (1 + Erf[z/Sqrt[2]]))),
    // {z, zn, \[Infinity]}, Assumptions -> {zn \[Element] Reals && mu \[Element] Reals}]

    Real zf = (v[v.size()-1] - mu)/sd;
    Real erfzf = erf(zf/sqrt(Real(2)));
    Real right_tail = (-1 + erfzf + 2*log(2/(1+erfzf)))/2;

    // Now we need each integral:
    // \int_{v_i}^{v_{i+1}} \frac{(i+1/n - F(x))^2F'(x)}{F(x)(1-F(x))}  \, \mathrm{d}x
    // Again we get an analytical evaluation via the following Mathematica command:
    // Integrate[((E^(-(z^2/2))/Sqrt[2 \[Pi]])*(k1 - F[z])^2)/(F[z]*(1 - F[z])),
    // {z, z1, z2}, Assumptions -> {z1 \[Element] Reals && z2 \[Element] Reals &&k1 \[Element] Reals}] // FullSimplify

    Real integrals = 0;
    int64_t N = v.size();
    Real s0 = (v[0]-mu)/(sd*sqrt(Real(2)));
    Real erfcs0 = erfc(s0);
    for (int64_t i = 0; i < N - 1; ++i) {
        Real k = (i+1)/Real(N);
        Real s1 = (v[i+1]-mu)/(sd*sqrt(Real(2)));
        Real erfcs1 = erfc(s1);

        Real term = -erfcs0 + erfcs1 + 2*log(erfcs0/erfcs1) + 2*k*(k*log(erfcs0*(-2 + erfcs1)/(erfcs1*(-2 + erfcs0))) + 2*log(erfcs1/erfcs0));
        term /= 2;
        if (v[i] > v[i+1]) {
            throw std::domain_error("Input data must be sorted in increasing order v[0] <= v[1] <= . . .  <= v[n-1]");
        }
        integrals += term;
        s0 = s1;
        erfcs0 = erfcs1;
    }

    return v.size()*(left_tail + right_tail + integrals);
}

template<class RandomAccessContainer>
auto anderson_darling_normality_step(RandomAccessContainer const & v, typename RandomAccessContainer::value_type mu = 0, typename RandomAccessContainer::value_type sd = 1)
{
    using Real = typename RandomAccessContainer::value_type;
    using std::log;
    using std::pow;

    auto normal = boost::math::normal_distribution(mu, sd);

    auto left_integrand = [&normal](Real x)->Real {
        Real Fx = boost::math::cdf(normal, x);
        Real dmu = boost::math::pdf(normal, x);
        return Fx*dmu/(1-Fx);
    };
    auto es = boost::math::quadrature::exp_sinh<Real>();
    Real left_tail = es.integrate(left_integrand, -std::numeric_limits<Real>::infinity(), v[0]);

    auto right_integrand = [&normal](Real x)->Real {
        Real Fx = boost::math::cdf(normal, x);
        Real dmu = boost::math::pdf(normal, x);
        return (1-Fx)*dmu/Fx;
    };
    Real right_tail = es.integrate(right_integrand, v[v.size()-1], std::numeric_limits<Real>::infinity());


    auto integrator = boost::math::quadrature::gauss<Real, 7>();
    Real integrals = 0;
    int64_t N = v.size();
    for (int64_t i = 0; i < N - 1; ++i) {
        auto integrand = [&normal, &i, &N](Real x)->Real {
            Real Fx = boost::math::cdf(normal, x);
            Real dmu = boost::math::pdf(normal, x);
            Real Fn = (i+1)/Real(N);
            return (Fn - Fx)*(Fn-Fx)*dmu/(Fx*(1-Fx));
        };
        auto term = integrator.integrate(integrand, v[i], v[i+1]);
        if (v[i] > v[i+1]) {
            throw std::domain_error("Input data must be sorted in increasing order v[0] <= v[1] <= . . .  <= v[n-1]");
        }
        integrals += term;
    }

    //std::cout << "Step integrals = " << integrals << "\n";
    //std::cout << "Tails = " << left_tail + right_tail << "\n";

    return v.size()*(left_tail + right_tail + integrals);
}

}}
#endif
