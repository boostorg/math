// Copyright Nick Thompson, 2019
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)
#ifndef BOOST_MATH_QUADRATURE_DETAIL_OOURA_FOURIER_INTEGRALS_DETAIL_HPP
#define BOOST_MATH_QUADRATURE_DETAIL_OOURA_FOURIER_INTEGRALS_DETAIL_HPP
#include <utility> // for std::pair.
#include <mutex>
#include <boost/math/special_functions/expm1.hpp>
#include <boost/math/special_functions/sin_pi.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/math/tools/condition_numbers.hpp>

namespace boost { namespace math { namespace quadrature { namespace detail {

// Ooura and Mori, A robust double exponential formula for Fourier-type integrals,
// eta is the argument to the exponential in equation 3.3:
template<class Real>
std::pair<Real, Real> ooura_eta(Real x, Real alpha) {
    using std::expm1;
    using std::exp;
    using std::abs;
    Real expx = exp(x);
    Real eta_prime = 2 + alpha/expx + expx/4;
    Real eta;
    // This is the fast branch:
    if (abs(x) > 0.125) {
        eta = 2*x - alpha*(1/expx - 1) + (expx - 1)/4;
    }
    else {// this is the slow branch using expm1 for small x:
        eta = 2*x - alpha*expm1(-x) + expm1(x)/4;
    }
    return {eta, eta_prime};
}

// Ooura and Mori, A robust double exponential formula for Fourier-type integrals,
// equation 3.6:
template<class Real>
Real calculate_ooura_alpha(Real h)
{
    using boost::math::constants::pi;
    using std::log1p;
    using std::sqrt;
    Real x = sqrt(16 + 4*log1p(pi<Real>()/h)/h);
    return 1/x;
}

template<class Real>
std::pair<Real, Real> ooura_sin_node_and_weight(long n, Real h, Real alpha)
{
    using std::expm1;
    using std::exp;
    using boost::math::constants::pi;

    if (n == 0) {
        // Equation 44 of https://arxiv.org/pdf/0911.4796.pdf
        Real eta_prime_0 = Real(2) + alpha + Real(1)/Real(4);
        Real node = pi<Real>()/(eta_prime_0*h);
        Real weight = pi<Real>()*boost::math::sin_pi(1/(eta_prime_0*h));
        Real eta_dbl_prime = -alpha + Real(1)/Real(4);
        Real phi_prime_0 = (1 - eta_dbl_prime/(eta_prime_0*eta_prime_0))/2;
        weight *= phi_prime_0;
        return {node, weight};
    }
    Real x = n*h;
    auto [eta, eta_prime] = ooura_eta(x, alpha);

    Real expm1_meta = expm1(-eta);
    Real exp_meta = exp(-eta);
    Real node = -n*pi<Real>()/expm1_meta;

    // I do worry about whether this is the most accurate method of computing phi'(x):
    Real phi_prime = -(expm1_meta + x*exp_meta*eta_prime)/(expm1_meta*expm1_meta);


    Real s = pi<Real>();
    Real arg;
    if(eta > 1) {
        arg = n/( 1/exp_meta - 1 );
    }
    else {
        arg = -n*exp_meta/expm1_meta;
    }

    s *= boost::math::sin_pi<Real>(arg);
    if (n&1) {
        s *= -1;
    }

    Real weight = s*phi_prime;
    using std::isfinite;
    if (!isfinite(weight)) {
        throw std::domain_error("Weight is not finite.");
    }
    if (node <= 0) {
        throw std::domain_error("Computed a non-positive quadrature node; this is a major problem.");
    }
    return {node, weight};
}

template<class Real>
class ooura_fourier_sin_detail {
public:
    ooura_fourier_sin_detail(const Real relative_error_tolerance, size_t levels = sizeof(Real)) {
        using std::abs;
        rel_err_ = relative_error_tolerance;
        Real unit_roundoff = std::numeric_limits<Real>::epsilon()/2;
        big_nodes_.resize(levels);
        bweights_.resize(levels);
        little_nodes_.resize(levels);
        lweights_.resize(levels);

        // h0 = 1. Then all further levels have h_i = 1/2^i.
        Real h = 1;
        for (size_t i = 0; i < big_nodes_.size(); ++i) {
            auto& bnode_row = big_nodes_[i];
            auto& bweight_row = bweights_[i];
            bnode_row.reserve((1<<i)*sizeof(Real));
            bweight_row.reserve((1<<i)*sizeof(Real));

            auto& lnode_row = little_nodes_[i];
            auto& lweight_row = lweights_[i];
            lnode_row.reserve((1<<i)*sizeof(Real));
            lweight_row.reserve((1<<i)*sizeof(Real));

            Real max_weight = 1;

            Real alpha = calculate_ooura_alpha(h);
            long n = 0;
            Real w;
            do {
                auto [node, weight] = ooura_sin_node_and_weight(n, h, alpha);
                w = weight;
                bnode_row.push_back(node);
                bweight_row.push_back(weight);
                if (abs(weight) > max_weight) {
                    max_weight = abs(weight);
                }
                ++n;
                // f(t)->0 as t->infty, which is why the weights are computed up to the unit roundoff.
            } while(abs(w) > unit_roundoff*max_weight);

            // It's reasonable to ask why we are splitting the nodes into regimes where t_n >> 1 and t_n << 1.
            // This is because it will create the opportunity to sensibly truncate the quadrature sum to significant terms.
            n = -1;
            do {
                auto [node, weight] = ooura_sin_node_and_weight(n, h, alpha);
                w = weight;
                lnode_row.push_back(node);
                lweight_row.push_back(weight);
                if (abs(weight) > max_weight) {
                    max_weight = abs(weight);
                }
                --n;
                // f(t)->infty is possible as t->0, hence compute up to the min.
            } while(abs(w) > std::numeric_limits<Real>::min()*max_weight);

            // Since the nodes don't nest, we could conceivably divide h by (say) 1.5, or 3.
            // It's not clear how much benefit (or loss) would be obtained from this.
            h /= 2;
        }
    }

    template<class F>
    Real integrate(F const & f, Real omega) {
        using std::abs;
        using std::max;
        using boost::math::constants::pi;

        if (omega == 0) {
            return Real(0);
        }
        if (omega < 0) {
            return -this->integrate(f, -omega);
        }
        //Real I0 = 0;
        Real I1 = std::numeric_limits<Real>::quiet_NaN();
        size_t i = 1;
        Real inv_omega = Real(1)/omega;
        do {
            auto cond = boost::math::tools::summation_condition_number<Real, true>(0);
            size_t calls = 0;
            auto& b_nodes = big_nodes_[i];
            auto& b_weights = bweights_[i];
            for(size_t j = 0 ; j < b_nodes.size(); ++j) {
                //I0 += f(b_nodes[j]*inv_omega)*b_weights[j];
                cond += f(b_nodes[j]*inv_omega)*b_weights[j];
                ++calls;
            }

            auto& l_nodes = little_nodes_[i];
            auto& l_weights = lweights_[i];
            // If f decays rapidly as |t|->infty, not all of these calls are necessary.
            for (size_t j = 0; j < l_nodes.size(); ++j) {
                //I0 += f(l_nodes[j]*inv_omega)*l_weights[j];
                cond += f(l_nodes[j]*inv_omega)*l_weights[j];
                ++calls;
            }

            std::cout << "I0 = " << cond.sum()/omega << ", calls = " << calls << ", err est = " << abs(cond.sum()-I1)  << ",  cond = " << cond() << "\n";
            if (abs(cond.sum()-I1) < rel_err_*max(abs(cond.sum()), abs(I1))) {
                return cond.sum()/omega;
            }
            I1 = cond.sum();
            //I0 = 0;
            ++i;

        } while(i < big_nodes_.size());

        std::cout << "Warning: Used all available levels.\n";
        return I1/omega;
    }

private:
    // Nodes for n >= 0, giving t_n = pi*phi(nh)/h. Generally t_n >> 1.
    std::vector<std::vector<Real>> big_nodes_;
    // The term bweights_ will indicate that these are weights corresponding
    // to the big nodes:
    std::vector<std::vector<Real>> bweights_;

    // Nodes for n < 0: Generally t_n << 1, and an invariant is that t_n > 0.
    std::vector<std::vector<Real>> little_nodes_;
    std::vector<std::vector<Real>> lweights_;
    Real rel_err_;
};

}}}}
#endif
