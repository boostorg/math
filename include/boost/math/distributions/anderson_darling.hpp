#ifndef BOOST_MATH_DISTRIBUTIONS_ANDERSON_DARLING_HPP
#define BOOST_MATH_DISTRIBUTIONS_ANDERSON_DARLING_HPP

#include <cmath>
#include <algorithm>

namespace boost { namespace math { namespace distributions {

template<class RandomAccessContainer, class CDF>
auto anderson_darling(RandomAccessContainer const & v, CDF& F) {
    using Real = typename RandomAccessContainer::value_type;
    using std::log;
    using std::pow;
    if (!std::is_sorted(v.begin(), v.end())) {
        throw std::domain_error("The input vector must be sorted in non-decreasing order v[0] <= v[1] <= ... <= v[n-1].");
    }


    Real S = 0;
    int64_t N = v.size();
    for (int64_t i = 0; i < N; ++i) {
        Real arg = F(v[i]);
        if (arg < 0 || arg > 1) {
            throw std::domain_error("CDF does not meet requirements of cumulative density function; range should be [0,1] but evaluated to " + std::to_string(arg));
        }
        S += (2*i + 1)*log(F(v[i]))+ (2*v.size()+1-2*(i+1))*log(1-F(v[i]));
        //S += log(pow(arg, 2*i+1)) + log(pow(1-arg, 2*(N-i) -1));
        /*Real arg = F(v[i])*(1-F(v[v.size()-1-i]));
        if(arg <= 0) {
            if (F(v[i]) < 0) {
                throw std::domain_error("The cumulative distribution function passed to the Anderson-Darling test returned a negative value!");
            }

            if (F(v[i]) > 1) {
                throw std::domain_error("The cumulative distribution function passed to the Anderson-Darling test returned a value greater than one!");
            }

            if (F(v[v.size()-1-i]) < 0) {
                throw std::domain_error("The cumulative distribution function passed to the Anderson-Darling test returned a negative value!");
            }

            if (F(v[v.size()-1-i]) > 1) {
                throw std::domain_error("The cumulative distribution function passed to the Anderson-Darling test returned a value greater than one!");
            }

            throw std::domain_error("Argument to the logarithm of the Anderson-Darling test is negative; does the CDF satisfy all the properties of a cumulative distribution function?");
        }
        // TODO: Check if S *= pow(arg, 2*i+1) and then take a log at the end is faster.
        S += (2*i+1)*log(arg);*/
    }
    S /= v.size();
    return -static_cast<Real>(v.size()) - S;
}

}}}
#endif
