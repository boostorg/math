//  (C) Copyright Nick Thompson 2018.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_TOOLS_NORMS_HPP
#define BOOST_MATH_TOOLS_NORMS_HPP

#include <algorithm>
#include <iterator>
#include <boost/type_traits.hpp>
#include <boost/assert.hpp>
#include <boost/multiprecision/detail/number_base.hpp>


namespace boost{ namespace math{ namespace tools {

// Mallat, "A Wavelet Tour of Signal Processing", equation 2.60:
template<class ForwardIterator>
auto total_variation(ForwardIterator first, ForwardIterator last)
{
    using Real = typename std::iterator_traits<ForwardIterator>::value_type;
    using std::abs;
    BOOST_ASSERT_MSG(first != last && std::next(first) != last, "At least two samples are required to compute the total variation.");
    Real tv = 0;
    auto it = first;
    Real tmp = *it;

    if constexpr (std::is_unsigned<Real>::value)
    {
        while (++it != last)
        {
            if (*it > tmp)
            {
                tv += *it - tmp;
            }
            else
            {
                tv += tmp - *it;
            }
            tmp = *it;
        }
        return tv;
    }
    else
    {
        while (++it != last)
        {
            tv += abs(*it - tmp);
            tmp = *it;
        }
        return tv;
    }
}

template<class Container>
inline auto total_variation(Container const & v)
{
    return total_variation(v.cbegin(), v.cend());
}


template<class ForwardIterator>
auto sup_norm(ForwardIterator first, ForwardIterator last)
{
    BOOST_ASSERT_MSG(first != last, "At least one value is required to compute the sup norm.");
    using RealOrComplex = typename std::iterator_traits<ForwardIterator>::value_type;
    using std::abs;
    if constexpr (boost::is_complex<RealOrComplex>::value ||
                  boost::multiprecision::number_category<RealOrComplex>::value == boost::multiprecision::number_kind_complex)
    {
        auto it = std::max_element(first, last, [](RealOrComplex a, RealOrComplex b) { return abs(b) > abs(a); });
        return abs(*it);
    }
    else
    {
        auto pair = std::minmax_element(first, last);
        if (abs(*pair.first) > abs(*pair.second))
        {
            return abs(*pair.first);
        }
        else
        {
            return abs(*pair.second);
        }
    }
}

template<class Container>
inline auto sup_norm(Container const & v)
{
    return sup_norm(v.cbegin(), v.cend());
}

template<class ForwardIterator>
auto l1_norm(ForwardIterator first, ForwardIterator last)
{
    using std::abs;
    decltype(abs(*first)) l1 = 0;
    for (auto it = first; it != last; ++it)
    {
        l1 += abs(*first);
    }
    return l1;
}

template<class Container>
inline auto l1_norm(Container const & v)
{
    return l1_norm(v.cbegin(), v.cend());
}


template<class ForwardIterator>
auto l2_norm(ForwardIterator first, ForwardIterator last)
{
    using RealOrComplex = typename std::iterator_traits<ForwardIterator>::value_type;
    using std::abs;
    using std::norm;
    using std::sqrt;
    using std::is_floating_point;
    if constexpr (boost::is_complex<RealOrComplex>::value ||
                  boost::multiprecision::number_category<RealOrComplex>::value == boost::multiprecision::number_kind_complex)
    {
        typedef typename RealOrComplex::value_type Real;
        Real l2 = 0;
        for (auto it = first; it != last; ++it)
        {
            l2 += norm(*it);
        }
        Real result = sqrt(l2);
        if (!isfinite(result))
        {
            Real a = sup_norm(first, last);
            l2 = 0;
            for (auto it = first; it != last; ++it)
            {
                l2 += norm(*it/a);
            }
            return a*sqrt(l2);
        }
        return result;
    }
    else if constexpr (is_floating_point<RealOrComplex>::value ||
                       boost::multiprecision::number_category<RealOrComplex>::value == boost::multiprecision::number_kind_floating_point)
    {
        RealOrComplex l2 = 0;
        for (auto it = first; it != last; ++it)
        {
            l2 += (*it)*(*it);
        }
        RealOrComplex result = sqrt(l2);
        if (!isfinite(result))
        {
            RealOrComplex a = sup_norm(first, last);
            l2 = 0;
            for (auto it = first; it != last; ++it)
            {
                RealOrComplex tmp = *it/a;
                l2 += tmp*tmp;
            }
            return a*sqrt(l2);
        }
        return result;
    }
}

template<class Container>
inline auto l2_norm(Container const & v)
{
    return l2_norm(v.cbegin(), v.cend());
}

template<class ForwardIterator>
size_t l0_pseudo_norm(ForwardIterator first, ForwardIterator last)
{
    using RealOrComplex = typename std::iterator_traits<ForwardIterator>::value_type;
    size_t count = 0;
    for (auto it = first; it != last; ++it)
    {
        if (*it != RealOrComplex(0))
        {
            ++count;
        }
    }
    return count;
}

template<class Container>
inline size_t l0_pseudo_norm(Container const & v)
{
    return l0_pseudo_norm(v.cbegin(), v.cend());
}

template<class ForwardIterator>
auto lp_norm(ForwardIterator first, ForwardIterator last, typename std::remove_const<typename std::remove_reference<decltype(*std::declval<ForwardIterator>())>::type>::type p)
{
    using std::pow;
    using std::is_floating_point;
    using std::isfinite;
    using RealOrComplex = typename std::iterator_traits<ForwardIterator>::value_type;
    if constexpr (boost::is_complex<RealOrComplex>::value ||
                  boost::multiprecision::number_category<RealOrComplex>::value == boost::multiprecision::number_kind_complex)
    {
        BOOST_ASSERT_MSG(p.real() >= 0, "For p < 0, the lp norm is not a norm.");
        BOOST_ASSERT_MSG(p.imag() == 0, "For imaginary p, the lp norm is not a norm.");
        using std::norm;
        decltype(p.real()) lp = 0;
        for (auto it = first; it != last; ++it)
        {
            lp += pow(norm(*it), p.real()/2);
        }

        auto result = pow(lp, 1/p.real());
        if (!isfinite(result))
        {
            auto a = boost::math::tools::sup_norm(first, last);
            decltype(p.real()) lp = 0;
            for (auto it = first; it != last; ++it)
            {
                lp += pow(abs(*it)/a, p.real());
            }
            result = a*pow(lp, 1/p.real());
        }
        return result;
    }
    else if constexpr (is_floating_point<RealOrComplex>::value ||
                       boost::multiprecision::number_category<RealOrComplex>::value == boost::multiprecision::number_kind_floating_point)
    {
        BOOST_ASSERT_MSG(p >= 0, "For p < 0, the lp norm is not a norm");
        RealOrComplex lp = 0;

        for (auto it = first; it != last; ++it)
        {
            lp += pow(abs(*it), p);
        }

        RealOrComplex result = pow(lp, 1/p);
        if (!isfinite(result))
        {
            RealOrComplex a = boost::math::tools::sup_norm(first, last);
            lp = 0;
            for (auto it = first; it != last; ++it)
            {
                lp += pow(abs(*it)/a, p);
            }
            result = a*pow(lp, 1/p);
        }
        return result;
    }
    else
    {
        BOOST_ASSERT_MSG(false, "Unable to determine if the input type is real or complex.");
    }
}

template<class Container>
inline auto lp_norm(Container const & v, typename Container::value_type p)
{
    return lp_norm(v.cbegin(), v.cend(), p);
}

}}}
#endif
