// Copyright Matt Borland, 2021
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Source: https://www.its.bldrdoc.gov/DownloadPublications.ashx?id=OT-75-70.pdf

#include <execution>
#include <type_traits>
#include <stdexcept>
#include <cstddef>

namespace boost { namespace math { namespace interpolators { namespace detail {

template <typename ExecutionPolicy, typename ForwardContainer>
class akima_impl
{
private:
    ExecutionPolicy exec_;
    const ForwardContainer& x_;
    const ForwardContainer& y_;
    const ForwardContainer& z_;
    const std::size_t n_;

    void validate_inputs();
public:
    akima_impl(ExecutionPolicy&& exec, const ForwardContainer& x, const ForwardContainer& y, const ForwardContainer& z, 
               const std::size_t n) : exec_{exec}, x_{x}, y_{y}, z_{z}, n_{n} 
    {
        validate_inputs();
    }
};

template <typename ExecutionPolicy, typename ForwardContainer>
void akima_impl<ExecutionPolicy, ForwardContainer>::validate_inputs()
{
    const auto x_size = std::size(x_);
    const auto y_size = std::size(y_);
    const auto z_size = std::size(z_);

    if(x_size != y_size || x_size != z_size || y_size != z_size)
    {
        throw std::domain_error("X, Y, and Z must all be of the same size");
    }

    if(x_size < 4)
    {
        throw std::domain_error("X, Y, and Z must all have at least 4 values");
    }

    if(n_ == 0)
    {
        throw std::domain_error("N must be at least 1");
    }
}

}}}} // namespaces
