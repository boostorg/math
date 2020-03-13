/*
 * Copyright Nick Thompson, 2020
 * Use, modification and distribution are subject to the
 * Boost Software License, Version 1.0. (See accompanying file
 * LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 */
#include <cstdint>
#include <boost/math/special_functions/daubechies_wavelet.hpp>
#include <boost/math/quadrature/trapezoidal.hpp>
#include <boost/gil.hpp>
#include <boost/gil/extension/io/png.hpp>
#include <boost/gil/extension/io/jpeg.hpp>
#include <boost/gil/image.hpp>

namespace bg = boost::gil;

template<class F, typename Real, int p>
class daubechies_wavelet_transform
{
public:
    daubechies_wavelet_transform(F f, int grid_refinements = -1, Real tol = boost::math::tools::root_epsilon<Real>(),
    int max_refinements = 12) : f_{f}, psi_(grid_refinements), tol_{tol}, max_refinements_{max_refinements}
    {}

    Real operator()(Real s, Real t) const
    {
        using std::sqrt;
        using boost::math::quadrature::trapezoidal;
        if (s <= 0)
        {
            throw std::domain_error("s > 0 is required.");
        }
        Real a = -s*p + s + t;
        Real b = s*p + t;
        Real Q = trapezoidal(f_, a, b, tol_, max_refinements_);
        return Q/sqrt(s);
    }

private:
    F f_;
    boost::math::daubechies_wavelet<Real, p> psi_;
    Real tol_;
    int max_refinements_;
};



int main(int argc, char *argv[])
{
    int height = 1000;
    int width = 1.618*height;
    auto img = bg::gray16_image_t{width, height, bg::gray16_pixel_t {0}};
    auto view = bg::view(img);
    auto count = std::uint16_t {0};
    for (auto it = view.begin(); it != view.end(); ++it)
    {
        *it = count++;
    }

    //bg::write_view("img.png", bg::const_view(img), bg::png_tag());
    bg::write_view("img.jpeg", bg::const_view(img), bg::png_tag());
    return 0;
}

/*int main()
{
    double a = 1.3;
    auto f = [&a](double t) {
        if(t==0) {
            return double(0);
        }
        return std::sin(a/t);
    };

    auto Wf = daubechies_wavelet_transform<decltype(f), double, 8>(f);

    double s = 7;
    double t = 9.2;
    double Wfst = Wf(s, t);
    std::cout << "W[f](s,t) = " << Wfst << "\n";
    return 0;
}*/