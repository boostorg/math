/*
 * Copyright Nick Thompson, 2019
 * Use, modification and distribution are subject to the
 * Boost Software License, Version 1.0. (See accompanying file
 * LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 */
#include <iostream>
#include <vector>
#include <numeric>
#include <list>
#include <cmath>
#include <cassert>
#include <fstream>
#include <Eigen/Eigenvalues>
#include <boost/hana/for_each.hpp>
#include <boost/hana/ext/std/integer_sequence.hpp>
#include <boost/core/demangle.hpp>
#include <boost/multiprecision/float128.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/math/filters/daubechies.hpp>
#include <boost/math/special_functions/factorials.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>



template<class Real, int p>
std::list<std::vector<Real>> integer_grid()
{
    std::cout << std::setprecision(std::numeric_limits<Real>::digits10 + 3);
    using std::abs;
    using std::sqrt;
    using std::pow;
    std::list<std::vector<Real>> grids;

    auto c = boost::math::filters::daubechies_scaling_filter<Real, p>();
    for (auto & x : c)
    {
        x *= boost::math::constants::root_two<Real>();
    }
    std::cout << "\n\nTaps in filter = " << c.size() << "\n";


    Eigen::Matrix<Real, 2*p - 2, 2*p-2> A;
    for (int j = 0; j < 2*p-2; ++j) {
        for (int k = 0; k < 2*p-2; ++k) {
            if ( (2*j-k + 1) < 0 || (2*j - k  + 1) >= 2*p)
            {
                A(j,k) = 0;
            }
            else {
                A(j,k) = c[2*j - k + 1];
            }
        }
    }

    Eigen::EigenSolver<decltype(A)> es(A);

    auto complex_eigs = es.eigenvalues();

    std::vector<Real> eigs(complex_eigs.size(), std::numeric_limits<Real>::quiet_NaN());

    std::cout << "Eigenvalues = {";
    for (long i = 0; i < complex_eigs.size(); ++i) {
        assert(abs(complex_eigs[i].imag()) < std::numeric_limits<Real>::epsilon());
        eigs[i] = complex_eigs[i].real();
        std::cout << eigs[i] << ", ";
    }
    std::cout << "}\n";

    // Eigen does not sort the eigenpairs by any criteria on the eigenvalues.
    // In any case, even if it did, some of the eigenpairs do not correspond to derivatives anyway.
    for (size_t j = 0; j < eigs.size(); ++j) {
        auto f = [&](Real x) {
                 return abs(x - Real(1)/Real(1 << j) ) < sqrt(std::numeric_limits<Real>::epsilon());
                 };
        auto it = std::find_if(eigs.begin(), eigs.end(), f);
        if (it == eigs.end()) {
            std::cout << "couldn't find eigenvalue " << Real(1)/Real(1 << j) << "\n";
            continue;
        }
        size_t idx = std::distance(eigs.begin(), it);
        std::cout << "Eigenvector for derivative " << j << " is at index " << idx << "\n";
        auto const & complex_eigenvec = es.eigenvectors().col(idx);
        std::vector<Real> eigenvec(complex_eigenvec.size() + 2, std::numeric_limits<Real>::quiet_NaN());
        eigenvec[0] = 0;
        eigenvec[eigenvec.size()-1] = 0;
        for (size_t i = 0; i < eigenvec.size() - 2; ++i) {
            assert(abs(complex_eigenvec[i].imag()) < std::numeric_limits<Real>::epsilon());
            eigenvec[i+1] = complex_eigenvec[i].real();
        }

        Real sum = 0;
        for(size_t k = 1; k < eigenvec.size(); ++k) {
            sum += pow(k, j)*eigenvec[k];
        }

        Real alpha = pow(-1, j)*boost::math::factorial<Real>(j)/sum;

        for (size_t i = 1; i < eigenvec.size(); ++i) {
            eigenvec[i] *= alpha;
        }


        std::cout << "Eigenvector = {";
        for (size_t i = 0; i < eigenvec.size() -1; ++i) {
            std::cout << eigenvec[i] << ", ";
        }
        std::cout << eigenvec[eigenvec.size()-1] << "}\n";

        sum = 0;
        for(size_t k = 1; k < eigenvec.size(); ++k) {
            sum += pow(k, j)*eigenvec[k];
        }

        std::cout << "Moment sum = " << sum << ", expected = " << pow(-1, j)*boost::math::factorial<Real>(j) << "\n";

        assert(abs(sum - pow(-1, j)*boost::math::factorial<Real>(j))/abs(pow(-1, j)*boost::math::factorial<Real>(j)) < sqrt(std::numeric_limits<Real>::epsilon()));

        grids.push_back(eigenvec);
    }


    return grids;
}

template<class Real, int p>
void write_grid(std::ofstream & fs)
{
    auto grids = integer_grid<Real, p>();
    fs << "    if constexpr (p == " << p << ") {\n";
    fs << "        if constexpr (std::is_same_v<Real, float>) {\n";
    size_t j = 0;
    for (auto it = grids.begin(); it != grids.end(); ++it) {

    fs << "            if constexpr (order == " << j << ") {\n";
    fs << "                return {";
        auto const & grid = *it;
        for (size_t i = 0; i < grid.size() -1; ++i) {
            fs << static_cast<float>(grid[i]) << "f, ";
        }
        fs << static_cast<float>(grid[grid.size()-1]) << "f};\n";
    fs << "            }\n";
        ++j;
    }
    fs << "        }\n";

    fs << "        if constexpr (std::is_same_v<Real, double>) {\n";
    j = 0;
    for (auto it = grids.begin(); it != grids.end(); ++it) {

    fs << "            if constexpr (order == " << j << ") {\n";
    fs << "                return {";
        auto const & grid = *it;
        for (size_t i = 0; i < grid.size() -1; ++i) {
            fs << static_cast<double>(grid[i]) << ", ";
        }
        fs << static_cast<double>(grid[grid.size()-1]) << "};\n";
    fs << "            }\n";
        ++j;
    }
    fs << "        }\n";


    fs << "        if constexpr (std::is_same_v<Real, long double>) {\n";
    j = 0;
    for (auto it = grids.begin(); it != grids.end(); ++it) {

    fs << "            if constexpr (order == " << j << ") {\n";
    fs << "                return {";
        auto const & grid = *it;
        for (size_t i = 0; i < grid.size() -1; ++i) {
            fs << static_cast<long double>(grid[i]) << "L, ";
        }
        fs << static_cast<long double>(grid[grid.size()-1]) << "L};\n";
    fs << "            }\n";
        ++j;
    }
    fs << "        }\n";

    fs << "        #ifdef BOOST_HAS_FLOAT128\n";
    fs << "        if constexpr (std::is_same_v<Real, boost::multiprecision::float128>) {\n";
    j = 0;
    for (auto it = grids.begin(); it != grids.end(); ++it) {

    fs << "            if constexpr (order == " << j << ") {\n";
    fs << "                return {";
        auto const & grid = *it;
        for (size_t i = 0; i < grid.size() -1; ++i) {
            fs << static_cast<boost::multiprecision::float128>(grid[i]) << "Q, ";
        }
        fs << static_cast<boost::multiprecision::float128>(grid[grid.size()-1]) << "Q};\n";
    fs << "            }\n";
        ++j;
    }
    fs << "        }\n";
    fs << "        #endif\n";

    fs << "        if constexpr (std::is_same_v<Real, boost::multiprecision::cpp_bin_float_oct>) {\n";
    j = 0;
    for (auto it = grids.begin(); it != grids.end(); ++it) {

    fs << std::setprecision(std::numeric_limits<boost::multiprecision::cpp_bin_float_oct>::digits10 + 3);
    fs << "            if constexpr (order == " << j << ") {\n";
    fs << "                return {";
        auto const & grid = *it;
        for (size_t i = 0; i < grid.size() -1; ++i) {
            fs << "boost::lexical_cast<boost::multiprecision::cpp_bin_float_oct>(\"" << static_cast<boost::multiprecision::cpp_bin_float_oct>(grid[i]) << "\"), ";
        }
        fs << "boost::lexical_cast<boost::multiprecision::cpp_bin_float_oct>(\"" << static_cast<boost::multiprecision::cpp_bin_float_oct>(grid[grid.size()-1]) << "\")};\n";
    fs << "            }\n";
        ++j;
    }
    fs << "        }\n";

    fs << "    }\n";

}

int main()
{
    constexpr const size_t p_max = 15;
    std::ofstream fs{"daubechies_scaling_integer_grid.hpp"};
    fs << "/*\n"
       << " * Copyright Nick Thompson, 2019\n"
       << " * Use, modification and distribution are subject to the\n"
       << " * Boost Software License, Version 1.0. (See accompanying file\n"
       << " * LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)\n"
       << " */\n"
       << "// THIS FILE GENERATED BY EXAMPLE/DAUBECHIES_SCALING_INTEGER_GRID.CPP, DO NOT EDIT.\n"
       << "#ifndef BOOST_MATH_DAUBECHIES_SCALING_INTEGER_GRID_HPP\n"
       << "#define BOOST_MATH_DAUBECHIES_SCALING_INTEGER_GRID_HPP\n"
       << "#include <array>\n"
       << "#ifdef BOOST_HAS_FLOAT128\n"
       << "#include <boost/multiprecision/float128.hpp>\n"
       << "#endif\n"
       << "namespace boost::math::detail {\n\n"
       << "template <typename Real, unsigned p, unsigned order>\n"
       << "constexpr std::array<Real, 2*p> daubechies_scaling_integer_grid()\n"
       << "{\n"
       << "    static_assert(sizeof(Real) <= 16, \"Integer grids only computed up to 128 bits of precision.\");\n"
       << "    static_assert(p <= " << p_max << ", \"Integer grids only implemented up to " << p_max << ".\");\n"
       << "    static_assert(p > 1, \"Integer grids only implemented for p >= 2.\");\n";


    fs << std::hexfloat;

    boost::hana::for_each(std::make_index_sequence<p_max>(), [&](auto idx){
        write_grid<boost::multiprecision::cpp_bin_float_oct, idx+2>(fs);
    });

    fs << "    std::array<Real, 2*p> m{};\n"
       << "    for (auto & x : m) {\n"
       << "        x = std::numeric_limits<Real>::quiet_NaN();\n"
       << "    }\n"
       << "    return m;\n";


    fs << "}\n\n";

    fs << "} // namespaces\n";
    fs << "#endif\n";
    fs.close();

    return 0;
}
