//           Copyright Maksym Zhelyenzyakov 2026.
// Distributed under the Boost Software License, Version 1.0.
//      (See accompanying file LICENSE_1_0.txt or copy at
//           https://www.boost.org/LICENSE_1_0.txt)
#include <boost/math/differentiation/autodiff_reverse.hpp>
#include <boost/math/optimization/gradient_descent.hpp>
#include <boost/math/optimization/minimizer.hpp>
#include <cmath>
#include <fstream>
#include <iostream>
#include <random>
#include <string>
namespace rdiff = boost::math::differentiation::reverse_mode;
namespace bopt  = boost::math::optimization;
double random_double(double min = 0.0, double max = 1.0)
{
    static thread_local std::mt19937       rng{std::random_device{}()};
    std::uniform_real_distribution<double> dist(min, max);
    return dist(rng);
}

template<typename S>
struct vec3
{
    /**
     * @brief R^3 coordinates of particle on Thomson Sphere
     */
    S x, y, z;
};

template<class S>
static inline vec3<S> sph_to_xyz(const S& theta, const S& phi)
{
    /**
     * convenience overload to convert from [theta,phi] -> x, y, z
     */
    return {sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta)};
}

template<typename T>
T thomson_energy(std::vector<T>& r)
{
    /* inverse square law
     */
    const size_t N    = r.size() / 2;
    const T      tiny = T(1e-12);

    T E = 0;
    for (size_t i = 0; i < N; ++i) {
        const T& theta_i = r[2 * i + 0];
        const T& phi_i   = r[2 * i + 1];
        auto     ri      = sph_to_xyz(theta_i, phi_i);

        for (size_t j = i + 1; j < N; ++j) {
            const T& theta_j = r[2 * j + 0];
            const T& phi_j   = r[2 * j + 1];
            auto     rj      = sph_to_xyz(theta_j, phi_j);

            T dx = ri.x - rj.x;
            T dy = ri.y - rj.y;
            T dz = ri.z - rj.z;

            T d2 = dx * dx + dy * dy + dz * dz + tiny;
            E += 1.0 / sqrt(d2);
        }
    }
    return E;
}

template<class T>
std::vector<rdiff::rvar<T, 1>> init_theta_phi_uniform(size_t N, unsigned seed = 12345)
{
    const T pi = T(3.1415926535897932384626433832795);

    std::mt19937                      rng(seed);
    std::uniform_real_distribution<T> unif01(T(0), T(1));
    std::uniform_real_distribution<T> unifm11(T(-1), T(1));

    std::vector<rdiff::rvar<T, 1>> u;
    u.reserve(2 * N);

    for (size_t i = 0; i < N; ++i) {
        T z     = unifm11(rng);
        T phi   = (T(2) * pi) * unif01(rng) - pi;
        T theta = std::acos(z);

        u.emplace_back(theta);
        u.emplace_back(phi);
    }
    return u;
}

int main(int argc, char* argv[])
{
    /*
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <N>\n";
        return 1;
    }*/

    const int    N      = 2;// std::stoi(argv[1]);
    const int    NSTEPS = 100000;
    const double lr     = 1e-3;

    auto u_ad = init_theta_phi_uniform<double>(N);

    auto gdopt = bopt::make_gradient_descent(&thomson_energy<rdiff::rvar<double, 1>>, u_ad, lr);

    // filenames
    std::string pos_filename    = "thomson_" + std::to_string(N) + ".csv";
    std::string energy_filename = "energy_" + std::to_string(N) + ".csv";

    std::ofstream pos_out(pos_filename);
    std::ofstream energy_out(energy_filename);

    pos_out << "step,particle,x,y,z\n";
    energy_out << "step,energy\n";

    /*
    for (int step = 0; step < NSTEPS; ++step) {
        gdopt.step();
        for (int pi = 0; pi < N; ++pi) {
            double theta = u_ad[2 * pi + 0].item();
            double phi   = u_ad[2 * pi + 1].item();
            auto   r     = sph_to_xyz(theta, phi);
            pos_out << step << "," << pi << "," << r.x << "," << r.y << "," << r.z << "\n";
        }
        auto E = gdopt.objective_value();
        energy_out << step << "," << E << "\n";
    }*/

    auto result = minimize(gdopt);
    for (int pi = 0; pi < N; ++pi) {
        double theta = u_ad[2 * pi + 0].item();
        double phi   = u_ad[2 * pi + 1].item();
        auto   r     = sph_to_xyz(theta, phi);
        pos_out << pi << "," << r.x << "," << r.y << "," << r.z << "\n";
    }
    auto E = gdopt.objective_value();
    energy_out << "," << E << "\n";

    pos_out.close();
    energy_out.close();

    return 0;
}
