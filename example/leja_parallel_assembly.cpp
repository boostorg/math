// Copyright Nick Thompson 2026.
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0. (See accompanying file
// LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

// Parallel assembly of a Hamiltonian with one shared Leja rule.
//
// The anharmonic oscillator H = -d^2/dx^2 + x^4, in the basis of Hermite functions psi_n(x / s), n < N. With y = x / s,
// H_ij = int psi_i(y) ((2j + 1 - y^2) / s^2 + s^4 y^4) psi_j(y) dy, since -psi_j'' = (2j + 1 - y^2) psi_j. The integrands
// are entire and negligible beyond |y| = L, and the threads integrate them over [-L, L] with one leja_legendre: whichever
// first needs more points computes them, and the others use them. Each entry takes as many points as it needs, and gets
// an error estimate from the points it has already used.
//
// The matrix elements are also known exactly, from x = (a + a^+) / sqrt(2), and the lowest eigenvalues agree with those
// of Hioe and Montroll (1975).
//
// A rule for exp(-y^2), leja_hermite, would be the wrong choice: h_i h_j grows like exp(y^2), the Leja points spread
// further than the Gauss points, and sum |w_k h_i h_j(x_k)| is 1e19 times H_ij for i = j = 40.

#include <boost/math/quadrature/leja.hpp>
#include <Eigen/Dense>
#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <thread>
#include <vector>

// psi_i(y) psi_j(y), from h_{n+1} = y sqrt(2 / (n + 1)) h_n - sqrt(n / (n + 1)) h_{n-1} times exp(-y^2 / 2).
double psi_product(unsigned i, unsigned j, double y)
{
   double previous = 0;
   double psi = std::exp(-y * y / 2) / std::sqrt(std::sqrt(boost::math::constants::pi<double>()));
   double psi_i = i == 0 ? psi : 0;
   double psi_j = j == 0 ? psi : 0;
   for (unsigned n = 0; n < (std::max)(i, j); ++n)
   {
      const double next = y * std::sqrt(2.0 / (n + 1)) * psi - std::sqrt(double(n) / (n + 1)) * previous;
      previous = psi;
      psi = next;
      if (n + 1 == i)
         psi_i = psi;
      if (n + 1 == j)
         psi_j = psi;
   }
   return psi_i * psi_j;
}

// H from the matrix X of y, tridiagonal with X_{n,n+1} = sqrt((n + 1) / 2): X^2 and X^4 are exact in the first N rows
// when X has N + 4.
Eigen::MatrixXd exact(unsigned N, double s)
{
   Eigen::MatrixXd X = Eigen::MatrixXd::Zero(N + 4, N + 4);
   for (unsigned n = 0; n < N + 3; ++n)
      X(n, n + 1) = X(n + 1, n) = std::sqrt((n + 1) / 2.0);
   const Eigen::MatrixXd X2 = X * X;
   const Eigen::MatrixXd D = Eigen::VectorXd::LinSpaced(N, 1, 2 * N - 1).asDiagonal();
   return (D - X2.topLeftCorner(N, N)) / (s * s) + s * s * s * s * (X2 * X2).topLeftCorner(N, N);
}

int main()
{
   const unsigned N = 60;
   const double s = 0.75;
   const double L = 14;
   const unsigned threads = (std::max)(1u, std::thread::hardware_concurrency());
   const boost::math::quadrature::leja_legendre<double> rule;

   Eigen::MatrixXd H(N, N);
   Eigen::MatrixXi calls(N, N);
   // The entries i <= j, taken in turn by the threads.
   auto assemble = [&] {
      std::atomic<unsigned> next{0};
      auto work = [&] {
         for (unsigned e = next++; e < N * N; e = next++)
         {
            const unsigned i = e / N, j = e % N;
            if (i > j)
               continue;
            unsigned count = 0;
            H(i, j) = H(j, i) = rule.integrate([&](double y) {
               ++count;
               return psi_product(i, j, y) * ((2 * j + 1 - y * y) / (s * s) + s * s * s * s * y * y * y * y);
            }, -L, L, 1e-13);
            calls(i, j) = calls(j, i) = count;
         }
      };
      const auto start = std::chrono::steady_clock::now();
      std::vector<std::thread> pool;
      for (unsigned t = 0; t < threads; ++t)
         pool.emplace_back(work);
      for (std::thread& t : pool)
         t.join();
      return std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
   };
   const double cold = assemble();
   const double warm = assemble();

   const Eigen::MatrixXd H_exact = exact(N, s);
   const double worst = (H - H_exact).cwiseAbs().maxCoeff(), largest = H_exact.cwiseAbs().maxCoeff();
   std::cout << std::setprecision(17) << N * (N + 1) / 2 << " entries on " << threads << " threads, with "
             << calls.minCoeff() << " to " << calls.maxCoeff()
             << " points: " << cold << " s computing the rule as it goes, " << warm << " s with it computed.\n"
             << "max |H - H_exact| = " << worst << ", " << worst / largest << " of the largest entry.\n";
   const Eigen::VectorXd E = Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd>(H, Eigen::EigenvaluesOnly).eigenvalues();
   const double reference[] = {1.0603620904841829, 3.7996730298013941, 7.4556979379867383, 11.644745511378162};
   for (unsigned n = 0; n < 4; ++n)
      std::cout << "E_" << n << " = " << E(n) << ", Hioe and Montroll " << reference[n] << '\n';
}
