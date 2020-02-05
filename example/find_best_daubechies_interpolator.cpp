#include <iostream>
#include <unordered_map>
#include <string>
#include <future>
#include <thread>
#include <boost/math/special_functions/daubechies_scaling.hpp>
#include <boost/math/special_functions/detail/daubechies_scaling_integer_grid.hpp>
#include <boost/math/interpolators/quintic_hermite.hpp>
#include <boost/math/interpolators/cubic_hermite.hpp>
#include <boost/math/interpolators/cardinal_quadratic_b_spline.hpp>
#include <boost/math/interpolators/cardinal_cubic_b_spline.hpp>
#include <boost/math/interpolators/cardinal_quintic_b_spline.hpp>
#include <boost/math/interpolators/whittaker_shannon.hpp>
#include <boost/math/interpolators/cardinal_trigonometric.hpp>
#include <boost/math/special_functions/next.hpp>
#include <boost/math/interpolators/makima.hpp>
#include <boost/math/interpolators/pchip.hpp>
#include <boost/multiprecision/float128.hpp>
#include <boost/core/demangle.hpp>
//#include <quicksvg/graph_fn.hpp>
//#include <quicksvg/ulp_plot.hpp>

using boost::multiprecision::float128;

/*template<class Real, int p>
void do_ulp()
{
    std::cout << "Creating ULP plot on type " << boost::core::demangle(typeid(Real).name()) << " and " << p << " vanishing moments.\n";  
    using std::floor;
    using std::ceil;
    using std::abs;
    int rmax = 14;
    std::cout << "Computing phi_dense\n";
    std::future<std::vector<long double>> f1 = std::async(std::launch::async, [&]{ return boost::math::detail::dyadic_grid<long double, p, 0>(rmax); });
    std::future<std::vector<long double>> f2 = std::async(std::launch::async, [&]{ return boost::math::detail::dyadic_grid<long double, p, 1>(rmax); });

    std::cout << "Computing phi and phi_prime\n";
    std::future<std::vector<long double>> f3 = std::async(std::launch::async, [&]{ return boost::math::detail::dyadic_grid<long double, p, 0>(rmax-2); });
    std::future<std::vector<long double>> f4 = std::async(std::launch::async, [&]{ return boost::math::detail::dyadic_grid<long double, p, 1>(rmax-2); });
    auto phi_dense = f1.get();
    auto phi_dense_prime = f2.get();
    auto phi_accurate = f3.get();
    auto phi_prime_accurate = f4.get();

    Real dx_dense = (2*p-1)/static_cast<Real>(phi_dense.size()-1);
    std::cout << "Done precomputing grids; downcasting now.\n";
    std::vector<Real> phi(phi_accurate.size());
    for (size_t i = 0; i < phi_accurate.size(); ++i) {
        phi[i] = Real(phi_accurate[i]);
    }
    std::vector<Real> phi_prime(phi_accurate.size());
    for (size_t i = 0; i < phi_prime_accurate.size(); ++i) {
        phi_prime[i] = Real(phi_prime_accurate[i]);
    }
    phi_accurate.resize(0);
    phi_prime_accurate.resize(0);

    std::vector<Real> x(phi.size());
    Real dx = (2*p-1)/static_cast<Real>(x.size()-1);
    std::cout << "dx = " << dx << "\n";
    for (size_t i = 0; i < x.size(); ++i) {
      x[i] = i*dx;
    }

    auto ch = boost::math::interpolators::cubic_hermite(std::move(x), std::move(phi), std::move(phi_prime));

    std::vector<long double> x_acc(phi_dense.size());
    
    for (size_t i = 0; i < x_acc.size(); ++i) {
      x_acc[i] = i*dx_dense;
    }

    auto acc = boost::math::interpolators::cubic_hermite(std::move(x_acc), std::move(phi_dense), std::move(phi_dense_prime));
    std::cout << "Writing ulp plot\n";
    std::string title = "daub" + std::to_string(p) + "_" + std::to_string(rmax-2) + "_" + boost::core::demangle(typeid(Real).name()) + ".svg";
    quicksvg::ulp_plot(ch, acc, Real(0), Real(2*p-1), 
                        "ULP plot of Daubechies with " + std::to_string(rmax-2) + " refinements on type " + boost::core::demangle(typeid(Real).name()),
                        title, 15000, 1100, 10);
    std::cout << "Done writing ulp plot\n";

}*/


template<class Real, int p>
void choose_refinement()
{
    using std::abs;
    int rmax = 22;
    auto phi_dense = boost::math::detail::dyadic_grid<long double, p, 0>(rmax);
    Real dx_dense = (2*p-1)/static_cast<Real>(phi_dense.size()-1);

    for (int r = 2; r <= rmax - 2; ++r) {
        auto phi_accurate = boost::math::detail::dyadic_grid<long double, p, 0>(r);
        std::vector<Real> phi(phi_accurate.size());
        for (size_t i = 0; i < phi_accurate.size(); ++i) {
            phi[i] = Real(phi_accurate[i]);
        }
        auto phi_prime_accurate = boost::math::detail::dyadic_grid<long double, p, 1>(r);
        std::vector<Real> phi_prime(phi_accurate.size());
        for (size_t i = 0; i < phi_prime_accurate.size(); ++i) {
            phi_prime[i] = Real(phi_prime_accurate[i]);
        }

        std::vector<Real> x(phi.size());
        Real dx = (2*p-1)/static_cast<Real>(x.size()-1);
        std::cout << "dx = " << dx << "\n";
        for (size_t i = 0; i < x.size(); ++i) {
          x[i] = i*dx;
        }

        if constexpr (p < 6 && p >= 3) {
            auto ch = boost::math::interpolators::cubic_hermite(std::move(x), std::move(phi), std::move(phi_prime));
            Real flt_distance = 0;
            Real sup  = 0;
            Real worst_abscissa = 0;
            Real worst_value = 0;
            Real worst_computed = 0;
            for (size_t i = 0; i < phi_dense.size(); ++i) {
              Real t = i*dx_dense;
              Real computed = ch(t);
              Real expected = Real(phi_dense[i]);
              if (std::abs(expected) < 100*std::numeric_limits<Real>::epsilon()) {
                continue;
              }
              
              Real diff = abs(computed - expected);
              Real distance = abs(boost::math::float_distance(computed, expected));
              if (distance > flt_distance) {
                flt_distance = distance;
                worst_abscissa = t;
                worst_value = expected;
                worst_computed = computed;
              }
              if (diff > sup) {
                sup = diff;
              }
            }
            std::cout << "Float distance at r = " << r << " is " << flt_distance << ", sup distance = " << sup << "\n";
            std::cout << "\tWorst abscissa = " << worst_abscissa << ", worst value = " << worst_value << ", computed = " << worst_computed << "\n"; 
            std::cout << "\tRAM = " << 3*phi_accurate.size()*sizeof(Real) << " bytes\n";
        }
        else if constexpr (p >= 6) {

            auto phi_dbl_prime = boost::math::detail::dyadic_grid<Real, p, 2>(r);
            auto qh = boost::math::interpolators::quintic_hermite(std::move(x), std::move(phi), std::move(phi_prime), std::move(phi_dbl_prime));
            Real flt_distance = 0;
            Real sup  = 0;
            Real worst_abscissa = 0;
            Real worst_value = 0;
            Real worst_computed = 0;
            for (size_t i = 0; i < phi_dense.size(); ++i) {
              Real t = i*dx_dense;
              Real computed = qh(t);
              Real expected = Real(phi_dense[i]);
              if (std::abs(expected) < 100*std::numeric_limits<Real>::epsilon()) {
                continue;
              }
              
              Real diff = abs(computed - expected);
              Real distance = abs(boost::math::float_distance(computed, expected));
              if (distance > flt_distance) {
                flt_distance = distance;
                worst_abscissa = t;
                worst_value = expected;
                worst_computed = computed;
              }
              if (diff > sup) {
                sup = diff;
              }
            }
            std::cout << "Float distance at r = " << r << " is " << flt_distance << ", sup distance = " << sup << "\n";
            std::cout << "\tWorst abscissa = " << worst_abscissa << ", worst value = " << worst_value << ", computed = " << worst_computed << "\n"; 
            std::cout << "\tRAM = " << 3*phi_accurate.size()*sizeof(Real) << " bytes\n";
        }
    }
}

template<class Real, int p>
void find_best_interpolator()
{
    using std::abs;
    int rmax = 15;
    auto phi_dense = boost::math::detail::dyadic_grid<Real, p, 0>(rmax);
    Real dx_dense = (2*p-1)/static_cast<Real>(phi_dense.size()-1);
    for (int r = 2; r < rmax-1; ++r)
    {
        std::map<Real, std::string> m;
        auto phi = boost::math::detail::dyadic_grid<Real, p, 0>(r);
        auto phi_prime = boost::math::detail::dyadic_grid<Real, p, 1>(r);

        std::vector<Real> x(phi.size());
        Real dx = (2*p-1)/static_cast<Real>(x.size()-1);
        std::cout << "dx = 1/" << (1 << r) << " = " << dx << "\n";
        for (size_t i = 0; i < x.size(); ++i) {
            x[i] = i*dx;
        }

        {
          auto phi_copy = phi;
          auto phi_prime_copy = phi_prime;
          auto mh = boost::math::detail::matched_holder(std::move(phi_copy), std::move(phi_prime_copy), r);
          Real sup = 0;
          for (size_t i = 0; i < phi_dense.size(); ++i) {
            Real x = i*dx_dense;
            Real diff = abs(phi_dense[i] - mh(x));
            if (diff > sup) {
              sup = diff;
            }
          }
          m.insert({sup, "matched_holder"});
        }


        {
            auto linear = [&phi, &dx, &r](Real x)->Real {
              if (x <= 0 || x >= 2*p-1) {
                return Real(0);
              }
              using std::floor;

              Real y = (1<<r)*x;
              Real k = floor(y);

              size_t kk = static_cast<size_t>(k);

              Real t = y - k;
              return (1-t)*phi[kk] + t*phi[kk+1];
            };
            Real sup = 0;
            for (size_t i = 0; i < phi_dense.size(); ++i) {
              Real x = i*dx_dense;
              Real diff = abs(phi_dense[i] - linear(x));
              if (diff > sup) {
                sup = diff;
              }
            }
        
          m.insert({sup, "linear interpolation"});
        }
        

        {
            auto qbs = boost::math::interpolators::cardinal_quadratic_b_spline(phi.data(), phi.size(), Real(0), dx, phi_prime.front(), phi_prime.back());
            Real qbs_sup = 0;
            for (size_t i = 0; i < phi_dense.size(); ++i) {
              Real x = i*dx_dense;
              Real diff = abs(phi_dense[i] - qbs(x));
              if (diff > qbs_sup) {
                qbs_sup = diff;
              }
          }
        
          m.insert({qbs_sup, "quadratic_b_spline"});
        }

        {
          auto cbs = boost::math::interpolators::cardinal_cubic_b_spline(phi.data(), phi.size(), Real(0), dx, phi_prime.front(), phi_prime.back());
          Real cbs_sup = 0;
          for (size_t i = 0; i < phi_dense.size(); ++i) {
            Real x = i*dx_dense;
            Real diff = abs(phi_dense[i] - cbs(x));
            if (diff > cbs_sup) {
              cbs_sup = diff;
            }
          }
          
          m.insert({cbs_sup, "cubic_b_spline"});
        }


        // Whittaker-Shannon interpolation has linear complexity; test over all points and it's quadratic.
        // I ran this a couple times and found it's not competitive; so comment out for now.
        /*{
            auto phi_copy = phi;
            auto ws = boost::math::interpolators::whittaker_shannon(std::move(phi_copy), Real(0), dx);
            Real sup = 0;
            for (size_t i = 0; i < phi_dense.size(); ++i) {
              Real x = i*dx_dense;
              using std::abs;
              Real diff = abs(phi_dense[i] - ws(x));
              if (diff > sup) {
                sup = diff;
              }
            }
          
            m.insert({sup, "whittaker_shannon"});
        }*/

        {
          auto qbs = boost::math::interpolators::cardinal_quintic_b_spline(phi.data(), phi.size(), Real(0), dx, {0,0}, {0,0});
          Real sup = 0;
          for (size_t i = 0; i < phi_dense.size(); ++i) {
            Real x = i*dx_dense;
            Real diff = abs(phi_dense[i] - qbs(x));
            if (diff > sup) {
              sup = diff;
            }
          }
          m.insert({sup, "quintic_b_spline"});
        }

        {
          auto phi_copy = phi;
          auto x_copy = x;
          auto phi_prime_copy = phi_prime;
          auto ch = boost::math::interpolators::cubic_hermite(std::move(x_copy), std::move(phi_copy), std::move(phi_prime_copy));
          Real sup = 0;
          for (size_t i = 0; i < phi_dense.size(); ++i) {
            Real x = i*dx_dense;
            Real diff = abs(phi_dense[i] - ch(x));
            if (diff > sup) {
              sup = diff;
            }
          }
          m.insert({sup, "cubic_hermite_spline"});
        }

        {
          auto phi_copy = phi;
          auto x_copy = x;
          auto phi_prime_copy = phi_prime;
          auto pc = boost::math::interpolators::pchip(std::move(x_copy), std::move(phi_copy));
          Real sup = 0;
          for (size_t i = 0; i < phi_dense.size(); ++i) {
            Real x = i*dx_dense;
            Real diff = abs(phi_dense[i] - pc(x));
            if (diff > sup) {
              sup = diff;
            }
          }
          m.insert({sup, "pchip"});
        }

        {
          auto phi_copy = phi;
          auto x_copy = x;
          auto pc = boost::math::interpolators::makima(std::move(x_copy), std::move(phi_copy));
          Real sup = 0;
          for (size_t i = 0; i < phi_dense.size(); ++i) {
            Real x = i*dx_dense;
            Real diff = abs(phi_dense[i] - pc(x));
            if (diff > sup) {
              sup = diff;
            }
          }
          m.insert({sup, "makima"});
        }

        // Again, linear complexity of evaluation => quadratic complexity of exhaustive checking.
        /*{
            auto trig = boost::math::interpolators::cardinal_trigonometric(phi, Real(0), dx);
            Real sup = 0;
            for (size_t i = 0; i < phi_dense.size(); ++i) {
              Real x = i*dx_dense;
              using std::abs;
              Real diff = abs(phi_dense[i] - trig(x));
              if (diff > sup) {
                sup = diff;
              }
            }
            m.insert({sup, "trig"});
        }*/

        {
            auto fotaylor = [&phi, &phi_prime, &r](Real x)->Real {
                if (x <= 0 || x >= 2*p-1) {
                    return 0;
                }
                using std::floor;

                Real y = (1<<r)*x;
                Real k = floor(y);

                size_t kk = static_cast<size_t>(k);
                if (y - k < k + 1 - y)
                {
                    Real eps = (y-k)/(1<<r);
                    return phi[kk] + eps*phi_prime[kk];
                }
                else {
                    Real eps = (y-k-1)/(1<<r);
                    return phi[kk+1] + eps*phi_prime[kk+1];
                }

            };
            Real sup = 0;
            for (size_t i = 0; i < phi_dense.size(); ++i) {
              Real x = i*dx_dense;
              Real diff = abs(phi_dense[i] - fotaylor(x));
              if (diff > sup) {
                sup = diff;
              }
            }
        
          m.insert({sup, "First-order Taylor"});
        }

        if constexpr (p > 2) {
            auto phi_dbl_prime = boost::math::detail::dyadic_grid<Real, p, 2>(r);

            {
                auto phi_copy = phi;
                auto x_copy = x;
                auto phi_prime_copy = phi_prime;
                auto phi_dbl_prime_copy = phi_dbl_prime;
                auto qh = boost::math::interpolators::quintic_hermite(std::move(x_copy), std::move(phi_copy), std::move(phi_prime_copy), std::move(phi_dbl_prime_copy));
                Real sup = 0;
                for (size_t i = 0; i < phi_dense.size(); ++i) {
                  Real x = i*dx_dense;
                  Real diff = abs(phi_dense[i] - qh(x));
                  if (diff > sup) {
                    sup = diff;
                  }
                }
                m.insert({sup, "quintic_hermite_spline"});
            }

            {
                auto sotaylor = [&phi, &phi_prime, &phi_dbl_prime, &r](Real x)->Real {
                      if (x <= 0 || x >= 2*p-1) {
                          return 0;
                      }
                      using std::floor;

                      Real y = (1<<r)*x;
                      Real k = floor(y);

                      size_t kk = static_cast<size_t>(k);
                      if (y - k < k + 1 - y)
                      {
                          Real eps = (y-k)/(1<<r);
                          return phi[kk] + eps*phi_prime[kk] + eps*eps*phi_dbl_prime[kk]/2;
                      }
                      else {
                          Real eps = (y-k-1)/(1<<r);
                          return phi[kk+1] + eps*phi_prime[kk+1] + eps*eps*phi_dbl_prime[kk+1]/2;
                      }
                };
                Real sup = 0;
                for (size_t i = 0; i < phi_dense.size(); ++i) {
                  Real x = i*dx_dense;
                  Real diff = abs(phi_dense[i] - sotaylor(x));
                  if (diff > sup) {
                    sup = diff;
                  }
                }
            
              m.insert({sup, "Second-order Taylor"});
            }
        }

        if constexpr (p > 3) {
            auto phi_dbl_prime = boost::math::detail::dyadic_grid<Real, p, 2>(r);
            auto phi_triple_prime = boost::math::detail::dyadic_grid<Real, p, 3>(r);

            {
                auto totaylor = [&phi, &phi_prime, &phi_dbl_prime, &phi_triple_prime, &r](Real x)->Real {
                      if (x <= 0 || x >= 2*p-1) {
                          return 0;
                      }
                      using std::floor;

                      Real y = (1<<r)*x;
                      Real k = floor(y);

                      size_t kk = static_cast<size_t>(k);
                      if (y - k < k + 1 - y)
                      {
                          Real eps = (y-k)/(1<<r);
                          return phi[kk] + eps*phi_prime[kk] + eps*eps*phi_dbl_prime[kk]/2 + eps*eps*eps*phi_triple_prime[kk]/6;
                      }
                      else {
                          Real eps = (y-k-1)/(1<<r);
                          return phi[kk+1] + eps*phi_prime[kk+1] + eps*eps*phi_dbl_prime[kk+1]/2 + eps*eps*eps*phi_triple_prime[kk]/6;
                      }
                };
                Real sup = 0;
                for (size_t i = 0; i < phi_dense.size(); ++i) {
                  Real x = i*dx_dense;
                  Real diff = abs(phi_dense[i] - totaylor(x));
                  if (diff > sup) {
                    sup = diff;
                  }
                }
            
              m.insert({sup, "Third-order Taylor"});
            }


        }
        std::string best = "none";
        Real best_sup = 1000000000;
        std::cout << std::setprecision(20) << std::fixed;
        for (auto & e : m) {
          
          std::cout << "\t" << e.first << " is error of " << e.second << "\n";
          if (e.first < best_sup) {
            best = e.second;
            best_sup = e.first;
          }
        }
        std::cout << "\tThe best method for p = " << p << " is the " << best << "\n";
    }
}

int main() {

    //do_ulp<float, 5>();

    //choose_refinement<float, 5>();
    //choose_refinement<double, 15>();
    // Says linear interpolation is the best:
    find_best_interpolator<long double, 2>();

    // Says linear interpolation is the best:
    find_best_interpolator<long double, 3>();

    // Says cubic_hermite_spline is best:
    find_best_interpolator<long double, 4>();

    // Says cubic_hermite_spline is best:
    find_best_interpolator<long double, 5>();

    // Says quintic_hermite_spline is best:
    find_best_interpolator<long double, 6>();

    // Says quintic_hermite_spline is best:
    find_best_interpolator<long double, 7>();

    // Says quintic_hermite_spline is best:
    find_best_interpolator<long double, 15>();
}
