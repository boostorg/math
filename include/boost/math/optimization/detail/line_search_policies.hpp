#ifndef LINE_SEARCH_POLICIES_HPP
#define LINE_SEARCH_POLICIES_HPP

#include <cmath>
#include <iostream>
#include <numeric>
#include <vector>
namespace boost {
namespace math {
namespace optimization {

/**
 * @brief> Armijo condition backtracking line search
 * https://en.wikipedia.org/wiki/Backtracking_line_search
 *
 *  f(x+alpha p) <= f(x) + alpha * c * grad(f)^T p
 * */
template <typename RealType> class armijo_line_search_policy {
private:
  RealType alpha0_; // initial step size
  RealType c_;      // sufficient decrease constant
  RealType rho_;    // backtracking factor
  int max_iter_;    // maximum backtracking steps

public:
  armijo_line_search_policy(RealType alpha0 = 1.0, RealType c = 1e-4,
                            RealType rho = 0.5, int max_iter = 20)
      : alpha0_(alpha0), c_(c), rho_(rho), max_iter_(max_iter) {}

  template <class Objective, class ObjectiveEvalPolicy,
            class GradientEvalPolicy, class ArgumentContainer>
  RealType operator()(Objective &objective, ObjectiveEvalPolicy &obj_eval,
                      GradientEvalPolicy &grad_eval, ArgumentContainer &x,
                      const std::vector<RealType> &g,
                      const std::vector<RealType> &p, RealType f_x) const {
    /** @brief> line search
     * */
    RealType alpha = alpha0_;
    ArgumentContainer x_trial = x; // copy
    const RealType gTp = dot(g, p);

    for (int iter = 0; iter < max_iter_; ++iter) {
      for (size_t i = 0; i < x.size(); ++i)
        x_trial[i] = x[i] + alpha * p[i];
      auto f_trial = obj_eval(objective, x_trial);
      if (f_trial <=
          f_x + c_ * alpha * gTp) // check if armijo condition is satisfied
        return alpha;
      alpha *= rho_; // half by default
    }
    return alpha;
  }
};

} // namespace optimization
} // namespace math
} // namespace boost
#endif // LINE_SEARCH_POLICIES_HPP
