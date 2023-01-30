#include <boost/math/tools/rational.hpp>
#include <boost/math/tools/estrin.hpp>
#include <boost/math/tools/ulps_plot.hpp>

using boost::math::tools::evaluate_polynomial;
using boost::math::tools::estrin;
using boost::math::tools::ulps_plot;

void estrin_vs_horner(size_t n) {
   std::random_device rd;
   auto seed = rd();
   std::mt19937_64 gen(seed);
   std::uniform_real_distribution<float> dis(-1, 1);
   std::vector<float> coeffs_float(n);
   std::vector<double> coeffs_double(n);
   for (size_t i = 0; i < n; ++i) {
        coeffs_float[i] = dis(gen);
	coeffs_double[i] = coeffs_float[i];
   }

   auto hi_acc = [&](double x) {
       return evaluate_polynomial(coeffs_double.data(), x, coeffs_double.size());
   };
   auto estrin_float = [&](float x) {
       return estrin(coeffs_float, x);
   };
   auto horner_float = [&](float x) {
       return evaluate_polynomial(coeffs_float.data(), x, coeffs_float.size());
   };
   auto plot = ulps_plot<decltype(hi_acc), double, float>(hi_acc, float(-2.0), float(2.0), 20000);
   plot.ulp_envelope(false);
   plot.add_fn(estrin_float, "steelblue");
   plot.add_fn(horner_float, "orange");
   plot.clip(10);
   plot.title("Horner (orange) vs Estrin (blue) accuracy");
   plot.write("horner_vs_estrin.svg");
}


int main() {
   size_t n = 30;
   estrin_vs_horner(n);
   return 0;
}
