//  (C) Copyright John Maddock 2006.
//  (C) Copyright Paul A. Bristow 2006.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_STATS_STUDENTS_T_HPP
#define BOOST_STATS_STUDENTS_T_HPP

// http://en.wikipedia.org/wiki/Student%27s_t_distribution
// http://www.itl.nist.gov/div898/handbook/eda/section3/eda3664.htm

#include <boost/math/special_functions/beta.hpp> // for ibeta(a, b, x).
#include <boost/math/distributions/complement.hpp>

namespace boost{ namespace math{

template <class RealType>
class students_t_distribution
{
public:
   students_t_distribution(RealType i) : m_df(i)
   {
       if(m_df <= 0)
       {  // Degrees of freedom must be > 0!
          // may or may not throw.
          tools::domain_error<RealType>(BOOST_CURRENT_FUNCTION, "degrees of freedom argument is %1%, but must be > 0 !", m_df);
       }
   }

   RealType degrees_of_freedom()const
   {
      return m_df;
   }

   // Parameter estimation:
   static inline RealType estimate_degrees_of_freedom(
      RealType M,            // true mean
      RealType Sm,           // sample mean
      RealType Sd,           // sample std deviation
      RealType p,            // required significance level
      RealType hint = 100)   // hint to where the answer lies.
   {
      return estimate_degrees_of_freedom_imp(M, Sm, Sd, p, 1-p, hint);
   }

   template <class R1, class R2, class R3, class R4>
   static inline RealType estimate_degrees_of_freedom(
      const complemented4_type<R1, R2, R3, R4>& c)
   {
      return estimate_degrees_of_freedom_imp(c.dist, c.param1, c.param2, 1 - c.param3, c.param3);
   }
   template <class R1, class R2, class R3, class R4, class R5>
   static inline RealType estimate_degrees_of_freedom(
      const complemented5_type<R1, R2, R3, R4, R5>& c)
   {
      return estimate_degrees_of_freedom_imp(c.dist, c.param1, c.param2, 1 - c.param3, c.param3, c.param4);
   }

   static RealType estimate_two_equal_degrees_of_freedom(
      RealType Sm1,           // sample 1 mean
      RealType Sd1,           // sample 1 std deviation
      RealType Sm2,           // sample 2 mean
      RealType Sd2,           // sample 2 std deviation
      RealType p,             // required probablity
      RealType hint = 100)    // hint to where the answer lies.
   {
      return estimate_two_equal_degrees_of_freedom_imp(Sm1, Sd1, Sm2, Sd2, p, 1 - p, hint);
   }
   template <class R1, class R2, class R3, class R4, class R5>
   static inline RealType estimate_two_equal_degrees_of_freedom(
      const complemented5_type<R1, R2, R3, R4, R5>& c)
   {
      return estimate_two_equal_degrees_of_freedom_imp(c.dist, c.param1, c.param2, c.param3, 1 - c.param4, c.param4);
   }
   template <class R1, class R2, class R3, class R4, class R5, class R6>
   static inline RealType estimate_two_equal_degrees_of_freedom(
      const complemented6_type<R1, R2, R3, R4, R5, R6>& c)
   {
      return estimate_two_equal_degrees_of_freedom_imp(c.dist, c.param1, c.param2, c.param3, 1 - c.param4, c.param4, c.param5);
   }

   static RealType estimate_two_unequal_degrees_of_freedom(
      RealType Sm1,           // sample 1 mean
      RealType Sd1,           // sample 1 std deviation
      unsigned Sn1,           // sample 1 fixed size
      RealType Sm2,           // sample 2 mean
      RealType Sd2,           // sample 2 std deviation
      RealType p,             // required probablity
      RealType hint = 100)    // hint to where the answer lies.
   {
      return estimate_two_unequal_degrees_of_freedom_imp(Sm1, Sd1, Sn1, Sm2, Sd2, p, 1 - p, hint);
   }
   template <class R1, class R2, class R3, class R4, class R5, class R6>
   static inline RealType estimate_two_unequal_degrees_of_freedom(
      const complemented6_type<R1, R2, R3, R4, R5, R6>& c)
   {
      return estimate_two_unequal_degrees_of_freedom_imp(c.dist, c.param1, c.param2, c.param3, c.param4, 1 - c.param5, c.param5);
   }
   template <class R1, class R2, class R3, class R4, class R5, class R6, class R7>
   static inline RealType estimate_two_unequal_degrees_of_freedom(
      const complemented7_type<R1, R2, R3, R4, R5, R6, R7>& c)
   {
      return estimate_two_unequal_degrees_of_freedom_imp(c.dist, c.param1, c.param2, c.param3, c.param4, 1 - c.param5, c.param5, c.param6);
   }

private:
   //
   // Private implementation methods for parameter estimation:
   //
   static RealType estimate_degrees_of_freedom_imp(
      RealType M,            // true mean
      RealType Sm,           // sample mean
      RealType Sd,           // sample std deviation
      RealType p,            // required probablity
      RealType q,            // 1 - required probablity
      RealType hint = 100);  // hint to where the answer lies.

   static RealType estimate_two_equal_degrees_of_freedom_imp(
      RealType Sm1,           // sample 1 mean
      RealType Sd1,           // sample 1 std deviation
      RealType Sm2,           // sample 2 mean
      RealType Sd2,           // sample 2 std deviation
      RealType p,             // required probablity
      RealType q,             // 1 - required probablity
      RealType hint = 100);   // hint to where the answer lies.

   static RealType estimate_two_unequal_degrees_of_freedom_imp(
      RealType Sm1,           // sample 1 mean
      RealType Sd1,           // sample 1 std deviation
      unsigned Sn1,           // sample 1 size (fixed)
      RealType Sm2,           // sample 2 mean
      RealType Sd2,           // sample 2 std deviation
      RealType p,             // required probablity
      RealType q,             // 1 - required probablity
      RealType hint = 100);   // hint to where the answer lies.
   //
   // Data members:
   //
   RealType m_df;  // degrees of freedom are a real number.
};

typedef students_t_distribution<double> students_t;

template <class RealType>
RealType pdf(const students_t_distribution<RealType>& dist, const RealType& t)
{
   using namespace std;  // for ADL of std functions
   RealType degrees_of_freedom = dist.degrees_of_freedom();
   if(degrees_of_freedom <= 0)
   { // Degrees of freedom must be > 0!
      return tools::domain_error<RealType>(BOOST_CURRENT_FUNCTION, "degrees of freedom argument is %1%, but must be > 0 !", degrees_of_freedom);
   }
   // TODO fails for t == 0 and df >=1e16 for ALL fp types.
   RealType basem1 = t * t / degrees_of_freedom;
   RealType result;
   if(basem1 < 0.125)
   {
      result = exp(-boost::math::log1p(basem1) * (1+degrees_of_freedom) / 2);
   }
   else
   {
      result = pow(1 / (1 + basem1), (degrees_of_freedom + 1) / 2);
   }
   result /= sqrt(degrees_of_freedom) * boost::math::beta(degrees_of_freedom / 2, RealType(0.5f));

   return result;
}

template <class RealType>
RealType cdf(const students_t_distribution<RealType>& dist, const RealType& t)
{
   RealType degrees_of_freedom = dist.degrees_of_freedom();
   if(degrees_of_freedom <= 0)
   { // Degrees of freedom must be > 0!
      return tools::domain_error<RealType>(BOOST_CURRENT_FUNCTION, "degrees of freedom argument is %1%, but must be > 0 !", degrees_of_freedom);
   }
   if (t == 0)
   {
     return 0.5;
   }
   RealType z = degrees_of_freedom / (degrees_of_freedom + t * t);
   // Calculate probability of Student's t using the incomplete beta function.
   // probability = ibeta(degrees_of_freedom/2, 1/2, degrees_of_freedom/ (degrees_of_freedom + t*t))
   RealType probability = ibeta(degrees_of_freedom / 2, static_cast<RealType>(0.5), z) / 2;
   // Check 0 <= probability probability <= 1.
   // Numerical errors might cause probability to be slightly outside the range < 0 or > 1.
   // This might cause trouble downstream, so warn, possibly throw exception, but constrain to the limits.
   if (probability < static_cast<RealType>(0.))
   {
      tools::logic_error<RealType>(BOOST_CURRENT_FUNCTION, "probability %1% is < 0, so has been constrained to zero !", probability);
      return static_cast<RealType>(0.); // Constrain to zero if logic_error does not throw.
   }
   if(probability > static_cast<RealType>(1.))
   {
      tools::logic_error<RealType>(BOOST_CURRENT_FUNCTION, "probability %1% is > 1, so has been constrained to unity!", probability);
      return static_cast<RealType>(1.); // Constrain to unity if logic_error does not throw.
   }
   return (t > 0 ? 1	- probability : probability);
}

template <class RealType>
RealType quantile(const students_t_distribution<RealType>& dist, const RealType& p)
{
   using namespace std; // for ADL of std functions
   // Error check.
   RealType degrees_of_freedom = dist.degrees_of_freedom();
   RealType probability = p;
   if(degrees_of_freedom <= 0)
   { // Degrees of freedom must be > 0!
      return tools::domain_error<RealType>(BOOST_CURRENT_FUNCTION, "degrees of freedom argument is %1%, but must be > 0 !", degrees_of_freedom);
   }
   if((probability < 0) || (probability > 1))
   { // probability must be >= 0 and <= 1!
      return tools::domain_error<RealType>(BOOST_CURRENT_FUNCTION, "probability argument is %1%, but must be >= 0 and <= 1 !", probability);
   }
   // Special cases, regardless of degrees_of_freedom.
   if (probability == 0)
      return numeric_limits<RealType>::has_infinity ? -numeric_limits<RealType>::infinity() : -tools::max_value<RealType>();
   if (probability == 1)
     return numeric_limits<RealType>::has_infinity ? numeric_limits<RealType>::infinity() : tools::max_value<RealType>();
   if (probability == static_cast<RealType>(0.5))
     return 0;
   // Calculate quantile of Student's t using the incomplete beta function.
   if ((probability > RealType(0.25)) && (probability < RealType(0.75)) )
   { // probability is middling.
     RealType z = 1 - 2 * probability;
     z = ibeta_inv(RealType(0.5), degrees_of_freedom / 2 , fabs(z));
     RealType t = sqrt(degrees_of_freedom * z / (1 - z));
     return (probability < RealType(0.5)) ? -t : t;
   }
   else
   { // Probability is small or large, we may get overflow unless we're careful.
     int sign = -1;
     if (probability >= 0.5)
     { // large and > 0.75
	     probability = 1 - probability;
	     sign = +1;
     }
     RealType z = ibeta_inv(degrees_of_freedom / 2, static_cast<RealType>(0.5), 2 * probability);
     if ((tools::max_value<RealType>() * z) < degrees_of_freedom)
     {
        return sign * tools::max_value<RealType>();
     }
     else
     {
	     return sign * sqrt(degrees_of_freedom/ z - degrees_of_freedom);
     }
   }
}

template <class RealType>
RealType cdf(const complemented2_type<students_t_distribution<RealType>, RealType>& c)
{
   return cdf(c.dist, -c.param);
}

template <class RealType>
RealType quantile(const complemented2_type<students_t_distribution<RealType>, RealType>& c)
{
   return -quantile(c.dist, c.param);
}

//
// Parameter estimation follows:
//
namespace detail{
//
// Functors for finding degrees of freedom:
//
template <class RealType>
struct single_case_df
{
   single_case_df(RealType md, RealType sd, RealType c, bool inv)
      : MeanDiff(md), Sd(sd), p(c), invert(inv) {}

   RealType operator()(const RealType df)
   {
      using namespace std;  // for ADL
      RealType N = df + 1;
      RealType t = MeanDiff * sqrt(double(N)) / Sd;
      students_t_distribution<RealType> dist(df);
      return invert ? p - cdf(complement(dist, t)) : cdf(dist, t) - p;
   }
private:
   RealType MeanDiff;      // Difference in means.
   RealType Sd;            // Sample std deviation.
   RealType p;             // Probability.
   bool invert;            // Whether p has been complemented.
};

template <class RealType>
struct two_case_equal_N_df
{
   two_case_equal_N_df(RealType m1, RealType sd1, RealType m2, RealType sd2, RealType c, bool inv)
      : M1(m1), M2(m2), Sd1(sd1), Sd2(sd2), p(c), invert(inv) {}

   RealType operator()(const RealType df)
   {
      using namespace std;  // for ADL
      RealType N = (df + 2) / 2;
      RealType sp = sqrt(((N-1) * Sd1 * Sd1 + (N-1) * Sd2 * Sd2) / df);
      RealType t = fabs(M1 - M2) / (sp * sqrt(2/N));
      students_t_distribution<RealType> dist(df);
      return invert ? p - cdf(complement(dist, t)) : cdf(dist, t) - p;
   }
private:
   RealType M1;            // Sample 1 Mean.
   RealType Sd1;           // Sample 1 std deviation.
   RealType M2;            // Sample 2 Mean.
   RealType Sd2;           // Sample 2 std deviation.
   RealType p;             // Probability.
   bool invert;            // Whether p has been complemented.
};

template <class RealType>
struct two_case_unequal_N_df
{
   two_case_unequal_N_df(RealType m1, RealType sd1, unsigned sn1, RealType m2, RealType sd2, RealType c, bool inv)
      : M1(m1), Sd1(sd1), Sn1(sn1), M2(m2), Sd2(sd2), p(c), invert(inv) {}

   RealType operator()(const RealType df)
   {
      using namespace std;  // for ADL
      RealType Sn2 = (df + 2) - Sn1;
      if(Sn2 <= 0)
      {
         // we can't have a negative sample size
         // the root must be above this point so
         // return a false value:
         return -1;
      }
      RealType sp = sqrt(((Sn1-1) * Sd1 * Sd1 + (Sn2-1) * Sd2 * Sd2) / df);
      RealType t = fabs(M1 - M2) / (sp * sqrt(1/Sn1 + 1/Sn2));
      students_t_distribution<RealType> dist(df);
      return invert ? p - cdf(complement(dist, t)) : cdf(dist, t) - p;
   }
private:
   RealType M1;            // Sample 1 Mean.
   RealType Sd1;           // Sample 1 std deviation.
   unsigned Sn1;           // Sample 1 fixed size.
   RealType M2;            // Sample 2 Mean.
   RealType Sd2;           // Sample 2 std deviation.
   RealType p;             // Probability.
   bool invert;            // Whether p has been complemented.
};

}  // namespace detail

template <class RealType>
RealType students_t_distribution<RealType>::estimate_degrees_of_freedom_imp(
   RealType M,                // true mean
   RealType Sm,               // sample mean
   RealType Sd,               // sample std deviation
   RealType p,                // required significance level p
   RealType q,                // required significance level q
   RealType hint /*= 100*/)   // hint to where the answer lies.
{
   using namespace std;
   detail::single_case_df<RealType> f(
      fabs(M-Sm),
      Sd,
      (p < q ? p : q),
      (p < q ? false : true));
   tools::eps_tolerance<RealType> tol(tools::digits<RealType>());
   boost::uintmax_t max_iter = 10000;
   std::pair<RealType, RealType> r = tools::bracket_and_solve_root(f, hint, RealType(2), true, tol, max_iter);
   RealType result = r.first + (r.second - r.first) / 2;
   if(max_iter == 10000)
   {
      tools::logic_error<RealType>(BOOST_CURRENT_FUNCTION, "Unable to locate solution in a reasonable time:"
         " either there is no answer to how many degrees of freedom are required"
         " or the answer is infinite.  Current best guess is %1%", result);
   }
   return result;
}

template <class RealType>
RealType students_t_distribution<RealType>::estimate_two_equal_degrees_of_freedom_imp(
   RealType Sm1,           // sample 1 mean
   RealType Sd1,           // sample 1 std deviation
   RealType Sm2,           // sample 2 mean
   RealType Sd2,           // sample 2 std deviation
   RealType p,             // required probablity
   RealType q,             // 1 - required probablity
   RealType hint /*= 100*/)    // hint to where the answer lies.
{
   using namespace std;
   RealType result;
   try{
      detail::two_case_equal_N_df<RealType> f(
         Sm1,
         Sd1,
         Sm2,
         Sd2,
         (p < q ? p : q),
         (p < q ? false : true));
      tools::eps_tolerance<RealType> tol(tools::digits<RealType>());
      boost::uintmax_t max_iter = 10000;
      std::pair<RealType, RealType> r = tools::bracket_and_solve_root(f, hint, RealType(2), true, tol, max_iter);
      result = r.first + (r.second - r.first) / 2;
      if(max_iter == 10000)
      {
         if(result > hint)
            result = tools::max_value<RealType>();
         else
            result = 0;
      }
   }
   catch(const std::logic_error&)
   {
      detail::two_case_equal_N_df<RealType> f(
         Sm1,
         Sd1,
         Sm2,
         Sd2,
         (p < q ? p : q),
         (p < q ? false : true));
      if(f(hint) > 0)
         result = 0;
      else
         result = tools::max_value<RealType>();
   }
   return result;
}

template <class RealType>
RealType students_t_distribution<RealType>::estimate_two_unequal_degrees_of_freedom_imp(
   RealType Sm1,           // sample 1 mean
   RealType Sd1,           // sample 1 std deviation
   unsigned Sn1,           // sample 1 fixed size
   RealType Sm2,           // sample 2 mean
   RealType Sd2,           // sample 2 std deviation
   RealType p,             // required probablity
   RealType q,             // 1 - required probablity
   RealType hint /*= 100*/)    // hint to where the answer lies.
{
   using namespace std;
   RealType result;
   try{
      if(hint < Sn1 - 1)
         hint = Sn1 - 1;
      detail::two_case_unequal_N_df<RealType> f(
         Sm1,
         Sd1,
         Sn1,
         Sm2,
         Sd2,
         (p < q ? p : q),
         (p < q ? false : true));
      tools::eps_tolerance<RealType> tol(tools::digits<RealType>());
      boost::uintmax_t max_iter = 1000;
      std::pair<RealType, RealType> r = tools::bracket_and_solve_root(f, hint, RealType(2), true, tol, max_iter);
      result = r.first + (r.second - r.first) / 2;
      if(max_iter == 10000)
      {
         if(result > hint)
            result = tools::max_value<RealType>();
         else
            result = Sn1 - 1;
      }
   }
   catch(const std::logic_error&)
   {
      detail::two_case_unequal_N_df<RealType> f(
         Sm1,
         Sd1,
         Sn1,
         Sm2,
         Sd2,
         (p < q ? p : q),
         (p < q ? false : true));
      if(f(hint) > 0)
         result = Sn1 - 1;
      else
         result = tools::max_value<RealType>();
   }
   return result;
}

} // namespace math
} // namespace boost

#endif // BOOST_STATS_STUDENTS_T_HPP

