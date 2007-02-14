// boost_math.h

#pragma once

using namespace System;

namespace boost_math {

   class any_imp
   {
   public:
     // Distribution properties.
      virtual double mean()const = 0;
      virtual double mode()const = 0;
      virtual double median()const = 0;
      virtual double variance()const = 0;
      virtual double standard_deviation()const = 0;
      virtual double skewness()const = 0;
      virtual double kurtosis()const = 0;
      virtual double kurtosis_excess()const = 0;
      virtual double coefficient_of_variation()const = 0;
      // Values computed from random variate x.
      virtual double hazard(double x)const = 0;
      virtual double chf(double x)const = 0;
      virtual double cdf(double x)const = 0;
      virtual double ccdf(double x)const = 0;
      virtual double pdf(double x)const = 0;
      virtual double quantile(double x)const = 0;
      virtual double quantile_c(double x)const = 0;
      // Range & support of x
      virtual double lowest()const = 0;
      virtual double uppermost()const = 0;
      virtual double lower()const = 0;
      virtual double upper()const = 0;
   };

   template <class Distribution>
   class concrete_distribution : public any_imp
   {
   public:
      concrete_distribution(const Distribution& d) : m_dist(d) {}
      // Distribution properties.
      virtual double mean()const
      {
         return boost::math::mean(m_dist);
      }
      virtual double median()const
      {
         return boost::math::median(m_dist);
      }
      virtual double mode()const
      {
         return boost::math::mode(m_dist);
      }
      virtual double variance()const
      {
         return boost::math::variance(m_dist);
      }
      virtual double skewness()const
      {
         return boost::math::skewness(m_dist);
      }
      virtual double standard_deviation()const
      {
         return boost::math::standard_deviation(m_dist);
      }
      virtual double coefficient_of_variation()const
      {
         return boost::math::coefficient_of_variation(m_dist);
      }
      virtual double kurtosis()const
      {
         return boost::math::kurtosis(m_dist);
      }
      virtual double kurtosis_excess()const
      {
         return boost::math::kurtosis_excess(m_dist);
      }
      // Range of x for the distribution.
      virtual double lowest()const
      {
         return boost::math::range(m_dist).first;
      }
      virtual double uppermost()const
      {
         return boost::math::range(m_dist).second;
      }
      // Support of x for the distribution.
      virtual double lower()const
      {
         return boost::math::support(m_dist).first;
      }
      virtual double upper()const
      {
         return boost::math::support(m_dist).second;
      }

      // Values computed from random variate x.
      virtual double hazard(double x)const
      {
         return boost::math::hazard(m_dist, x);
      }
      virtual double chf(double x)const
      {
         return boost::math::chf(m_dist, x);
      }
      virtual double cdf(double x)const
      {
         return boost::math::cdf(m_dist, x);
      }
      virtual double ccdf(double x)const
      {
         return boost::math::cdf(complement(m_dist, x));
      }
      virtual double pdf(double x)const
      {
         return boost::math::pdf(m_dist, x);
      }
      virtual double quantile(double x)const
      {
         return boost::math::quantile(m_dist, x);
      }
      virtual double quantile_c(double x)const
      {
         return boost::math::quantile(complement(m_dist, x));
      }
   private:
      Distribution m_dist;
   };

	public ref class any_distribution
	{
	  public:
		// Added methods for this class here.
      any_distribution(int t, double arg1, double arg2, double arg3);
      ~any_distribution()
      {
         reset(0);
      }
      // Is it OK for these to be inline?
      // Distribution properties as 'pointer-to-implementions'.
      double mean()
      {
         return pimpl->mean();
      }
      double median()
      {
         return pimpl->median();
      }
      double mode()
      {
         return pimpl->mode();
      }
      double variance()
      {
         return pimpl->variance();
      }
      double standard_deviation()
      {
         return pimpl->standard_deviation();
      }
      double coefficient_of_variation()
      { // aka Relative Standard deviation.
         return pimpl->coefficient_of_variation();
      }
      double skewness()
      {
         return pimpl->skewness();
      }
      double kurtosis()
      {
         return pimpl->kurtosis();
      }
      double kurtosis_excess()
      {
         return pimpl->kurtosis_excess();
      }
      // Values computed from random variate x.
      double hazard(double x)
      {
         return pimpl->hazard(x);
      }
      double chf(double x)
      {
         return pimpl->chf(x);
      }
      double cdf(double x)
      {
         return pimpl->cdf(x);
      }
      double ccdf(double x)
      {
         return pimpl->ccdf(x);
      }
      double pdf(double x)
      {
         return pimpl->pdf(x);
      }
      double quantile(double x)
      {
         return pimpl->quantile(x);
      }
      double quantile_c(double x)
      {
         return pimpl->quantile_c(x);
      }

      double lowest()
      {
         return pimpl->lowest();
      }

      double uppermost()
      {
         return pimpl->uppermost();
      }

      double lower()
      {
         return pimpl->lower();
      }
      double upper()
      {
         return pimpl->upper();
      }

      // How many distributions are supported:
      static int size();
      // Display name of i'th distribution:
      static System::String^ distribution_name(int i);
      // Name of first distribution parameter, or null if not supported:
      static System::String^ first_param_name(int i);
      // Name of second distribution parameter, or null if not supported:
      static System::String^ second_param_name(int i);
      // Name of third distribution parameter, or null if not supported:
      static System::String^ third_param_name(int i);
      // Default value for first parameter:
      static double first_param_default(int i);
      // Default value for second parameter:
      static double second_param_default(int i);
      // Default value for third parameter:
      static double third_param_default(int i);

   private:
      any_distribution(const any_distribution^)
      { // Constructor is private.
      }
      const any_distribution^ operator=(const any_distribution^ d)
      { // Copy Constructor is private too.
         return d;
      }
      // We really should use a shared_ptr here, 
      // but apparently it's not allowed in a managed class like this :-(
      void reset(any_imp* p)
      {
         if(pimpl)
         { // Exists already, so
            delete pimpl;
         }
         pimpl = p;
      }
      any_imp* pimpl;
	};
}
