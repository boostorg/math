//  Copyright Jeremy Murphy 2017.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifdef _MSC_VER
#  pragma warning (disable : 4224)
#endif

#include <boost/math/common_factor_rt.hpp>
#include <boost/math/special_functions/prime.hpp>
#include <boost/math/tools/polynomial.hpp>
#include <boost/math/tools/polynomial_gcd.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/multiprecision/cpp_int.hpp>
#include <boost/multiprecision/integer.hpp>
#include <boost/random.hpp>
#include <boost/array.hpp>
#include <boost/type_traits.hpp>

#include <iostream>
#include <algorithm>
#include <numeric>
#include <string>
#include <tuple>
#include <type_traits>
#include <vector>
#include <functional>
#include <typeinfo>

#include "fibonacci.hpp"
#include "../../test/table_type.hpp"
#include "table_helper.hpp"
#include "performance.hpp"


using namespace std;
using namespace boost::math::tools;
using boost::math::gcd_traits;

polynomial<double> total_sum(0);
polynomial<long int> total_sum_int(0);
typedef boost::multiprecision::cpp_bin_float_quad mp_float_type;
polynomial<mp_float_type> total_sum_mp(0);


template <typename T>
typename boost::enable_if_c<boost::is_floating_point<T>::value, void>::type
add_total(polynomial<T> const &x)
{
    total_sum += x;
}


template <typename T>
typename boost::enable_if_c<boost::is_integral<T>::value, void>::type
add_total(polynomial<T> const &x)
{
    total_sum += x;
}

void add_total(polynomial<mp_float_type> const &x)
{
    total_sum_mp += x;
}

template <typename SteinDomain>
SteinDomain Stein_gcd_Joux(SteinDomain m, SteinDomain n)
{
    using std::swap;
    BOOST_ASSERT(m || n);
    if (!m)
        return n;
    if (!n)
        return m;
    unsigned shifts = std::min(gcd_traits<SteinDomain>::make_odd(m), gcd_traits<SteinDomain>::make_odd(n));
    // odd(m) && odd(n)
    while (m != n)
    {
        if (gcd_traits<SteinDomain>::less(m, n))
            swap(n, m);
        gcd_traits<SteinDomain>::Joux_subtraction(m, n);
        BOOST_ASSERT(even(m));
        // With polynomials for example, it is possible that m is now zero.
        if (!m)
        {
            n <<= shifts;
            return n;
        }
        gcd_traits<SteinDomain>::make_odd(m);
    }
    // m == n
    m <<= shifts;
    return m;
}


template <typename SteinDomain>
SteinDomain Stein_gcd_factored_Joux(SteinDomain m, SteinDomain n)
{
    using std::swap;
    BOOST_ASSERT(m || n);
    if (!m)
        return n;
    if (!n)
        return m;
    unsigned shifts = std::min(gcd_traits<SteinDomain>::make_odd(m), gcd_traits<SteinDomain>::make_odd(n));
    // odd(m) && odd(n)
    while (m != n)
    {
        if (gcd_traits<SteinDomain>::less(m, n))
            swap(n, m);
        gcd_traits<SteinDomain>::factored_Joux_subtraction(m, n);
        BOOST_ASSERT(even(m));
        // With polynomials for example, it is possible that m is now zero.
        if (!m)
        {
            n <<= shifts;
            return n;
        }
        gcd_traits<SteinDomain>::make_odd(m);
    }
    // m == n
    m <<= shifts;
    return m;
}


template <typename Func, class Table>
double exec_timed_test_foo(Func f, const Table& data, double min_elapsed = 0.5)
{
    double t = 0;
    unsigned repeats = 1;
    typename Table::value_type::first_type sum{0};
    stopwatch<boost::chrono::high_resolution_clock> w;
    do
    {
       for(unsigned count = 0; count < repeats; ++count)
       {
          for(typename Table::size_type n = 0; n < data.size(); ++n)
            sum += f(data[n].first, data[n].second);
       }

        t = boost::chrono::duration_cast<boost::chrono::duration<double>>(w.elapsed()).count();
        if(t < min_elapsed)
            repeats *= 2;
    }
    while(t < min_elapsed);
    add_total(sum);
    return t / repeats;
}


template <typename T>
struct test_function_template
{
    vector<pair<T, T> > const & data;
    const char* data_name;
    
    test_function_template(vector<pair<T, T> > const &data, const char* name) : data(data), data_name(name)
    {
        cerr << "Testing: " << name << endl;
    }
    
    template <typename Function>
    void operator()(pair<Function, string> const &f) const
    {
        cerr << "Algorithm: " << f.second << endl;
        auto result = exec_timed_test_foo(f.first, data);
        auto table_name = string("polynomial gcd method comparison with ") + compiler_name() + string(" on ") + platform_name();

        report_execution_time(result, 
                            table_name,
                            string(data_name), 
                            string(f.second) + "\n" + boost_name());
    }
};

boost::random::mt19937 rng;
boost::random::uniform_int_distribution<> d_0_6(0, 6);
boost::random::uniform_int_distribution<> d_1_20(1, 20);

template <class T>
T get_prime_products()
{
   int n_primes = d_0_6(rng);
   switch(n_primes)
   {
   case 0:
      // Generate a power of 2:
      // return static_cast<T>(1u) << d_1_20(rng);
      return std::pow(2, d_1_20(rng));
   case 1:
      // prime number:
      return boost::math::prime(d_1_20(rng) + 3);
   }
   T result = 1;
   for(int i = 0; i < n_primes; ++i)
      result *= boost::math::prime(d_1_20(rng) + 3) * boost::math::prime(d_1_20(rng) + 3) * boost::math::prime(d_1_20(rng) + 3) * boost::math::prime(d_1_20(rng) + 3) * boost::math::prime(d_1_20(rng) + 3);
   return result;
}


// T is integral.
template <class T>
typename boost::enable_if_c<boost::is_integral<T>::value, T>::type
get_uniform_random()
{
    static boost::random::uniform_int_distribution<T> minmax(0, 255);
    return minmax(rng);
}


// T is floating point, limit to small values so that IVS does not crash test.
template <class T>
typename boost::enable_if_c<boost::is_floating_point<T>::value, T>::type
get_uniform_random()
{
    using std::round;
    static boost::random::uniform_real_distribution<T> minmax(0, 3);
    return round(minmax(rng));
}


// T is not a POD.
template <class T>
typename boost::enable_if_c<!boost::is_pod<T>::value, T>::type
get_uniform_random()
{
    using std::round;
    static boost::random::uniform_real_distribution<T> minmax(0, 3);
    return round(minmax(rng));
}


template <class T>
inline bool even(T const& val)
{
   return !(val & 1u);
}


template <class T, class OutputIterator>
typename boost::enable_if_c<boost::is_integral<T>::value, void>::type
gcd_algorithms(OutputIterator output)
{
    *output++ = std::make_pair(boost::math::detail::Stein_gcd< polynomial<T> >, "Stein_gcd");
}

template <class T, class OutputIterator>
typename boost::enable_if_c<boost::is_same<T, boost::multiprecision::cpp_int>::value, void>::type
gcd_algorithms(OutputIterator output)
{
    *output++ = std::make_pair(boost::math::detail::Stein_gcd< polynomial<T> >, "Stein_gcd");
    // *output++ = std::make_pair(boost::math::detail::subresultant_gcd< polynomial<T> >, "subresultant gcd");    
}


template <class T, class OutputIterator>
typename boost::enable_if_c<boost::is_floating_point<T>::value || !std::numeric_limits<T>::is_integer, void>::type
gcd_algorithms(OutputIterator output)
{
    // *output++ = std::make_pair(boost::math::detail::Euclid_gcd< polynomial<T> >, "Euclid_gcd");
    *output++ = std::make_pair(boost::math::detail::Stein_gcd< polynomial<T> >, "Stein_gcd (Stepanov-factored Joux)");
    *output++ = std::make_pair(Stein_gcd_factored_Joux< polynomial<T> >, "Stein_gcd (factored Joux)");
    *output++ = std::make_pair(Stein_gcd_Joux< polynomial<T> >, "Stein_gcd (Joux)");
}


template <class T>
void test_type(const std::string name)
{
   using namespace boost::math::detail;
   std::vector<pair< polynomial<T>, polynomial<T> > > data;

   for(unsigned i = 0; i < 10; ++i)
   {
       data.push_back(pair< polynomial<T>, polynomial<T> >());
       for (unsigned j = 0; j != 5; j++)
       {
           data.back().first.data().push_back(get_uniform_random<T>());
           data.back().second.data().push_back(get_uniform_random<T>());
       }
       data.back().first.normalize();
       data.back().second.normalize();
   }
   std::string row_name("gcd<");
   row_name += name;
   row_name += "> (random prime number products)";

   typedef pair< function<polynomial<T>(polynomial<T>, polynomial<T>)>, string> f_test;
   vector<f_test> test_functions;
   gcd_algorithms<T>(back_inserter(test_functions));
   for_each(begin(test_functions), end(test_functions), test_function_template< polynomial<T> >(data, row_name.c_str()));

   data.clear();
#if 0
   for(unsigned i = 0; i < 1000; ++i)
   {
      data.push_back(std::make_pair(get_uniform_random<T>(), get_uniform_random<T>()));
   }
   row_name.erase();
   row_name += "gcd<";
   row_name += name;
   row_name += "> (uniform random numbers)";
   for_each(begin(test_functions), end(test_functions), test_function_template<int_type>(data, row_name.c_str()));

   // Fibonacci number tests:
   row_name.erase();
   row_name += "gcd<";
   row_name += name;
   row_name += "> (adjacent Fibonacci numbers)";
   for_each(begin(test_functions), end(test_functions), test_function_template<int_type>(fibonacci_numbers_permution_1<T>(), row_name.c_str()));

   row_name.erase();
   row_name += "gcd<";
   row_name += name;
   row_name += "> (permutations of Fibonacci numbers)";
   for_each(begin(test_functions), end(test_functions), test_function_template<int_type>(fibonacci_numbers_permution_2<T>(), row_name.c_str()));

   row_name.erase();
   row_name += "gcd<";
   row_name += name;
   row_name += "> (Trivial cases)";
   for_each(begin(test_functions), end(test_functions), test_function_template<int_type>(trivial_gcd_test_cases<T>(), row_name.c_str()));
#endif
}

/*******************************************************************************************************************/

template <class T>
T generate_random(unsigned bits_wanted)
{
   static boost::random::mt19937 gen;
   typedef boost::random::mt19937::result_type random_type;

   T max_val;
   unsigned digits;
   if(std::numeric_limits<T>::is_bounded && (bits_wanted == (unsigned)std::numeric_limits<T>::digits))
   {
      max_val = (std::numeric_limits<T>::max)();
      digits = std::numeric_limits<T>::digits;
   }
   else
   {
      max_val = T(1) << bits_wanted;
      digits = bits_wanted;
   }

   unsigned bits_per_r_val = std::numeric_limits<random_type>::digits - 1;
   while((random_type(1) << bits_per_r_val) > (gen.max)()) --bits_per_r_val;

   unsigned terms_needed = digits / bits_per_r_val + 1;

   T val = 0;
   for(unsigned i = 0; i < terms_needed; ++i)
   {
      val *= (gen.max)();
      val += gen();
   }
   val %= max_val;
   return val;
}



int main()
{
    // test_type<int>("polynomial<int>");
    test_type<float>("polynomial<float>");
    test_type<double>("polynomial<double>");
    test_type<mp_float_type>("polynomial<cpp_bin_float_quad>");
}
