//  Copyright Jeremy Murphy 2016.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifdef _MSC_VER
#  pragma warning (disable : 4224)
#endif

#include "../../test/table_type.hpp"
#include "table_helper.hpp"
#include "performance.hpp"

#include <boost/math/common_factor_rt.hpp>

#include <boost/array.hpp>
#include <iostream>
#include <algorithm>
#include <numeric>
#include <string>
#include <tuple>
#include <type_traits>
#include <vector>
#include <functional>

using namespace std;

template <typename Func, typename Result = BOOST_DEDUCED_TYPENAME std::result_of<Func()>::type>
pair<double, Result> exec_timed_test_foo(Func f, double min_elapsed = 0.5)
{
    double t = 0;
    unsigned repeats = 1;
    Result sum{0};
    stopwatch<boost::chrono::high_resolution_clock> w;
    do
    {
        for(unsigned count = 0; count < repeats; ++count)
            sum += f();

        t = boost::chrono::duration_cast<boost::chrono::duration<double>>(w.elapsed()).count();
        if(t < min_elapsed)
            repeats *= 2;
    }
    while(t < min_elapsed);
    return {t / repeats, sum};
}


template <typename T>
struct test_function_template
{
    pair<T, T> data;
    
    test_function_template(pair<T, T> const &data) : data(data) {}
    
    template <typename Function>
    void operator()(pair<Function, string> const &f) const
    {
        auto result = exec_timed_test_foo(bind(f.first, data.first, data.second));
        report_execution_time(result.first, 
                            string("gcd method comparison with ") + compiler_name() + string(" on ") + platform_name(), 
                            f.second, 
                            boost_name());
    }
};


int main()
{
    using namespace boost::math::detail;
    
    typedef unsigned int_type;
    pair<int_type, int_type> test_data{1836311903, 2971215073}; // 46th and 47th Fibonacci numbers. 47th is prime.
    typedef pair< function<int_type(int_type, int_type)>, string> f_test;
    array<f_test, 2> test_functions{{{gcd_euclidean<int_type>, "gcd_euclidean"}, {gcd_Stein<int_type>, "gcd_Stein"}}};
    for_each(begin(test_functions), end(test_functions), test_function_template<int_type>(test_data));
}
