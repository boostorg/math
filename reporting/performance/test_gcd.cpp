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
    void operator()(Function f) const
    {
        auto result = exec_timed_test_foo(bind(f, data.first, data.second));
        cout << f.target_type().name() << ": " << result.first << " (" << result.second << ")" << endl;
        report_execution_time(result.first, 
                            string("gcd method comparison with ") + compiler_name() + string(" on ") + platform_name(), 
                            typeid(decltype(data.first)).name(), 
                            typeid(f).name());
    }
};


int main()
{
    using namespace boost::math::detail;
    
    // pair<unsigned, unsigned> test_data{24140, 40902};
    pair<unsigned, unsigned> test_data{1836311903, 2971215073}; // 46th and 47th Fibonacci numbers. 47th is prime.
    array<function<unsigned(unsigned, unsigned)>, 2> test_functions{{gcd_euclidean<unsigned>, gcd_binary<unsigned>}};
    for_each(begin(test_functions), end(test_functions), test_function_template<unsigned>(test_data));
}
