import sys

def generate_header(data):
    """
    Generates the C++ header file content for autodiff expression templates.

    Args:
        data (list): A list of dictionaries, where each dictionary
                     describes a function, its arguments, and its derivatives.

    Returns:
        str: The complete C++ header file content.
    """
    header_content = """#ifndef REVERSE_MODE_AUTODIFF_GENERATED_H
#define REVERSE_MODE_AUTODIFF_GENERATED_H

#include <boost/math/differentiation/detail/reverse_mode_autodiff_basic_operator_overloads.hpp>
#include <boost/math/differentiation/detail/reverse_mode_autodiff_expression_template_base.hpp>
#include <boost/math/differentiation/detail/reverse_mode_autodiff_stl_overloads.hpp>
#include <boost/math/differentiation/detail/reverse_mode_autodiff_utilities.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <string>

namespace boost {
namespace math {
namespace differentiation {
namespace reverse_mode {

// Forward declarations for all expression structs
"""
    for func in data:
        func_name = func['name']
        num_args = func['num_args']
        if num_args == 1:
            header_content += f"template<typename T, size_t order, typename ARG>\nstruct {func_name}_expr;\n\n"
        elif num_args == 2:
            header_content += f"template<typename T, size_t order, typename LHS, typename RHS>\nstruct {func_name}_expr;\n\n"
            header_content += f"template<typename T, size_t order, typename ARG>\nstruct {func_name}_left_float_expr;\n\n"
            header_content += f"template<typename T, size_t order, typename ARG>\nstruct {func_name}_right_float_expr;\n\n"
    
    header_content += "\n/*******************************************************************************************/\n\n"

    # Function overloads
    for func in data:
        func_name = func['name']
        num_args = func['num_args']

        if num_args == 1:
            header_content += f"""template<typename T, size_t order, typename ARG>
{func_name}_expr<T, order, ARG> {func_name}(const expression<T, order, ARG> &arg)
{{
    return {func_name}_expr<T, order, ARG>(arg, 0.0);
}}
"""
        elif num_args == 2:
            header_content += f"""template<typename T, size_t order, typename LHS, typename RHS>
{func_name}_expr<T, order, LHS, RHS> {func_name}(const expression<T, order, LHS> &left,
                                                     const expression<T, order, RHS> &right)
{{
    return {func_name}_expr<T, order, LHS, RHS>(left, right);
}}
template<typename T, size_t order, typename ARG>
{func_name}_left_float_expr<T, order, ARG> {func_name}(const T                         &left,
                                                          const expression<T, order, ARG> &arg)
{{
    return {func_name}_left_float_expr<T, order, ARG>(arg, left);
}}

template<typename T, size_t order, typename ARG>
{func_name}_right_float_expr<T, order, ARG> {func_name}(const expression<T, order, ARG> &arg,
                                                           const T                         &right)
{{
    return {func_name}_right_float_expr<T, order, ARG>(arg, right);
}}
"""

    header_content += "\n/*******************************************************************************************/\n"
    
    # Expression template structs
    for func in data:
        func_name = func['name']
        num_args = func['num_args']
        
        if num_args == 1:
            header_content += f"""
template<typename T, size_t order, typename ARG>
struct {func_name}_expr : public abstract_unary_expression<T, order, ARG, {func_name}_expr<T, order, ARG>>
{{
    /** @brief {func_name}(x)
     * d/dx {func_name}(x) = {func['derivative_single_arg']} 
      * */
    using arg_type   = ARG;
    using value_type = T;
    using inner_t    = rvar_t<T, order - 1>;

    explicit {func_name}_expr(const expression<T, order, ARG> &arg_expr, const T v)
        : abstract_unary_expression<T, order, ARG, {func_name}_expr<T, order, ARG>>(arg_expr, v){{}};

    inner_t evaluate() const
    {{
        return detail::if_functional_dispatch<(order > 1)>(
            [this](auto &&x) {{ return reverse_mode::{func_name}(std::forward<decltype(x)>(x)); }},
            [this](auto &&x) {{ return boost::math::{func_name}(std::forward<decltype(x)>(x)); }},
            this->arg.evaluate());
    }}
    static const inner_t derivative(const inner_t &argv, const inner_t &v, const T &constant)
    {{
        return detail::if_functional_dispatch<(
            order > 1)>([{func['lambda_capture']}](auto &&x) {{ return {func['derivative_code_single']}; }},
                        [{func['lambda_capture']}](auto &&x) {{ return {func['derivative_boost_single']}; }},
                        argv);
    }}
}};
"""
        elif num_args == 2:
            header_content += f"""
template<typename T, size_t order, typename LHS, typename RHS>
struct {func_name}_expr
    : public abstract_binary_expression<T, order, LHS, RHS, {func_name}_expr<T, order, LHS, RHS>>
{{
    using lhs_type   = LHS;
    using rhs_type   = RHS;
    using value_type = T;
    using inner_t    = rvar_t<T, order - 1>;
    // Explicitly define constructor to forward to base class
    explicit {func_name}_expr(const expression<T, order, LHS> &left_hand_expr,
                               const expression<T, order, RHS> &right_hand_expr)
        : abstract_binary_expression<T, order, LHS, RHS, {func_name}_expr<T, order, LHS, RHS>>(
              left_hand_expr, right_hand_expr)
    {{}}

    inner_t evaluate() const
    {{
        return detail::if_functional_dispatch<(order > 1)>(
            [this](auto &&x, auto &&y) {{
                return reverse_mode::{func_name}(std::forward<decltype(x)>(x),
                                                  std::forward<decltype(y)>(y));
            }},
            [this](auto &&x, auto &&y) {{
                return boost::math::{func_name}(std::forward<decltype(x)>(x),
                                                 std::forward<decltype(y)>(y));
            }},
            this->lhs.evaluate(),
            this->rhs.evaluate());
    }};
    static const inner_t left_derivative(const inner_t &l, const inner_t &r, const inner_t &v)
    {{
        return detail::if_functional_dispatch<(order > 1)>(
            [](auto &&x, auto &&y) {{
                return {func['derivative_code_left']};
            }},
            [](auto &&x, auto &&y) {{
                return {func['derivative_boost_left']};
            }},
            l,
            r);
    }};
    static const inner_t right_derivative(const inner_t &l, const inner_t &r, const inner_t &v)
    {{
        return detail::if_functional_dispatch<(order > 1)>(
            [](auto &&x, auto &&y) {{
                return {func['derivative_code_right']};
            }},
            [](auto &&x, auto &&y) {{
                return {func['derivative_boost_right']};
            }},
            l,
            r);
    }};
}};
"""
            header_content += f"""
template<typename T, size_t order, typename ARG>
struct {func_name}_left_float_expr
    : public abstract_unary_expression<T, order, ARG, {func_name}_left_float_expr<T, order, ARG>>
{{
    using arg_type   = ARG;
    using value_type = T;
    using inner_t    = rvar_t<T, order - 1>;

    explicit {func_name}_left_float_expr(const expression<T, order, ARG> &arg_expr, const T v)
        : abstract_unary_expression<T, order, ARG, {func_name}_left_float_expr<T, order, ARG>>(
              arg_expr, v){{}};

    inner_t evaluate() const
    {{
        return detail::if_functional_dispatch<(order > 1)>(
            [this](auto &&x) {{
                return reverse_mode::{func_name}(this->constant, std::forward<decltype(x)>(x));
            }},
            [this](auto &&x) {{
                return boost::math::{func_name}(this->constant, std::forward<decltype(x)>(x));
            }},
            this->arg.evaluate());
    }}
    static const inner_t derivative(const inner_t &argv, const inner_t &v, const T &constant)
    {{
        return detail::if_functional_dispatch<(order > 1)>(
            [constant](auto &&x) {{
                return {func['derivative_code_right'].replace('y', 'constant')};
            }},
            [constant](auto &&x) {{
                return {func['derivative_boost_right'].replace('y', 'constant')};
            }},
            argv);
    }}
}};

template<typename T, size_t order, typename ARG>
struct {func_name}_right_float_expr
    : public abstract_unary_expression<T, order, ARG, {func_name}_right_float_expr<T, order, ARG>>
{{
    using arg_type   = ARG;
    using value_type = T;
    using inner_t    = rvar_t<T, order - 1>;

    explicit {func_name}_right_float_expr(const expression<T, order, ARG> &arg_expr, const T v)
        : abstract_unary_expression<T, order, ARG, {func_name}_right_float_expr<T, order, ARG>>(
              arg_expr, v){{}};

    inner_t evaluate() const
    {{
        return detail::if_functional_dispatch<(order > 1)>(
            [this](auto &&x) {{
                return reverse_mode::{func_name}(std::forward<decltype(x)>(x), this->constant);
            }},
            [this](auto &&x) {{
                return boost::math::{func_name}(std::forward<decltype(x)>(x), this->constant);
            }},
            this->arg.evaluate());
    }}
    static const inner_t derivative(const inner_t &argv, const inner_t &v, const T &constant)
    {{
        return detail::if_functional_dispatch<(order > 1)>(
            [constant](auto &&x) {{
                return {func['derivative_code_left'].replace('y', 'constant')};
            }},
            [constant](auto &&x) {{
                return {func['derivative_boost_left'].replace('y', 'constant')};
            }},
            argv);
    }}
}};
"""

    header_content += """
} // namespace reverse_mode
} // namespace differentiation
} // namespace math
} // namespace boost

#endif
"""
    return header_content

def main():
    """
    Main function to define the functions to be generated and call the generator.
    """
    # A list of dictionaries defining the functions. Add more functions here to generate more code.
    functions_to_generate = [
        {
            "name": "tgamma",
            "num_args": 1,
            "derivative_single_arg": "gamma(x)*digamma(x)",
            "lambda_capture": "",
            "derivative_code_single": "reverse_mode::tgamma(x) * reverse_mode::digamma(x)",
            "derivative_boost_single": "boost::math::tgamma(x) * boost::math::digamma(x)"
        },
        {
            "name": "lgamma",
            "num_args": 1,
            "derivative_single_arg": "digamma(x)",
            "lambda_capture": "",
            "derivative_code_single": "reverse_mode::digamma(x)",
            "derivative_boost_single": "boost::math::digamma(x)"
        },
        {
            "name": "digamma",
            "num_args": 1,
            "derivative_single_arg": "polygamma(1,x)",
            "lambda_capture": "",
            "derivative_code_single": "reverse_mode::polygamma(1, x)",
            "derivative_boost_single": "boost::math::polygamma(1, x)"
        },
        {
            "name": "trigamma",
            "num_args": 1,
            "derivative_single_arg": "polygamma(2,x)",
            "lambda_capture": "",
            "derivative_code_single": "reverse_mode::polygamma(2, x)",
            "derivative_boost_single": "boost::math::polygamma(2, x)"
        },
        {
            "name": "polygamma",
            "num_args": 1,
            "derivative_single_arg": "polygamma(n+1,x)",
            "lambda_capture": "constant",
            "derivative_code_single": "reverse_mode::polygamma(static_cast<int>(constant) + 1, x)",
            "derivative_boost_single": "boost::math::polygamma(static_cast<int>(constant) + 1, x)"
        },
        {
            "name": "tgamma_ratio",
            "num_args": 2,
            "derivative_code_left": "reverse_mode::tgamma_ratio(x, y) * reverse_mode::digamma(x)",
            "derivative_boost_left": "boost::math::tgamma_ratio(x, y) * boost::math::digamma(x)",
            "derivative_code_right": "-reverse_mode::tgamma_ratio(x, y) * reverse_mode::digamma(y)",
            "derivative_boost_right": "-boost::math::tgamma_ratio(x, y) * boost::math::digamma(y)"
        },
        {
            "name": "tgamma_delta_ratio",
            "num_args": 2,
            "derivative_code_left": "reverse_mode::tgamma_delta_ratio(x, y) * (reverse_mode::digamma(x) - reverse_mode::digamma(x + y))",
            "derivative_boost_left": "boost::math::tgamma_delta_ratio(x, y) * (boost::math::digamma(x) - boost::math::digamma(x + y))",
            "derivative_code_right": "-reverse_mode::tgamma_delta_ratio(x, y) * (reverse_mode::digamma(x + y))",
            "derivative_boost_right": "-boost::math::tgamma_delta_ratio(x, y) * (boost::math::digamma(x + y))"
        },
        {
            "name": "gamma_p",
            "num_args": 2,
            "derivative_code_left": "reverse_mode::tgamma_delta_ratio(x, y) * (reverse_mode::digamma(x) - reverse_mode::digamma(x + y))",
            "derivative_boost_left": "boost::math::tgamma_delta_ratio(x, y) * (boost::math::digamma(x) - boost::math::digamma(x + y))",
            "derivative_code_right": "-reverse_mode::tgamma_delta_ratio(x, y) * (reverse_mode::digamma(x + y))",
            "derivative_boost_right": "-boost::math::tgamma_delta_ratio(x, y) * (boost::math::digamma(x + y))"
        }
    ]

    generated_code = generate_header(functions_to_generate)
    print(generated_code)

if __name__ == "__main__":
    main()
