//           Copyright Matthew Pulver 2018 - 2019.
// Distributed under the Boost Software License, Version 1.0.
//      (See accompanying file LICENSE_1_0.txt or copy at
//           https://www.boost.org/LICENSE_1_0.txt)

/*!
\mainpage

\section synopsis Synopsis
\dontinclude synopsis.cpp \skip #include
\until **

\section description Description

Autodiff is a header-only C++ library that facilitates the
[automatic differentiation](https://en.wikipedia.org/wiki/Automatic_differentiation) (forward mode)
of mathematical functions in single and multiple variables.

The formula central to this implementation of automatic differentiation is the following [Taylor
series](https://en.wikipedia.org/wiki/Taylor_series) expansion of an analytic function \f$f\f$ at the point \f$x_0\f$:
\f{align*}
f(x_0+\varepsilon) &= f(x_0) + f'(x_0)\varepsilon + \frac{f''(x_0)}{2!}\varepsilon^2 + \frac{f'''(x_0)}{3!}\varepsilon^3 + \cdots \\
  &= \sum_{n=0}^N\frac{f^{(n)}(x_0)}{n!}\varepsilon^n + O\left(\varepsilon^{N+1}\right).
\f}

The essential idea of autodiff is the substitution of numbers with polynomials in the evaluation of \f$f\f$. By
selecting the proper polynomial \f$x_0+\varepsilon\f$ as input, the resulting polynomial contains the function's
derivatives within the polynomial coefficients. One then multiplies by a factorial term to obtain the desired
derivative of any order.

Assume one is interested in the first \f$N\f$ derivatives of \f$f\f$ at \f$x_0\f$. Then without any loss of
precision to the calculation of the derivatives, all terms \f$O\left(\varepsilon^{N+1}\right)\f$ that include
powers of \f$\varepsilon\f$ greater than \f$N\f$ can be discarded, and under these truncation rules, \f$f\f$
provides a polynomial-to-polynomial transformation:
\f[
f\quad:\quad x_0+\varepsilon\quad\mapsto\quad\sum_{n=0}^N\frac{f^{(n)}(x_0)}{n!}\varepsilon^n.
\f]

C++ includes the ability to overload operators and functions, and thus when \f$f\f$ is written as a template
function that can receive and return a generic type, then that is sufficient to perform automatic differentiation:
Create a class that models polynomials, and overload all of the arithmetic operators to model polynomial arithmetic
that drop all terms in \f$O\left(\varepsilon^{N+1}\right)\f$. The derivatives are then found in the coefficients of
the return value. This is essentially what the autodiff library does (generalizing to multiple independent variables.)

\subsection requirements Requirements

1. C++11 compiler, but optimized for C++17. Visual Studio 2015 is not supported.
2. Maximum derivative orders are set at compile-time. This allows for compile-time calculation of memory
   requirements for use with `std::array<>`, so that use of dynamic memory allocation is avoided.
3. Mathematical functions, whose derivatives are desired, should accept generic types (template variables) for
   the parameters that derivatives are calculated with respect to, and internal function calls should allow for
   [argument-dependent lookup](https://en.cppreference.com/w/cpp/language/adl) (ADL) on these parameters.

\subsection quick-install Quick Install

Since this is a header-only library, there is an easy way to use/test this library without having
to install it system-wide or compile libraries ahead of time:

1. Download the latest version of boost from https://www.boost.org/ (minimum 1.70.0) and uncompress it.
2. A single directory is created, for example, `./boost_1_70_0`.
3. Add the new directory to your include path when compiling. Example: `g++ -I./boost_1_70_0 ./boost_1_70_0/libs/math/example/autodiff_fourth_power.cpp`

Long term it is recommended to install boost via a package manager specific to your operating system, or with
the `INSTALL` directions provided in the download.

\section examples Examples

\subsection example-single-variable Example 1: Single-variable derivatives
<h3>Calculate derivatives of \f$f(x)=x^4\f$ at \f$x=2.0\f$.</h3>

In this example, `autodiff::variable<double,5>` is a data type that can hold a polynomial of up to degree 5,
and the initialization `autodiff::variable<double,5> x(2.0)` represents the polynomial \f$2.0 + \varepsilon\f$.
Internally, this is modeled by a `std::array<double,6>` whose elements correspond to the 6 coefficients of the
polynomial: `{2, 1, 0, 0, 0, 0}` upon initialization.

\dontinclude fourth_power.cpp \skip #include
\until **

The above calculates \f{alignat*}{{3}
{\tt y.derivative(0)} &=& f(2) =&& \left.x^4\right|_{x=2} &= 16\\
{\tt y.derivative(1)} &=& f'(2) =&& \left.4\cdot x^3\right|_{x=2} &= 32\\
{\tt y.derivative(2)} &=& f''(2) =&& \left.4\cdot 3\cdot x^2\right|_{x=2} &= 48\\
{\tt y.derivative(3)} &=& f'''(2) =&& \left.4\cdot 3\cdot2\cdot x\right|_{x=2} &= 48\\
{\tt y.derivative(4)} &=& f^{(4)}(2) =&& 4\cdot 3\cdot2\cdot1 &= 24\\
{\tt y.derivative(5)} &=& f^{(5)}(2) =&& 0 &
\f}

<hr>

\subsection example-multiprecision Example 2: Multi-variable mixed partial derivatives with multi-precision data type
<h3>Calculate \f$\frac{\partial^{12}f}{\partial w^{3}\partial x^{2}\partial y^{4}\partial z^{3}}(11,12,13,14)\f$
with a precision of about 100 decimal digits, where
\f$f(w,x,y,z)=\exp\left(w\sin\left(\frac{x\log(y)}{z}\right)+\sqrt{\frac{wz}{xy}}\right)+\frac{w^2}{\tan(z)}\f$.</h3>

In this example, the data type `autodiff::variable<cpp_dec_float_100,Nw,Nx,Ny,Nz>` represents a polynomial in 4
independent variables, where the highest powers of each are `Nw`, `Nx`, `Ny` and `Nz`. The underlying arithmetic
data type, referred to as `root_type`, is `boost::multiprecision::cpp_dec_float_100`. The internal data type is
`std::array<std::array<std::array<std::array<cpp_dec_float_100,Nz+1>,Ny+1>,Nx+1>,Nw+1>`. In general, the `root_type`
is always the first template parameter to `autodiff::variable<>` followed by the maximum derivative order that is
to be calculated for each independent variable.

Note that when variables are initialized, the position of the last derivative order given in the parameter pack
determines which variable is taken to be independent. In other words, it determines which of the 4 different
polynomial variables \f$\varepsilon_w,\varepsilon_x,\varepsilon_y,\f$ or \f$\varepsilon_z\f$ are to be added:
\f{align*}
\texttt{autodiff::variable<cpp_dec_float_100,Nw>(11)} &= 11+\varepsilon_w \\
\texttt{autodiff::variable<cpp_dec_float_100,0,Nx>(12)} &= 12+\varepsilon_x \\
\texttt{autodiff::variable<cpp_dec_float_100,0,0,Ny>(13)} &= 13+\varepsilon_y \\
\texttt{autodiff::variable<cpp_dec_float_100,0,0,0,Nz>(14)} &= 14+\varepsilon_z
\f}

Instances of different types are automatically promoted to the smallest multi-variable type that accommodates
both when they are arithmetically combined (added, subtracted, multiplied, divided.)

\dontinclude multiprecision.cpp \skip #include
\until **

<hr>

\subsection example-black-scholes Example 3: Black-Scholes option pricing
<h3>Using the standard Black-Scholes model for pricing European options, calculate call/put prices and greeks.</h3>
https://en.wikipedia.org/wiki/Greeks_(finance)#Formulas_for_European_option_Greeks

One of the primary benefits of using automatic differentiation is the elimination of additional functions to
calculate derivatives, which is a form of code redundancy.
\dontinclude black_scholes.cpp \skip #include
\until **

<hr>

\subsection example-multi-variable Example 4: Multi-variable mixed partial derivatives
<h3>Calculate mixed partial derivatives of
\f$f(w,x,y,z)=\exp\left(w\sin\left(\frac{x\log(y)}{z}\right)+\sqrt{\frac{wz}{xy}}\right)+\frac{w^2}{\tan(z)}\f$
at \f$(w,x,y,z)=(11,12,13,14)\f$.</h3>

This example calculates
\f[
{\tt v.derivative(iw,ix,iy,iz)} =
\frac{\partial^{i_w+i_x+i_y+i_z}f}{\partial w^{i_w}\partial x^{i_x}\partial y^{i_y}\partial z^{i_z}}(11,12,13,14)
\qquad\text{for each}\qquad (i_w,i_x,i_y,i_z)\in \mathbb{N}_4\times\mathbb{N}_3\times\mathbb{N}_5\times\mathbb{N}_4
\f]

where \f$\mathbb{N}_i=\{0,1,2,...,i-1\}\f$.
For testing purposes, the \f$4\times3\times5\times4=240\f$-element `answers[]` array was calculated independently
by Mathematica in 2 steps:

 1. <a href="https://reference.wolfram.com/language/tutorial/SymbolicComputation.html">Symbolic differentiation</a>
 2. <a href="https://reference.wolfram.com/language/tutorial/ArbitraryPrecisionNumbers.html">Arbitrary-precision computation</a> to calculate the answers to 20 decimal places.

Out of the 240 calculated values, the maximum relative error between the values calculated by Mathematica vs.
the Boost Autodiff library is found to be about \f$6.82\times10^{-13}\f$ using the standard IEEE double precision
floating point data type. Since the data type is a template variable, the error can be reduced arbitrarily by
using a data type with greater precision.

\dontinclude mixed_partials.cpp \skip #include
\until **

<hr>

\section mathematics Mathematics

In order for the usage of the autodiff library to make sense, a basic understanding of the mathematics will help.

\subsection taylor-series Truncated Taylor Series

Basic calculus courses teach that a real [analytic function](https://en.wikipedia.org/wiki/Analytic_function)
\f$f : D\rightarrow\mathbb{R}\f$ is one which can be expressed as a Taylor series at a point
\f$x_0\in D\subseteq\mathbb{R}\f$:

\f[
f(x) = f(x_0) + f'(x_0)(x-x_0) + \frac{f''(x_0)}{2!}(x-x_0)^2 + \frac{f'''(x_0)}{3!}(x-x_0)^3 + \cdots \f]

One way of thinking about this form is that given the value of an analytic function \f$f(x_0)\f$ and its derivatives
\f$f'(x_0), f''(x_0), f'''(x_0), ...\f$ evaluated at a point \f$x_0\f$, then the value of the function
\f$f(x)\f$ can be obtained at any other point \f$x\in D\f$ using the above formula.

Let us make the substitution \f$x=x_0+\varepsilon\f$ and rewrite the above equation to get:

\f[
f(x_0+\varepsilon) = f(x_0) + f'(x_0)\varepsilon + \frac{f''(x_0)}{2!}\varepsilon^2 + \frac{f'''(x_0)}{3!}\varepsilon^3 + \cdots
\f]

Now consider \f$\varepsilon\f$ as **an abstract algebraic entity that never acquires a numeric value**, much like
one does in basic algebra with variables like \f$x\f$ or \f$y\f$. For example, we can still manipulate entities
like \f$xy\f$ and \f$(1+2x+3x^2)\f$ without having to assign specific numbers to them.

Using this formula, autodiff goes in the other direction. Given a general formula/algorithm for calculating
\f$f(x_0+\varepsilon)\f$, the derivatives are obtained from the coefficients of the powers of \f$\varepsilon\f$
in the resulting computation. The general coefficient for \f$\varepsilon^n\f$ is

\f[\frac{f^{(n)}(x_0)}{n!}.\f]

Thus to obtain \f$f^{(n)}(x_0)\f$, the coefficient of \f$\varepsilon^n\f$ is multiplied by \f$n!\f$.

\subsubsection taylor-series-example Example

Apply the above technique to calculate the derivatives of \f$f(x)=x^4\f$ at \f$x_0=2\f$.

The first step is to evaluate \f$f(x_0+\varepsilon)\f$ and simply go through the calculation/algorithm, treating
\f$\varepsilon\f$ as an abstract algebraic entity:

\f{align*}
f(x_0+\varepsilon) &= f(2+\varepsilon) \\
 &= (2+\varepsilon)^4 \\
 &= \left(4+4\varepsilon+\varepsilon^2\right)^2 \\
 &= 16+32\varepsilon+24\varepsilon^2+8\varepsilon^3+\varepsilon^4.
\f}

Equating the powers of \f$\varepsilon\f$ from this result with the above \f$\varepsilon\f$-taylor expansion
yields the following equalities:
\f[
f(2) = 16, \qquad
f'(2) = 32, \qquad
\frac{f''(2)}{2!} = 24, \qquad
\frac{f'''(2)}{3!} = 8, \qquad
\frac{f^{(4)}(2)}{4!} = 1, \qquad
\frac{f^{(5)}(2)}{5!} = 0.
\f]

Multiplying both sides by the respective factorials gives
\f[
f(2) = 16, \qquad
f'(2) = 32, \qquad
f''(2) = 48, \qquad
f'''(2) = 48, \qquad
f^{(4)}(2) = 24, \qquad
f^{(5)}(2) = 0.
\f]

These values can be directly confirmed by the [power rule](https://en.wikipedia.org/wiki/Power_rule)
applied to \f$f(x)=x^4\f$.

<hr>

\subsection arithmetic Arithmetic

What was essentially done above was to take a formula/algorithm for calculating \f$f(x_0)\f$ from a number \f$x_0\f$,
and instead apply the same formula/algorithm to a polynomial \f$x_0+\varepsilon\f$. Intermediate steps operate on
values of the form

\f[ {\bf x} = x_0 + x_1\varepsilon + x_2\varepsilon^2 +\cdots+ x_N\varepsilon^N \f]

and the final return value is of this polynomial form as well. In other words, the normal arithmetic operators
\f$+,-,\times,\div\f$ applied to numbers \f$x\f$ are instead applied to polynomials \f$\bf x\f$. Through the
overloading of C++ operators and functions, floating point data types are replaced with data types that represent
these polynomials. More specifically, C++ types such as `double` are replaced with `std::array<double,N+1>`, which
hold the above \f$N+1\f$ coefficients \f$x_i\f$, and are wrapped in a `class` that overloads all of the arithmetic
operators.

The logic of these arithmetic operators simply mirror that which is applied to polynomials. We'll look at
each of the 4 arithmetic operators in detail.

\subsubsection arithmetic-addition Addition

Given polynomials \f$\bf x\f$ and \f$\bf y\f$, how is \f$\bf z=x+y\f$ calculated?

To answer this, one simply expands \f$\bf x\f$ and \f$\bf y\f$ into their polynomial forms and add them together:

\f{align*}
{\bf z} &= {\bf x} + {\bf y} \\
 &= \left(\sum_{i=0}^Nx_i\varepsilon^i\right) + \left(\sum_{i=0}^Ny_i\varepsilon^i\right) \\
 &= \sum_{i=0}^N(x_i+y_i)\varepsilon^i \\
z_i &= x_i + y_i \qquad \text{for}\; i\in\{0,1,2,...,N\}.
\f}

\subsubsection arithmetic-subtraction Subtraction

Subtraction follows the same form as addition:

\f{align*}
{\bf z} &= {\bf x} - {\bf y} \\
 &= \left(\sum_{i=0}^Nx_i\varepsilon^i\right) - \left(\sum_{i=0}^Ny_i\varepsilon^i\right) \\
 &= \sum_{i=0}^N(x_i-y_i)\varepsilon^i \\
z_i &= x_i - y_i \qquad \text{for}\; i\in\{0,1,2,...,N\}.
\f}

\subsubsection arithmetic-multiplication Multiplication

Multiplication is a bit more interesting:

\f{align*}
{\bf z} &= {\bf x} \times {\bf y} \\
 &= \left(\sum_{i=0}^Nx_i\varepsilon^i\right) \left(\sum_{i=0}^Ny_i\varepsilon^i\right) \\
 &= x_0y_0 + (x_0y_1+x_1y_0)\varepsilon + (x_0y_2+x_1y_1+x_2y_0)\varepsilon^2 + \cdots +
    \left(\sum_{j=0}^Nx_jy_{N-j}\right)\varepsilon^N + O\left(\varepsilon^{N+1}\right) \\
 &= \sum_{i=0}^N\sum_{j=0}^ix_jy_{i-j}\varepsilon^i + O\left(\varepsilon^{N+1}\right) \\
z_i &= \sum_{j=0}^ix_jy_{i-j} \qquad \text{for}\; i\in\{0,1,2,...,N\}.
\f}

In the case of multiplication, terms involving powers of \f$\varepsilon\f$ greater than \f$N\f$, collectively
denoted by \f$O\left(\varepsilon^{N+1}\right)\f$, are simply discarded. Fortunately, the values of \f$z_i\f$
for \f$i\le N\f$ do not depend on any of these discarded terms, so there is no loss of precision in the final
answer. The only information that is lost are the values of higher order derivatives, which we are not interested
in anyway. If we were, then we would have simply chosen a larger value of \f$N\f$ to begin with.

\subsubsection arithmetic-division Division

Division is not directly calculated as are the others. Instead, to find the components of
\f${\bf z}={\bf x}\div{\bf y}\f$ we require that \f${\bf x}={\bf y}\times{\bf z}\f$. This yields
a recursive formula for the components \f$z_i\f$:

\f{align*}
x_i &= \sum_{j=0}^iy_jz_{i-j} \\
 &= y_0z_i + \sum_{j=1}^iy_jz_{i-j} \\
z_i &= \frac{1}{y_0}\left(x_i - \sum_{j=1}^iy_jz_{i-j}\right) \qquad \text{for}\; i\in\{0,1,2,...,N\}.
\f}

In the case of division, the values for \f$z_i\f$ must be calculated sequentially, since \f$z_i\f$
depends on the previously calculated values \f$z_0, z_1, ..., z_{i-1}\f$.

<hr>

\subsection general-functions General Functions

When calling standard mathematical functions such as `log()`, `cos()`, etc. how should these be written in order
to support autodiff variable types? That is, how should they be written to provide accurate derivatives?

To simplify notation, for a given polynomial \f${\bf x} = x_0 + x_1\varepsilon + x_2\varepsilon^2 +\cdots+
x_N\varepsilon^N\f$ define

\f[
{\bf x}_\varepsilon = x_1\varepsilon + x_2\varepsilon^2 +\cdots+ x_N\varepsilon^N = \sum_{i=1}^Nx_i\varepsilon^i
\f]

This allows for a concise expression of a general function \f$f\f$ of \f$\bf x\f$:

\f{align*}
f({\bf x}) &= f(x_0 + {\bf x}_\varepsilon) \\
 & = f(x_0) + f'(x_0){\bf x}_\varepsilon + \frac{f''(x_0)}{2!}{\bf x}_\varepsilon^2 + \frac{f'''(x_0)}{3!}{\bf x}_\varepsilon^3 + \cdots + \frac{f^{(N)}(x_0)}{N!}{\bf x}_\varepsilon^N + O\left(\varepsilon^{N+1}\right) \\
 & = \sum_{i=0}^N\frac{f^{(i)}(x_0)}{i!}{\bf x}_\varepsilon^i + O\left(\varepsilon^{N+1}\right)
\f}

where \f$\varepsilon\f$ has been substituted with \f${\bf x}_\varepsilon\f$ in the \f$\varepsilon\f$-taylor series
for \f$f(x)\f$. This form gives a recipe for calculating \f$f({\bf x})\f$ in general from regular numeric calculations
\f$f(x_0)\f$, \f$f'(x_0)\f$, \f$f''(x_0)\f$, ... and successive powers of the epsilon terms \f${\bf x}_\varepsilon\f$.

For an application in which we are interested in up to \f$N\f$ derivatives in \f$x\f$ the data structure to hold
this information is an \f$(N+1)\f$-element array `v` whose general element is

\f[ {\tt v[i]} = \frac{f^{(i)}(x_0)}{i!} \qquad \text{for}\; i\in\{0,1,2,...,N\}. \f]

<hr>

\subsection multiple-variables Multiple Variables

In C++, the generalization to mixed partial derivatives with multiple independent variables is conveniently achieved
with recursion. To begin to see the recursive pattern, consider a two-variable function \f$f(x,y)\f$. Since \f$x\f$
and \f$y\f$ are independent, they require their own independent epsilons \f$\varepsilon_x\f$ and \f$\varepsilon_y\f$,
respectively.

Expand \f$f(x,y)\f$ for \f$x=x_0+\varepsilon_x\f$:
\f{align*}
f(x_0+\varepsilon_x,y) &= f(x_0,y)
+ \frac{\partial f}{\partial x}(x_0,y)\varepsilon_x
+ \frac{1}{2!}\frac{\partial^2 f}{\partial x^2}(x_0,y)\varepsilon_x^2
+ \frac{1}{3!}\frac{\partial^3 f}{\partial x^3}(x_0,y)\varepsilon_x^3
+ \cdots
+ \frac{1}{M!}\frac{\partial^M f}{\partial x^M}(x_0,y)\varepsilon_x^M
+ O\left(\varepsilon_x^{M+1}\right) \\
&= \sum_{i=0}^M\frac{1}{i!}\frac{\partial^i f}{\partial x^i}(x_0,y)\varepsilon_x^i + O\left(\varepsilon_x^{M+1}\right).
\f}

Next, expand \f$f(x_0+\varepsilon_x,y)\f$ for \f$y=y_0+\varepsilon_y\f$:
\f{align*}
f(x_0+\varepsilon_x,y_0+\varepsilon_y) &= \sum_{j=0}^N\frac{1}{j!}\frac{\partial^j}{\partial y^j}
    \left(\sum_{i=0}^M\frac{1}{i!}\frac{\partial^if}{\partial x^i}\right)(x_0,y_0)\varepsilon_x^i\varepsilon_y^j
    + O\left(\varepsilon_x^{M+1}\right) + O\left(\varepsilon_y^{N+1}\right) \\
&= \sum_{i=0}^M\sum_{j=0}^N\frac{1}{i!j!}\frac{\partial^{i+j}f}{\partial x^i\partial y^j}(x_0,y_0)
   \varepsilon_x^i\varepsilon_y^j + O\left(\varepsilon_x^{M+1}\right) + O\left(\varepsilon_y^{N+1}\right).
\f}

Similarly to the single-variable case, for an application in which we are interested in up to \f$M\f$ derivatives in
\f$x\f$ and \f$N\f$ derivatives in \f$y\f$, the data structure to hold this information is an \f$(M+1)\times(N+1)\f$
array `v` whose element at \f$(i,j)\f$ is

\f[ {\tt v[i][j]} = \frac{1}{i!j!}\frac{\partial^{i+j}f}{\partial x^i\partial y^j}(x_0,y_0)
    \qquad \text{for}\; (i,j)\in\{0,1,2,...,M\}\times\{0,1,2,...,N\}. \f]

The generalization to additional independent variables follows the same pattern. This is made more concrete with
C++ code in the next section.

<hr>

\section Usage

\subsection usage-single-variable Single Variable

To calculate derivatives of a single variable \f$x\f$, at a particular value \f$x_0\f$, the following must be
specified at compile-time:

1. The numeric data type `T` of \f$x_0\f$. Examples: `double`, `boost::multiprecision::cpp_dec_float_100`, etc.
2. The **maximum** derivative order \f$M\f$ that is to be calculated with respect to \f$x\f$.

Note that both of these requirements are entirely analogous to declaring and using a `std::array<T,N>`. `T` and
`N` must be set as compile-time, but which elements in the array are accessed can be determined at run-time,
just as the choice of what derivatives to query in autodiff can be made during run-time.

To declare and initialize \f$x\f$:

\code{.cpp}
boost::math::differentiation::autodiff::variable<T,M> x(x0);
\endcode

where `x0` is a run-time value of type `T`. Assuming `0 < M`, this represents the polynomial \f$ x_0 + \varepsilon
\f$. Internally, the member variable of type `std::array<T,M>` is `v = { x0, 1, 0, 0, ... }`, consistent with the
above mathematical treatise.

To find the derivatives \f$f^{(n)}(x_0)\f$ for \f$0\le n\le M\f$ of a function
\f$f : \mathbb{R}\rightarrow\mathbb{R}\f$, the function can be represented as a template

\code{.cpp}
template<typename T>
T f(T x);
\endcode

Using a generic type `T` allows for `x` to be of a regular type such as `double`, but also allows for
`boost::math::differentiation::autodiff::variable<>` types.

Internal calls to mathematical functions must allow for [argument-dependent
lookup](https://en.cppreference.com/w/cpp/language/adl) (ADL). Many standard library functions are overloaded in
the `boost::math::differentiation::autodiff` namespace. For example, instead of calling `std::cos(x)` from within
`f`, include the line `using std::cos;` and call `cos(x)` without a namespace prefix.

Calling \f$f\f$ and retrieving the calculated value and derivatives:

\code{.cpp}
boost::math::differentiation::autodiff::variable<T,M> x(x0);
boost::math::differentiation::autodiff::variable<T,M> y = f(x);
for (int n=0 ; n<=M ; ++n)
    std::cout << "y.derivative("<<n<<") == " << y.derivative(n) << std::endl;
\endcode

`y.derivative(0)` returns the undifferentiated value \f$f(x_0)\f$, and `y.derivative(n)` returns \f$f^{(n)}(x_0)\f$.
Casting `y` to type `T` also gives the undifferentiated value. In other words, the following 3 values are equal:
1. `f(x0)`
2. `y.derivative(0)`
3. `static_cast<T>(y)`

\subsection usage-multiple-variables Multiple Variables

Independent variables are represented in autodiff as independent dimensions within a multi-dimensional array.
This is perhaps best illustrated with examples.

The following instantiates a variable of \f$x=13\f$ with up to 3 orders of derivatives:
\code{.cpp}
boost::math::differentiation::autodiff::variable<double,3> x(13);
\endcode

This instantiates *an independent* value of \f$y=14\f$ with up to 4 orders of derivatives:
\code{.cpp}
boost::math::differentiation::autodiff::variable<double,0,4> y(14);
\endcode

Combining them together *promotes* their data type automatically to the smallest multidimensional array that
accommodates both.
\code{.cpp}
auto z = 10*x*x + 50*x*y + 100*y*y; // z is promoted to boost::math::differentiation::autodiff::variable<double,3,4>
\endcode

The object `z` holds a 2-dimensional array, thus `derivative(...)` is a 2-parameter method:

\f[
{\tt z.derivative(i,j)} = \frac{\partial^{i+j}f}{\partial x^i\partial y^j}(13,14)
    \qquad \text{for}\; (i,j)\in\{0,1,2,3\}\times\{0,1,2,3,4\}.
\f]

A few values of the result can be confirmed through inspection:
\code{.cpp}
z.derivative(2,0) == 20
z.derivative(1,1) == 50
z.derivative(0,2) == 200
\endcode

Note how the position of the parameters in `derivative(..)` match how `x` and `y` were declared.
This will be clarified next.

\subsubsection two-rules Two Rules of Variable Initialization

In general, there are two rules to keep in mind when dealing with multiple variables:

1. Independent variables correspond to parameter position, in both the declaration
   `boost::math::differentiation::autodiff::variable<T,...>` and calls to `derivative(...)`.
2. The last template position in a value-initialized instance determines which variable a derivative will be
   taken with respect to.

Both rules are illustrated with an example in which there are 3 independent variables \f$x,y,z\f$ and 1 dependent
variable \f$w=f(x,y,z)\f$, though the following code readily generalizes to any number of independent variables,
limited only by the C++ compiler/memory/platform. The maximum derivative order of each variable is `Nx`, `Ny`,
and `Nz`, respectively. Then the type for `w` is `boost::math::differentiation::autodiff::variable<T,Nx,Ny,Nz>` and
all possible mixed partial derivatives are available via
\f[
{\tt w.derivative(nx,ny,nz)} =
    \frac{\partial^{n_x+n_y+n_z}f}{\partial x^{n_x}\partial y^{n_y}\partial z^{n_z} }(x_0,y_0,z_0)
    \qquad \text{for}\; (n_x,n_y,n_z)\in\{0,1,2,...,N_x\}\times\{0,1,2,...,N_y\}\times\{0,1,2,...,N_z\}
\f]

where \f$x_0, y_0, z_0\f$ are the numerical values at which the function \f$f\f$ and its derivatives are evaluated.

In code:
\code{.cpp}
using namespace boost::math::differentiation;

using var = autodiff::variable<double,Nx,Ny,Nz>; // Nx, Ny, Nz are constexpr size_t.

var x = autodiff::variable<double,Nx>(x0);       // x0 is of type double
var y = autodiff::variable<double,Nx,Ny>(y0);    // y0 is of type double
var z = autodiff::variable<double,Nx,Ny,Nz>(z0); // z0 is of type double

var w = f(x,y,z);

for (size_t nx=0 ; nx<=Nx ; ++nx)
    for (size_t ny=0 ; ny<=Ny ; ++ny)
        for (size_t nz=0 ; nz<=Nz ; ++nz)
            std::cout << "w.derivative("<<nx<<','<<ny<<','<<nz<<") == " << w.derivative(nx,ny,nz) << std::endl;
\endcode

Note how `x`, `y`, and `z` are initialized: the last template parameter determines which variable
\f$x, y,\f$ or \f$z\f$ a derivative is taken with respect to. In terms of the \f$\varepsilon\f$-polynomials
above, this determines whether to add \f$\varepsilon_x, \varepsilon_y,\f$ or \f$\varepsilon_z\f$ to
\f$x_0, y_0,\f$ or \f$z_0\f$, respectively.

In contrast, the following initialization of `x` would be INCORRECT:
\code{.cpp}
var x = autodiff::variable<T,Nx,0>(x0); // WRONG
\endcode

Mathematically, this represents \f$x_0+\varepsilon_y\f$, since the last template parameter corresponds to the
\f$y\f$ variable, and thus the resulting value will be invalid.

\subsubsection type-promotion Type Promotion

The previous example can be optimized to save some unnecessary computation, by declaring smaller arrays,
and relying on autodiff's automatic type-promotion:
\code{.cpp}
using namespace boost::math::differentiation;

autodiff::variable<double,Nx> x(x0);
autodiff::variable<double,0,Ny> y(y0);
autodiff::variable<double,0,0,Nz> z(z0);

autodiff::variable<double,Nx,Ny,Nz> w = f(x,y,z);

for (size_t nx=0 ; nx<=Nx ; ++nx)
    for (size_t ny=0 ; ny<=Ny ; ++ny)
        for (size_t nz=0 ; nz<=Nz ; ++nz)
            std::cout << "w.derivative("<<nx<<','<<ny<<','<<nz<<") == " << w.derivative(nx,ny,nz) << std::endl;
\endcode

For example, if one of the first steps in the computation of \f$f\f$ was `z*z`, then a significantly less
number of multiplications and additions may occur if `z` is declared as `variable<double,0,0,Nz>` as opposed to
`variable<double,Nx,Ny,Nz>`. There is no loss of precision with the former, since the extra dimensions represent
0 values. Once `z` is combined with `x` and `y` during the computation, the types will be promoted as necessary.
This is the recommended way to initialize variables in autodiff.

\section acknowledgments Acknowledgments

- Kedar Bhat - C++11 compatibility, codecov integration, and feedback.
- Nick Thompson - Initial feedback and help with Boost integration.
- John Maddock - Initial feedback and help with Boost integration.

\copyright
<center>
       Copyright &copy; Matthew Pulver 2018 - 2019.<br/>
Distributed under the Boost Software License, Version 1.0.<br/>
    (See accompanying file LICENSE_1_0.txt or copy at<br/>
          https://www.boost.org/LICENSE_1_0.txt)
</center>
*/

#ifndef BOOST_MATH_AUTODIFF_HPP
#define BOOST_MATH_AUTODIFF_HPP

#include <boost/config.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/factorials.hpp>
#include <boost/math/tools/promotion.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <initializer_list>
#include <limits>
#include <numeric>
#include <ostream>
#include <type_traits>

// Automatic Differentiation v1
namespace boost { namespace math { namespace differentiation { namespace autodiff { inline namespace v1 {

// Use variable<> instead of dimension<> or nested_dimensions<>.
template<typename RealType,size_t Order>
class dimension;
template<typename RealType,size_t Order,size_t... Orders> // specialized for dimension<> below.
struct nested_dimensions { using type = dimension<typename nested_dimensions<RealType,Orders...>::type,Order>; };

// The variable<> alias is the primary template for instantiating autodiff variables.
// This nests one or more dimension<> classes together into a multi-dimensional type.
template<typename RealType,size_t Order,size_t... Orders>
using variable = typename nested_dimensions<RealType,Order,Orders...>::type;

////////////////////////////////////////////////////////////////////////////////

template<typename RealType0,typename RealType1,typename... RealTypes>
struct promote_args_n { using type = typename boost::math::tools::promote_args_2<RealType0,
    typename promote_args_n<RealType1,RealTypes...>::type>::type; };

template<typename RealType0,typename RealType1,typename... RealTypes>
using promote = typename promote_args_n<RealType0,RealType1,RealTypes...>::type;

// Get non-dimension<> root type T of variable<T,O0,O1,O2,...>.
template<typename RealType>
struct root_type_finder { using type = RealType; }; // specialized for dimension<> below.

template<typename RealType,size_t Depth>
struct type_at { using type = RealType; }; // specialized for dimension<> below.

// Satisfies Boost's Conceptual Requirements for Real Number Types.
template<typename RealType,size_t Order>
class dimension
{
    std::array<RealType,Order+1> v;
  public:
    using root_type = typename root_type_finder<RealType>::type; // RealType in the root dimension<RealType,Order>.
    dimension() = default;
    // RealType(cr) | RealType | RealType is copy constructible.
    dimension(const dimension<RealType,Order>&) = default;
    // Be aware of implicit casting from one dimension<> type to another by this copy constructor.
    template<typename RealType2,size_t Order2>
    dimension<RealType,Order>(const dimension<RealType2,Order2>&);
    // RealType(ca) | RealType | RealType is copy constructible from the arithmetic types.
    explicit dimension(const root_type&); // Initialize a variable of differentiation.
    explicit dimension(const std::initializer_list<root_type>&); // Initialize a constant.
    // r = cr | RealType& | Assignment operator.
    dimension<RealType,Order>& operator=(const dimension<RealType,Order>&) = default;
    // r = ca | RealType& | Assignment operator from the arithmetic types.
    dimension<RealType,Order>& operator=(const root_type&); // Set a constant.
    // r += cr | RealType& | Adds cr to r.
    template<typename RealType2,size_t Order2>
    dimension<RealType,Order>& operator+=(const dimension<RealType2,Order2>&);
    // r += ca | RealType& | Adds ar to r.
    dimension<RealType,Order>& operator+=(const root_type&);
    // r -= cr | RealType& | Subtracts cr from r.
    template<typename RealType2,size_t Order2>
    dimension<RealType,Order>& operator-=(const dimension<RealType2,Order2>&);
    // r -= ca | RealType& | Subtracts ca from r.
    dimension<RealType,Order>& operator-=(const root_type&);
    // r *= cr | RealType& | Multiplies r by cr.
    template<typename RealType2,size_t Order2>
    dimension<RealType,Order>& operator*=(const dimension<RealType2,Order2>&);
    // r *= ca | RealType& | Multiplies r by ca.
    dimension<RealType,Order>& operator*=(const root_type&);
    // r /= cr | RealType& | Divides r by cr.
    template<typename RealType2,size_t Order2>
    dimension<RealType,Order>& operator/=(const dimension<RealType2,Order2>&);
    // r /= ca | RealType& | Divides r by ca.
    dimension<RealType,Order>& operator/=(const root_type&);
    // -r | RealType | Unary Negation.
    dimension<RealType,Order> operator-() const;
    // +r | RealType& | Identity Operation.
    const dimension<RealType,Order>& operator+() const;
    // cr + cr2 | RealType | Binary Addition
    template<typename RealType2,size_t Order2>
    promote<dimension<RealType,Order>,dimension<RealType2,Order2>> operator+(const dimension<RealType2,Order2>&) const;
    // cr + ca | RealType | Binary Addition
    dimension<RealType,Order> operator+(const root_type&) const;
    // ca + cr | RealType | Binary Addition
    template<typename RealType2,size_t Order2>
    friend dimension<RealType2,Order2>
        operator+(const typename dimension<RealType2,Order2>::root_type&, const dimension<RealType2,Order2>&);
    // cr - cr2 | RealType | Binary Subtraction
    template<typename RealType2,size_t Order2>
    promote<dimension<RealType,Order>,dimension<RealType2,Order2>> operator-(const dimension<RealType2,Order2>&) const;
    // cr - ca | RealType | Binary Subtraction
    dimension<RealType,Order> operator-(const root_type&) const;
    // ca - cr | RealType | Binary Subtraction
    template<typename RealType2,size_t Order2>
    friend dimension<RealType2,Order2>
        operator-(const typename dimension<RealType2,Order2>::root_type&, const dimension<RealType2,Order2>&);
    // cr * cr2 | RealType | Binary Multiplication
    template<typename RealType2,size_t Order2>
    promote<dimension<RealType,Order>,dimension<RealType2,Order2>> operator*(const dimension<RealType2,Order2>&) const;
    // cr * ca | RealType | Binary Multiplication
    dimension<RealType,Order> operator*(const root_type&) const;
    // ca * cr | RealType | Binary Multiplication
    template<typename RealType2,size_t Order2>
    friend dimension<RealType2,Order2>
        operator*(const typename dimension<RealType2,Order2>::root_type&, const dimension<RealType2,Order2>&);
    // cr / cr2 | RealType | Binary Subtraction
    template<typename RealType2,size_t Order2>
    promote<dimension<RealType,Order>,dimension<RealType2,Order2>> operator/(const dimension<RealType2,Order2>&) const;
    // cr / ca | RealType | Binary Subtraction
    dimension<RealType,Order> operator/(const root_type&) const;
    // ca / cr | RealType | Binary Subtraction
    template<typename RealType2,size_t Order2>
    friend dimension<RealType2,Order2>
        operator/(const typename dimension<RealType2,Order2>::root_type&, const dimension<RealType2,Order2>&);
    // cr == cr2 | bool | Equality Comparison
    template<typename RealType2,size_t Order2> // This only compares the root term. All other terms are ignored.
    bool operator==(const dimension<RealType2,Order2>&) const;
    // cr == ca | bool | Equality Comparison
    bool operator==(const root_type&) const;
    // ca == cr | bool | Equality Comparison
    template<typename RealType2,size_t Order2> // This only compares the root term. All other terms are ignored.
    friend bool operator==(const typename dimension<RealType2,Order2>::root_type&, const dimension<RealType2,Order2>&);
    // cr != cr2 | bool | Inequality Comparison
    template<typename RealType2,size_t Order2>
    bool operator!=(const dimension<RealType2,Order2>&) const;
    // cr != ca | bool | Inequality Comparison
    bool operator!=(const root_type&) const;
    // ca != cr | bool | Inequality Comparison
    template<typename RealType2,size_t Order2>
    friend bool operator!=(const typename dimension<RealType2,Order2>::root_type&, const dimension<RealType2,Order2>&);
    // cr <= cr2 | bool | Less than equal to.
    template<typename RealType2,size_t Order2>
    bool operator<=(const dimension<RealType2,Order2>&) const;
    // cr <= ca | bool | Less than equal to.
    bool operator<=(const root_type&) const;
    // ca <= cr | bool | Less than equal to.
    template<typename RealType2,size_t Order2>
    friend bool operator<=(const typename dimension<RealType2,Order2>::root_type&, const dimension<RealType2,Order2>&);
    // cr >= cr2 | bool | Greater than equal to.
    template<typename RealType2,size_t Order2>
    bool operator>=(const dimension<RealType2,Order2>&) const;
    // cr >= ca | bool | Greater than equal to.
    bool operator>=(const root_type&) const;
    // ca >= cr | bool | Greater than equal to.
    template<typename RealType2,size_t Order2>
    friend bool operator>=(const typename dimension<RealType2,Order2>::root_type&, const dimension<RealType2,Order2>&);
    // cr < cr2 | bool | Less than comparison.
    template<typename RealType2,size_t Order2>
    bool operator<(const dimension<RealType2,Order2>&) const;
    // cr < ca | bool | Less than comparison.
    bool operator<(const root_type&) const;
    // ca < cr | bool | Less than comparison.
    template<typename RealType2,size_t Order2>
    friend bool operator<(const typename dimension<RealType2,Order2>::root_type&, const dimension<RealType2,Order2>&);
    // cr > cr2 | bool | Greater than comparison.
    template<typename RealType2,size_t Order2>
    bool operator>(const dimension<RealType2,Order2>&) const;
    // cr > ca | bool | Greater than comparison.
    bool operator>(const root_type&) const;
    // ca > cr | bool | Greater than comparison.
    template<typename RealType2,size_t Order2>
    friend bool operator>(const typename dimension<RealType2,Order2>::root_type&, const dimension<RealType2,Order2>&);

    // Will throw std::out_of_range if Order < order.
    template<typename... Orders>
    typename type_at<RealType,sizeof...(Orders)>::type at(size_t order, Orders... orders) const;
    template<typename... Orders>
    typename type_at<RealType,sizeof...(Orders)-1>::type derivative(Orders... orders) const;

    dimension<RealType,Order> inverse() const; // Multiplicative inverse.

    static constexpr size_t depth(); // Number of nested std::array<RealType,Order>.
    static constexpr size_t order_sum();

    explicit operator root_type() const; // Must be explicit, otherwise overloaded operators are ambiguous.
    dimension<RealType,Order>& set_root(const root_type&);

    // Use when function returns derivatives.
    dimension<RealType,Order> apply(const std::function<root_type(size_t)>&) const;
    // Use when function returns derivative(i)/factorial(i) (slightly more efficient than apply().)
    dimension<RealType,Order> apply_with_factorials(const std::function<root_type(size_t)>&) const;
    // Same as apply() but uses horner method. May be more accurate in some cases but not as good with inf derivatives.
    dimension<RealType,Order> apply_with_horner(const std::function<root_type(size_t)>&) const;
    // Same as apply_with_factorials() but uses horner method.
    dimension<RealType,Order> apply_with_horner_factorials(const std::function<root_type(size_t)>&) const;

private:
    RealType epsilon_inner_product(size_t z0, size_t isum0, size_t m0,
        const dimension<RealType,Order>& cr, size_t z1, size_t isum1, size_t m1, size_t j) const;
    dimension<RealType,Order> epsilon_multiply(size_t z0, size_t isum0,
        const dimension<RealType,Order>& cr, size_t z1, size_t isum1) const;
    dimension<RealType,Order> epsilon_multiply(size_t z0, size_t isum0, const root_type& ca) const;
    dimension<RealType,Order> inverse_apply() const;
    dimension<RealType,Order>& multiply_assign_by_root_type(bool is_root, const root_type&);

    template<typename RealType2,size_t Orders2>
    friend class dimension;
    template<typename RealType2,size_t Order2>
    friend std::ostream& operator<<(std::ostream&, const dimension<RealType2,Order2>&);

// C++11 Compatibility
#ifdef BOOST_NO_CXX17_IF_CONSTEXPR
    template<typename... Orders>
    typename type_at<RealType, sizeof...(Orders)>::type at_cpp11(std::true_type, size_t order, Orders... orders) const;
    template<typename... Orders>
    typename type_at<RealType, sizeof...(Orders)>::type at_cpp11(std::false_type, size_t order, Orders... orders) const;
    template<typename SizeType>
    dimension<RealType,Order> epsilon_multiply_cpp11(std::true_type, SizeType z0, size_t isum0,
        const dimension<RealType,Order>& cr, size_t z1, size_t isum1) const;
    template<typename SizeType>
    dimension<RealType,Order> epsilon_multiply_cpp11(std::false_type, SizeType z0, size_t isum0,
        const dimension<RealType,Order>& cr, size_t z1, size_t isum1) const;
    template<typename SizeType>
    dimension<RealType,Order> epsilon_multiply_cpp11(std::true_type, SizeType z0, size_t isum0,
        const root_type& ca) const;
    template<typename SizeType>
    dimension<RealType,Order> epsilon_multiply_cpp11(std::false_type, SizeType z0, size_t isum0,
        const root_type& ca) const;
    template<typename RootType>
    dimension<RealType,Order>& multiply_assign_by_root_type_cpp11(std::true_type, bool is_root, const RootType& ca);
    template<typename RootType>
    dimension<RealType,Order>& multiply_assign_by_root_type_cpp11(std::false_type, bool is_root, const RootType& ca);
    template<typename RootType>
    dimension<RealType,Order>& set_root_cpp11(std::true_type, const RootType& root);
    template<typename RootType>
    dimension<RealType,Order>& set_root_cpp11(std::false_type, const RootType& root);
#endif
};

// Standard Library Support Requirements

// fabs(cr1) | RealType
template<typename RealType,size_t Order>
dimension<RealType,Order> fabs(const dimension<RealType,Order>&);
// abs(cr1) | RealType
template<typename RealType,size_t Order>
dimension<RealType,Order> abs(const dimension<RealType,Order>&);
// ceil(cr1) | RealType
template<typename RealType,size_t Order>
dimension<RealType,Order> ceil(const dimension<RealType,Order>&);
// floor(cr1) | RealType
template<typename RealType,size_t Order>
dimension<RealType,Order> floor(const dimension<RealType,Order>&);
// exp(cr1) | RealType
template<typename RealType,size_t Order>
dimension<RealType,Order> exp(const dimension<RealType,Order>&);
// pow(cr, ca) | RealType
template<typename RealType,size_t Order>
dimension<RealType,Order> pow(const dimension<RealType,Order>&,const typename dimension<RealType,Order>::root_type&);
// pow(ca, cr) | RealType
template<typename RealType,size_t Order>
dimension<RealType,Order> pow(const typename dimension<RealType,Order>::root_type&,const dimension<RealType,Order>&);
// pow(cr1, cr2) | RealType
template<typename RealType1,size_t Order1,typename RealType2,size_t Order2>
promote<dimension<RealType1,Order1>,dimension<RealType2,Order2>>
    pow(const dimension<RealType1,Order1>&, const dimension<RealType2,Order2>&);
// sqrt(cr1) | RealType
template<typename RealType,size_t Order>
dimension<RealType,Order> sqrt(const dimension<RealType,Order>&);
// log(cr1) | RealType
template<typename RealType,size_t Order>
dimension<RealType,Order> log(const dimension<RealType,Order>&);
// frexp(cr1, &i) | RealType
template<typename RealType,size_t Order>
dimension<RealType,Order> frexp(const dimension<RealType,Order>&, int*);
// ldexp(cr1, i) | RealType
template<typename RealType,size_t Order>
dimension<RealType,Order> ldexp(const dimension<RealType,Order>&, int);
// cos(cr1) | RealType
template<typename RealType,size_t Order>
dimension<RealType,Order> cos(const dimension<RealType,Order>&);
// sin(cr1) | RealType
template<typename RealType,size_t Order>
dimension<RealType,Order> sin(const dimension<RealType,Order>&);
// asin(cr1) | RealType
template<typename RealType,size_t Order>
dimension<RealType,Order> asin(const dimension<RealType,Order>&);
// tan(cr1) | RealType
template<typename RealType,size_t Order>
dimension<RealType,Order> tan(const dimension<RealType,Order>&);
// atan(cr1) | RealType
template<typename RealType,size_t Order>
dimension<RealType,Order> atan(const dimension<RealType,Order>&);
// fmod(cr1) | RealType
template<typename RealType,size_t Order>
dimension<RealType,Order> fmod(const dimension<RealType,Order>&, const typename dimension<RealType,Order>::root_type&);
// round(cr1) | RealType
template<typename RealType,size_t Order>
dimension<RealType,Order> round(const dimension<RealType,Order>&);
// iround(cr1) | int
template<typename RealType,size_t Order>
int iround(const dimension<RealType,Order>&);
// trunc(cr1) | RealType
template<typename RealType,size_t Order>
dimension<RealType,Order> trunc(const dimension<RealType,Order>&);
// itrunc(cr1) | int
template<typename RealType,size_t Order>
int itrunc(const dimension<RealType,Order>&);

// Additional functions
template<typename RealType,size_t Order>
dimension<RealType,Order> acos(const dimension<RealType,Order>&);
template<typename RealType,size_t Order>
dimension<RealType,Order> erfc(const dimension<RealType,Order>&);
template<typename RealType,size_t Order>
long lround(const dimension<RealType,Order>&);
template<typename RealType,size_t Order>
long long llround(const dimension<RealType,Order>&);
template<typename RealType,size_t Order>
long double truncl(const dimension<RealType,Order>&);

template<typename RealType,size_t Order>
struct nested_dimensions<RealType,Order> { using type = dimension<RealType,Order>; };

template<typename RealType0,typename RealType1>
struct promote_args_n<RealType0,RealType1>
{
    using type = typename boost::math::tools::promote_args_2<RealType0,RealType1>::type;
};

template<typename RealType,size_t Order>
struct root_type_finder<dimension<RealType,Order>> { using type = typename root_type_finder<RealType>::type; };

// Specialization of type_at<> for dimension<>. Examples:
// * type_at<T,0>::type is T.
// * type_at<dimension<T,O1>,1>::type is T.
// * type_at<dimension<dimension<T,O2>,O1>,2>::type is T.
// * type_at<dimension<dimension<dimension<T,O3>,O2>,O1>,3>::type is T.
template<typename RealType,size_t Order,size_t Depth>
struct type_at<dimension<RealType,Order>,Depth>
{
    using type =
        typename std::conditional<Depth==0,dimension<RealType,Order>,typename type_at<RealType,Depth-1>::type>::type;
};

// Compile-time test for dimension<> type.
template<typename>
struct is_dimension : std::false_type {};
template<typename RealType,size_t Order>
struct is_dimension<dimension<RealType,Order>> : std::true_type {};

// C++11 compatibility
#ifdef BOOST_NO_CXX17_IF_CONSTEXPR
#  define BOOST_AUTODIFF_IF_CONSTEXPR
#else
#  define BOOST_AUTODIFF_IF_CONSTEXPR constexpr
#endif

template<typename RealType,size_t Order>
template<typename RealType2,size_t Order2>
dimension<RealType,Order>::dimension(const dimension<RealType2,Order2>& cr)
{
    if BOOST_AUTODIFF_IF_CONSTEXPR (is_dimension<RealType2>::value)
        for (size_t i=0 ; i<=std::min(Order,Order2) ; ++i)
            v[i] = RealType(cr.v[i]);
    else
        for (size_t i=0 ; i<=std::min(Order,Order2) ; ++i)
            v[i] = cr.v[i];
    if BOOST_AUTODIFF_IF_CONSTEXPR (Order2 < Order)
        std::fill(v.begin()+(Order2+1), v.end(), RealType{0});
}

// Note difference between arithmetic constructors and arithmetic assignment:
//  * non-initializer_list arithmetic constructor creates a variable dimension (epsilon coefficient = 1).
//  * initializer_list arithmetic constructor creates a constant dimension (epsilon coefficient = 0).
//  * arithmetic assignment creates a constant (epsilon coefficients = 0).
template<typename RealType,size_t Order>
dimension<RealType,Order>::dimension(const root_type& ca)
:    v{{static_cast<RealType>(ca)}}
{
    if BOOST_AUTODIFF_IF_CONSTEXPR (depth() == 1 && 0 < Order)
        v[1] = static_cast<root_type>(1); // Set epsilon coefficient = 1.
}

template<typename RealType,size_t Order>
dimension<RealType,Order>::dimension(const std::initializer_list<root_type>& list)
:    v{}
{
    for (size_t i=0 ; i<std::min(Order+1,list.size()) ; ++i)
        v[i] = *(list.begin()+i);
}

template<typename RealType,size_t Order>
dimension<RealType,Order>& dimension<RealType,Order>::operator=(const root_type& ca)
{
    v.front() = RealType{ca};
    if BOOST_AUTODIFF_IF_CONSTEXPR (0 < Order)
        std::fill(v.begin()+1, v.end(), RealType{0});
    return *this;
}

template<typename RealType,size_t Order>
template<typename RealType2,size_t Order2>
dimension<RealType,Order>& dimension<RealType,Order>::operator+=(const dimension<RealType2,Order2>& cr)
{
    for (size_t i=0 ; i<=std::min(Order,Order2) ; ++i)
        v[i] += cr.v[i];
    return *this;
}

template<typename RealType,size_t Order>
dimension<RealType,Order>& dimension<RealType,Order>::operator+=(const root_type& ca)
{
    v.front() += ca;
    return *this;
}

template<typename RealType,size_t Order>
template<typename RealType2,size_t Order2>
dimension<RealType,Order>& dimension<RealType,Order>::operator-=(const dimension<RealType2,Order2>& cr)
{
    for (size_t i=0 ; i<=Order ; ++i)
        v[i] -= cr.v[i];
    return *this;
}

template<typename RealType,size_t Order>
dimension<RealType,Order>& dimension<RealType,Order>::operator-=(const root_type& ca)
{
    v.front() -= ca;
    return *this;
}

template<typename RealType,size_t Order>
template<typename RealType2,size_t Order2>
dimension<RealType,Order>& dimension<RealType,Order>::operator*=(const dimension<RealType2,Order2>& cr)
{
    const promote<RealType,RealType2> zero{0};
    if BOOST_AUTODIFF_IF_CONSTEXPR (Order <= Order2)
        for (size_t i=0, j=Order ; i<=Order ; ++i, --j)
            v[j] = std::inner_product(v.cbegin(), v.cend()-i, cr.v.crbegin()+i, zero);
    else
    {
        for (size_t i=0, j=Order ; i<=Order-Order2 ; ++i, --j)
            v[j] = std::inner_product(cr.v.cbegin(), cr.v.cend(), v.crbegin()+i, zero);
        for (size_t i=Order-Order2+1, j=Order2-1 ; i<=Order ; ++i, --j)
            v[j] = std::inner_product(cr.v.cbegin(), cr.v.cbegin()+(j+1), v.crbegin()+i, zero);
    }
    return *this;
}

template<typename RealType,size_t Order>
dimension<RealType,Order>& dimension<RealType,Order>::operator*=(const root_type& ca)
{
    return multiply_assign_by_root_type(true, ca);
}

template<typename RealType,size_t Order>
template<typename RealType2,size_t Order2>
dimension<RealType,Order>& dimension<RealType,Order>::operator/=(const dimension<RealType2,Order2>& cr)
{
    const RealType zero{0};
    v.front() /= cr.v.front();
    if BOOST_AUTODIFF_IF_CONSTEXPR (Order < Order2)
        for (size_t i=1, j=Order2-1, k=Order ; i<=Order ; ++i, --j, --k)
            (v[i] -= std::inner_product(cr.v.cbegin()+1, cr.v.cend()-j, v.crbegin()+k, zero)) /= cr.v.front();
    else if BOOST_AUTODIFF_IF_CONSTEXPR (0 < Order2)
        for (size_t i=1, j=Order2-1, k=Order ; i<=Order ; ++i, j&&--j, --k)
            (v[i] -= std::inner_product(cr.v.cbegin()+1, cr.v.cend()-j, v.crbegin()+k, zero)) /= cr.v.front();
    else
        for (size_t i=1 ; i<=Order ; ++i)
            v[i] /= cr.v.front();
    return *this;
}

template<typename RealType,size_t Order>
dimension<RealType,Order>& dimension<RealType,Order>::operator/=(const root_type& ca)
{
    std::for_each(v.begin(), v.end(), [&ca](RealType& x) { x /= ca; });
    return *this;
}

template<typename RealType,size_t Order>
dimension<RealType,Order> dimension<RealType,Order>::operator-() const
{
    dimension<RealType,Order> retval;
    for (size_t i=0 ; i<=Order ; ++i)
        retval.v[i] = -v[i];
    return retval;
}

template<typename RealType,size_t Order>
const dimension<RealType,Order>& dimension<RealType,Order>::operator+() const
{
    return *this;
}

template<typename RealType,size_t Order>
template<typename RealType2,size_t Order2>
promote<dimension<RealType,Order>,dimension<RealType2,Order2>>
    dimension<RealType,Order>::operator+(const dimension<RealType2,Order2>& cr) const
{
    promote<dimension<RealType,Order>,dimension<RealType2,Order2>> retval;
    for (size_t i=0 ; i<=std::min(Order,Order2) ; ++i)
        retval.v[i] = v[i] + cr.v[i];
    if BOOST_AUTODIFF_IF_CONSTEXPR (Order < Order2)
        for (size_t i=Order+1 ; i<=Order2 ; ++i)
            retval.v[i] = cr.v[i];
    else if BOOST_AUTODIFF_IF_CONSTEXPR (Order2 < Order)
        for (size_t i=Order2+1 ; i<=Order ; ++i)
            retval.v[i] = v[i];
    return retval;
}

template<typename RealType,size_t Order>
dimension<RealType,Order> dimension<RealType,Order>::operator+(const root_type& ca) const
{
    dimension<RealType,Order> retval(*this);
    retval.v.front() += ca;
    return retval;
}

template<typename RealType,size_t Order>
dimension<RealType,Order>
    operator+(const typename dimension<RealType,Order>::root_type& ca, const dimension<RealType,Order>& cr)
{
    return cr + ca;
}

template<typename RealType,size_t Order>
template<typename RealType2,size_t Order2>
promote<dimension<RealType,Order>,dimension<RealType2,Order2>>
    dimension<RealType,Order>::operator-(const dimension<RealType2,Order2>& cr) const
{
    promote<dimension<RealType,Order>,dimension<RealType2,Order2>> retval;
    for (size_t i=0 ; i<=std::min(Order,Order2) ; ++i)
        retval.v[i] = v[i] - cr.v[i];
    if BOOST_AUTODIFF_IF_CONSTEXPR (Order < Order2)
        for (size_t i=Order+1 ; i<=Order2 ; ++i)
            retval.v[i] = -cr.v[i];
    else if BOOST_AUTODIFF_IF_CONSTEXPR (Order2 < Order)
        for (size_t i=Order2+1 ; i<=Order ; ++i)
            retval.v[i] = v[i];
    return retval;
}

template<typename RealType,size_t Order>
dimension<RealType,Order> dimension<RealType,Order>::operator-(const root_type& ca) const
{
    dimension<RealType,Order> retval(*this);
    retval.v.front() -= ca;
    return retval;
}

template<typename RealType,size_t Order>
dimension<RealType,Order>
    operator-(const typename dimension<RealType,Order>::root_type& ca, const dimension<RealType,Order>& cr)
{
    return -cr += ca;
}

template<typename RealType,size_t Order>
template<typename RealType2,size_t Order2>
promote<dimension<RealType,Order>,dimension<RealType2,Order2>>
    dimension<RealType,Order>::operator*(const dimension<RealType2,Order2>& cr) const
{
    const promote<RealType,RealType2> zero{0};
    promote<dimension<RealType,Order>,dimension<RealType2,Order2>> retval;
    if BOOST_AUTODIFF_IF_CONSTEXPR (Order < Order2)
        for (size_t i=0, j=Order, k=Order2 ; i<=Order2 ; ++i, j&&--j, --k)
            retval.v[i] = std::inner_product(v.cbegin(), v.cend()-j, cr.v.crbegin()+k, zero);
    else
        for (size_t i=0, j=Order2, k=Order ; i<=Order ; ++i, j&&--j, --k)
            retval.v[i] = std::inner_product(cr.v.cbegin(), cr.v.cend()-j, v.crbegin()+k, zero);
    return retval;
}

template<typename RealType,size_t Order>
dimension<RealType,Order> dimension<RealType,Order>::operator*(const root_type& ca) const
{
    return dimension<RealType,Order>(*this) *= ca;
}

template<typename RealType,size_t Order>
dimension<RealType,Order>
    operator*(const typename dimension<RealType,Order>::root_type& ca, const dimension<RealType,Order>& cr)
{
    return cr * ca;
}

template<typename RealType,size_t Order>
template<typename RealType2,size_t Order2>
promote<dimension<RealType,Order>,dimension<RealType2,Order2>>
    dimension<RealType,Order>::operator/(const dimension<RealType2,Order2>& cr) const
{
    const promote<RealType,RealType2> zero{0};
    promote<dimension<RealType,Order>,dimension<RealType2,Order2>> retval;
    retval.v.front() = v.front() / cr.v.front();
    if BOOST_AUTODIFF_IF_CONSTEXPR (Order < Order2)
    {
        for (size_t i=1, j=Order2-1 ; i<=Order ; ++i, --j)
            retval.v[i] = (v[i] -
                std::inner_product(cr.v.cbegin()+1, cr.v.cend()-j, retval.v.crbegin()+(j+1), zero)) / cr.v.front();
        for (size_t i=Order+1, j=Order2-Order-1 ; i<=Order2 ; ++i, --j)
            retval.v[i] =
                -std::inner_product(cr.v.cbegin()+1, cr.v.cend()-j, retval.v.crbegin()+(j+1), zero) / cr.v.front();
    }
    else if BOOST_AUTODIFF_IF_CONSTEXPR (0 < Order2)
        for (size_t i=1, j=Order2-1, k=Order ; i<=Order ; ++i, j&&--j, --k)
            retval.v[i] =
                (v[i] - std::inner_product(cr.v.cbegin()+1, cr.v.cend()-j, retval.v.crbegin()+k, zero)) / cr.v.front();
    else
        for (size_t i=1 ; i<=Order ; ++i)
            retval.v[i] = v[i] / cr.v.front();
    return retval;
}

template<typename RealType,size_t Order>
dimension<RealType,Order> dimension<RealType,Order>::operator/(const root_type& ca) const
{
    return dimension<RealType,Order>(*this) /= ca;
}

template<typename RealType,size_t Order>
dimension<RealType,Order>
    operator/(const typename dimension<RealType,Order>::root_type& ca, const dimension<RealType,Order>& cr)
{
    dimension<RealType,Order> retval;
    retval.v.front() = ca / cr.v.front();
    if BOOST_AUTODIFF_IF_CONSTEXPR (0 < Order)
    {
        const RealType zero{0};
        for (size_t i=1, j=Order-1 ; i<=Order ; ++i, --j)
            retval.v[i] = -std::inner_product(cr.v.cbegin()+1, cr.v.cend()-j, retval.v.crbegin()+(j+1), zero)
                / cr.v.front();
    }
    return retval;
}

template<typename RealType,size_t Order>
template<typename RealType2,size_t Order2>
bool dimension<RealType,Order>::operator==(const dimension<RealType2,Order2>& cr) const
{
    return v.front() == cr.v.front();
}

template<typename RealType,size_t Order>
bool dimension<RealType,Order>::operator==(const root_type& ca) const
{
    return v.front() == ca;
}

template<typename RealType,size_t Order>
bool operator==(const typename dimension<RealType,Order>::root_type& ca, const dimension<RealType,Order>& cr)
{
    return ca == cr.v.front();
}

template<typename RealType,size_t Order>
template<typename RealType2,size_t Order2>
bool dimension<RealType,Order>::operator!=(const dimension<RealType2,Order2>& cr) const
{
    return v.front() != cr.v.front();
}

template<typename RealType,size_t Order>
bool dimension<RealType,Order>::operator!=(const root_type& ca) const
{
    return v.front() != ca;
}

template<typename RealType,size_t Order>
bool operator!=(const typename dimension<RealType,Order>::root_type& ca, const dimension<RealType,Order>& cr)
{
    return ca != cr.v.front();
}

template<typename RealType,size_t Order>
template<typename RealType2,size_t Order2>
bool dimension<RealType,Order>::operator<=(const dimension<RealType2,Order2>& cr) const
{
    return v.front() <= cr.v.front();
}

template<typename RealType,size_t Order>
bool dimension<RealType,Order>::operator<=(const root_type& ca) const
{
    return v.front() <= ca;
}

template<typename RealType,size_t Order>
bool operator<=(const typename dimension<RealType,Order>::root_type& ca, const dimension<RealType,Order>& cr)
{
    return ca <= cr.v.front();
}

template<typename RealType,size_t Order>
template<typename RealType2,size_t Order2>
bool dimension<RealType,Order>::operator>=(const dimension<RealType2,Order2>& cr) const
{
    return v.front() >= cr.v.front();
}

template<typename RealType,size_t Order>
bool dimension<RealType,Order>::operator>=(const root_type& ca) const
{
    return v.front() >= ca;
}

template<typename RealType,size_t Order>
bool operator>=(const typename dimension<RealType,Order>::root_type& ca, const dimension<RealType,Order>& cr)
{
    return ca >= cr.v.front();
}

template<typename RealType,size_t Order>
template<typename RealType2,size_t Order2>
bool dimension<RealType,Order>::operator<(const dimension<RealType2,Order2>& cr) const
{
    return v.front() < cr.v.front();
}

template<typename RealType,size_t Order>
bool dimension<RealType,Order>::operator<(const root_type& ca) const
{
    return v.front() < ca;
}

template<typename RealType,size_t Order>
bool operator<(const typename dimension<RealType,Order>::root_type& ca, const dimension<RealType,Order>& cr)
{
    return ca < cr.v.front();
}

template<typename RealType,size_t Order>
template<typename RealType2,size_t Order2>
bool dimension<RealType,Order>::operator>(const dimension<RealType2,Order2>& cr) const
{
    return v.front() > cr.v.front();
}

template<typename RealType,size_t Order>
bool dimension<RealType,Order>::operator>(const root_type& ca) const
{
    return v.front() > ca;
}

template<typename RealType,size_t Order>
bool operator>(const typename dimension<RealType,Order>::root_type& ca, const dimension<RealType,Order>& cr)
{
    return ca > cr.v.front();
}

/*** Other methods and functions ***/

// f : order -> derivative(order)
template<typename RealType,size_t Order>
dimension<RealType,Order> dimension<RealType,Order>::apply(const std::function<root_type(size_t)>& f) const
{
    const dimension<RealType,Order> epsilon = dimension<RealType,Order>(*this).set_root(0);
    dimension<RealType,Order> epsilon_i = dimension<RealType,Order>{1}; // epsilon to the power of i
    dimension<RealType,Order> accumulator = dimension<RealType,Order>{f(0)};
    for (size_t i=1 ; i<=order_sum() ; ++i)
    {    // accumulator += (epsilon_i *= epsilon) * (f(i) / boost::math::factorial<root_type>(i));
        epsilon_i = epsilon_i.epsilon_multiply(i-1, 0, epsilon, 1, 0);
        accumulator += epsilon_i.epsilon_multiply(i, 0, f(i) / boost::math::factorial<root_type>(i));
    }
    return accumulator;
}

// f : order -> derivative(order)/factorial(order)
// Use this when the computation of the derivatives already includes the factorial terms. E.g. See atan().
template<typename RealType,size_t Order>
dimension<RealType,Order>
    dimension<RealType,Order>::apply_with_factorials(const std::function<root_type(size_t)>& f) const
{
    const dimension<RealType,Order> epsilon = dimension<RealType,Order>(*this).set_root(0);
    dimension<RealType,Order> epsilon_i = dimension<RealType,Order>{1}; // epsilon to the power of i
    dimension<RealType,Order> accumulator = dimension<RealType,Order>{f(0)};
    for (size_t i=1 ; i<=order_sum() ; ++i)
    {    // accumulator += (epsilon_i *= epsilon) * f(i);
        epsilon_i = epsilon_i.epsilon_multiply(i-1, 0, epsilon, 1, 0);
        accumulator += epsilon_i.epsilon_multiply(i, 0, f(i));
    }
    return accumulator;
}

// f : order -> derivative(order)
template<typename RealType,size_t Order>
dimension<RealType,Order> dimension<RealType,Order>::apply_with_horner(const std::function<root_type(size_t)>& f) const
{
    const dimension<RealType,Order> epsilon = dimension<RealType,Order>(*this).set_root(0);
    auto accumulator = dimension<RealType,Order>{f(order_sum())/boost::math::factorial<root_type>(order_sum())};
    for (size_t i=order_sum() ; i-- ;)
        (accumulator *= epsilon) += f(i) / boost::math::factorial<root_type>(i);
    return accumulator;
}

// f : order -> derivative(order)/factorial(order)
// Use this when the computation of the derivatives already includes the factorial terms. E.g. See atan().
template<typename RealType,size_t Order>
dimension<RealType,Order>
    dimension<RealType,Order>::apply_with_horner_factorials(const std::function<root_type(size_t)>& f) const
{
    const dimension<RealType,Order> epsilon = dimension<RealType,Order>(*this).set_root(0);
    auto accumulator = dimension<RealType,Order>{f(order_sum())};
    for (size_t i=order_sum() ; i-- ;)
        (accumulator *= epsilon) += f(i);
    return accumulator;
}

#ifndef BOOST_NO_CXX17_IF_CONSTEXPR
// Can throw "std::out_of_range: array::at: __n (which is 7) >= _Nm (which is 7)"
template<typename RealType,size_t Order>
template<typename... Orders>
typename type_at<RealType,sizeof...(Orders)>::type dimension<RealType,Order>::at(size_t order, Orders... orders) const
{
    if constexpr (0 < sizeof...(orders))
        return v.at(order).at(orders...);
    else
        return v.at(order);
}
#endif

#ifndef BOOST_NO_CXX17_IF_CONSTEXPR
template<typename RealType,size_t Order>
constexpr size_t dimension<RealType,Order>::depth()
{
    if constexpr (is_dimension<RealType>::value)
        return 1 + RealType::depth();
    else
        return 1;
}
#endif

#ifndef BOOST_NO_CXX17_IF_CONSTEXPR
template<typename RealType,size_t Order>
constexpr size_t dimension<RealType,Order>::order_sum()
{
    if constexpr (is_dimension<RealType>::value)
        return Order + RealType::order_sum();
    else
        return Order;
}
#endif

#ifndef BOOST_NO_CXX17_IF_CONSTEXPR
// Can throw "std::out_of_range: array::at: __n (which is 7) >= _Nm (which is 7)"
template<typename RealType,size_t Order>
template<typename... Orders>
typename type_at<RealType,sizeof...(Orders)-1>::type dimension<RealType,Order>::derivative(Orders... orders) const
{
    static_assert(sizeof...(orders) <= depth(),
        "Number of parameters to derivative(...) cannot exceed the number of dimensions in the dimension<...>.");
    return at(orders...) * (... * boost::math::factorial<root_type>(orders));
}
#endif

template<typename RealType,size_t Order>
RealType dimension<RealType,Order>::epsilon_inner_product(size_t z0, size_t isum0, size_t m0,
    const dimension<RealType,Order>& cr, size_t z1, size_t isum1, size_t m1, size_t j) const
{
    static_assert(is_dimension<RealType>::value, "epsilon_inner_product() must have 1 < depth().");
    RealType accumulator = RealType();
    const size_t i0_max = m1 < j ? j-m1 : 0;
    for (size_t i0=m0, i1=j-m0 ; i0<=i0_max ; ++i0, --i1)
        accumulator += v.at(i0).epsilon_multiply(z0, isum0+i0, cr.v.at(i1), z1, isum1+i1);
    return accumulator;
}

#ifndef BOOST_NO_CXX17_IF_CONSTEXPR
template<typename RealType,size_t Order>
dimension<RealType,Order> dimension<RealType,Order>::epsilon_multiply(size_t z0, size_t isum0,
    const dimension<RealType,Order>& cr, size_t z1, size_t isum1) const
{
    const RealType zero{0};
    const size_t m0 = order_sum() + isum0 < Order + z0 ? Order + z0 - (order_sum() + isum0) : 0;
    const size_t m1 = order_sum() + isum1 < Order + z1 ? Order + z1 - (order_sum() + isum1) : 0;
    const size_t i_max = m0 + m1 < Order ? Order - (m0 + m1) : 0;
    dimension<RealType,Order> retval = dimension<RealType,Order>();
    if constexpr (is_dimension<RealType>::value)
        for (size_t i=0, j=Order ; i<=i_max ; ++i, --j)
            retval.v[j] = epsilon_inner_product(z0, isum0, m0, cr, z1, isum1, m1, j);
    else
        for (size_t i=0, j=Order ; i<=i_max ; ++i, --j)
            retval.v[j] = std::inner_product(v.cbegin()+m0, v.cend()-(i+m1), cr.v.crbegin()+(i+m0), zero);
    return retval;
}
#endif

#ifndef BOOST_NO_CXX17_IF_CONSTEXPR
// When called from outside this method, z0 should be non-zero. Otherwise if z0=0 then it will give an
// incorrect result of 0 when the root value is 0 and ca=inf, when instead the correct product is nan.
// If z0=0 then use the regular multiply operator*() instead.
template<typename RealType,size_t Order>
dimension<RealType,Order> dimension<RealType,Order>::epsilon_multiply(size_t z0, size_t isum0,
    const root_type& ca) const
{
    dimension<RealType,Order> retval(*this);
    const size_t m0 = order_sum() + isum0 < Order + z0 ? Order + z0 - (order_sum() + isum0) : 0;
    if constexpr (is_dimension<RealType>::value)
        for (size_t i=m0 ; i<=Order ; ++i)
            retval.v[i] = retval.v[i].epsilon_multiply(z0, isum0+i, ca);
    else
        for (size_t i=m0 ; i<=Order ; ++i)
            if (retval.v[i] != static_cast<RealType>(0))
                retval.v[i] *= ca;
    return retval;
}
#endif

template<typename RealType,size_t Order>
dimension<RealType,Order> dimension<RealType,Order>::inverse() const
{
    return operator root_type() == 0 ? inverse_apply() : 1 / *this;
}

// This gives autodiff::log(0.0) = depth(1)(-inf,inf,-inf,inf,-inf,inf)
// 1 / *this: autodiff::log(0.0) = depth(1)(-inf,inf,-inf,-nan,-nan,-nan)
template<typename RealType,size_t Order>
dimension<RealType,Order> dimension<RealType,Order>::inverse_apply() const
{
    root_type derivatives[order_sum()+1]; // derivatives of 1/x
    const root_type x0 = static_cast<root_type>(*this);
    derivatives[0] = 1 / x0;
    for (size_t i=1 ; i<=order_sum() ; ++i)
        derivatives[i] = -derivatives[i-1] * i / x0;
    return apply([&derivatives](size_t j) { return derivatives[j]; });
}

#ifndef BOOST_NO_CXX17_IF_CONSTEXPR
template<typename RealType,size_t Order>
dimension<RealType,Order>& dimension<RealType,Order>::multiply_assign_by_root_type(bool is_root, const root_type& ca)
{
    auto itr = v.begin();
    if constexpr (is_dimension<RealType>::value)
    {
        itr->multiply_assign_by_root_type(is_root, ca);
        for (++itr ; itr!=v.end() ; ++itr)
            itr->multiply_assign_by_root_type(false, ca);
    }
    else
    {
        if (is_root || *itr != 0)
            *itr *= ca; // Skip multiplication of 0 by ca=inf to avoid nan. Exception: root value is always multiplied.
        for (++itr ; itr!=v.end() ; ++itr)
            if (*itr != 0)
                *itr *= ca;
    }
    return *this;
}
#endif

template<typename RealType,size_t Order>
dimension<RealType,Order>::operator root_type() const
{
    return static_cast<root_type>(v.front());
}

#ifndef BOOST_NO_CXX17_IF_CONSTEXPR
template<typename RealType,size_t Order>
dimension<RealType,Order>& dimension<RealType,Order>::set_root(const root_type& root)
{
    if constexpr (is_dimension<RealType>::value)
        v.front().set_root(root);
    else
        v.front() = root;
    return *this;
}
#endif

// Standard Library Support Requirements

template<typename RealType,size_t Order>
dimension<RealType,Order> fabs(const dimension<RealType,Order>& cr)
{
    const typename dimension<RealType,Order>::root_type zero{0};
    return zero < cr ? cr
        : cr < zero ? -cr
        : cr == zero ? dimension<RealType,Order>() // Canonical fabs'(0) = 0.
        : cr; // Propagate NaN.
}

template<typename RealType,size_t Order>
dimension<RealType,Order> abs(const dimension<RealType,Order>& cr)
{
    return fabs(cr);
}

template<typename RealType,size_t Order>
dimension<RealType,Order> ceil(const dimension<RealType,Order>& cr)
{
    using std::ceil;
    return dimension<RealType,Order>{ceil(static_cast<typename dimension<RealType,Order>::root_type>(cr))};
}

template<typename RealType,size_t Order>
dimension<RealType,Order> floor(const dimension<RealType,Order>& cr)
{
    using std::floor;
    return dimension<RealType,Order>{floor(static_cast<typename dimension<RealType,Order>::root_type>(cr))};
}

template<typename RealType,size_t Order>
dimension<RealType,Order> exp(const dimension<RealType,Order>& cr)
{
    using std::exp;
    using root_type = typename dimension<RealType,Order>::root_type;
    const root_type d0 = exp(static_cast<root_type>(cr));
    return cr.apply_with_horner([&d0](size_t) { return d0; });
}

template<typename RealType,size_t Order>
dimension<RealType,Order> pow(const dimension<RealType,Order>& x,const typename dimension<RealType,Order>::root_type& y)
{
    using std::pow;
    using root_type = typename dimension<RealType,Order>::root_type;
    constexpr size_t order = dimension<RealType,Order>::order_sum();
    std::array<root_type,order+1> derivatives; // array of derivatives
    const root_type x0 = static_cast<root_type>(x);
    size_t i = 0;
    root_type coef = 1;
    for (; i<=order && coef!=0 ; ++i)
    {
        derivatives[i] = coef * pow(x0, y-i);
        coef *= y - i;
    }
    return x.apply([&derivatives,i](size_t j) { return j < i ? derivatives[j] : 0; });
}

template<typename RealType,size_t Order>
dimension<RealType,Order> pow(const typename dimension<RealType,Order>::root_type& x,const dimension<RealType,Order>& y)
{
    using std::log;
    return exp(y*log(x));
}

template<typename RealType1,size_t Order1,typename RealType2,size_t Order2>
promote<dimension<RealType1,Order1>,dimension<RealType2,Order2>>
    pow(const dimension<RealType1,Order1>& x, const dimension<RealType2,Order2>& y)
{
    return exp(y*log(x));
}

template<typename RealType,size_t Order>
dimension<RealType,Order> sqrt(const dimension<RealType,Order>& cr)
{
    return pow(cr,0.5);
}

// Natural logarithm. If cr==0 then derivative(i) may have nans due to nans from inverse().
template<typename RealType,size_t Order>
dimension<RealType,Order> log(const dimension<RealType,Order>& cr)
{
    using std::log;
    using root_type = typename dimension<RealType,Order>::root_type;
    constexpr size_t order = dimension<RealType,Order>::order_sum();
    const root_type d0 = log(static_cast<root_type>(cr));
    if BOOST_AUTODIFF_IF_CONSTEXPR (order == 0)
        return dimension<RealType,0>(d0);
    else
    {
        const auto d1 = dimension<root_type,order-1>(static_cast<root_type>(cr)).inverse(); // log'(x) = 1 / x
        return cr.apply_with_factorials([&d0,&d1](size_t i) { return i ? d1.at(i-1)/i : d0; });
    }
}

template<typename RealType,size_t Order>
dimension<RealType,Order> frexp(const dimension<RealType,Order>& cr, int* exp)
{
    using std::exp2;
    using std::frexp;
    using root_type = typename dimension<RealType,Order>::root_type;
    frexp(static_cast<root_type>(cr), exp);
    return cr * exp2(-*exp);
}

template<typename RealType,size_t Order>
dimension<RealType,Order> ldexp(const dimension<RealType,Order>& cr, int exp)
{
    using std::exp2;
    return cr * exp2(exp);
}

template<typename RealType,size_t Order>
dimension<RealType,Order> cos(const dimension<RealType,Order>& cr)
{
    using std::cos;
    using std::sin;
    using root_type = typename dimension<RealType,Order>::root_type;
    const root_type d0 = cos(static_cast<root_type>(cr));
    if BOOST_AUTODIFF_IF_CONSTEXPR (dimension<RealType,Order>::order_sum() == 0)
        return dimension<RealType,0>(d0);
    else
    {
        const root_type d1 = -sin(static_cast<root_type>(cr));
        const root_type derivatives[] { d0, d1, -d0, -d1 };
        return cr.apply_with_horner([&derivatives](size_t i) { return derivatives[i&3]; });
    }
}

template<typename RealType,size_t Order>
dimension<RealType,Order> sin(const dimension<RealType,Order>& cr)
{
    using std::sin;
    using std::cos;
    using root_type = typename dimension<RealType,Order>::root_type;
    const root_type d0 = sin(static_cast<root_type>(cr));
    if BOOST_AUTODIFF_IF_CONSTEXPR (dimension<RealType,Order>::order_sum() == 0)
        return dimension<RealType,0>(d0);
    else
    {
        const root_type d1 = cos(static_cast<root_type>(cr));
        const root_type derivatives[] { d0, d1, -d0, -d1 };
        return cr.apply_with_horner([&derivatives](size_t i) { return derivatives[i&3]; });
    }
}

template<typename RealType,size_t Order>
dimension<RealType,Order> asin(const dimension<RealType,Order>& cr)
{
    using std::asin;
    using root_type = typename dimension<RealType,Order>::root_type;
    constexpr size_t order = dimension<RealType,Order>::order_sum();
    const root_type d0 = asin(static_cast<root_type>(cr));
    if BOOST_AUTODIFF_IF_CONSTEXPR (order == 0)
        return dimension<RealType,0>(d0);
    else
    {
        auto d1 = dimension<root_type,order-1>(static_cast<root_type>(cr)); // asin'(x) = 1 / sqrt(1-x*x).
        d1 = sqrt(1-(d1*=d1)).inverse(); // asin(1): d1 = depth(1)(inf,inf,-nan,-nan,-nan)
        //d1 = sqrt((1-(d1*=d1)).inverse()); // asin(1): d1 = depth(1)(inf,-nan,-nan,-nan,-nan)
        return cr.apply_with_factorials([&d0,&d1](size_t i) { return i ? d1.at(i-1)/i : d0; });
    }
}

template<typename RealType,size_t Order>
dimension<RealType,Order> tan(const dimension<RealType,Order>& cr)
{
    return sin(cr) / cos(cr);
}

template<typename RealType,size_t Order>
dimension<RealType,Order> atan(const dimension<RealType,Order>& cr)
{
    using std::atan;
    using root_type = typename dimension<RealType,Order>::root_type;
    constexpr size_t order = dimension<RealType,Order>::order_sum();
    const root_type d0 = atan(static_cast<root_type>(cr));
    if BOOST_AUTODIFF_IF_CONSTEXPR (order == 0)
        return dimension<RealType,0>(d0);
    else
    {
        auto d1 = dimension<root_type,order-1>(static_cast<root_type>(cr));
        d1 = ((d1*=d1)+=1).inverse(); // atan'(x) = 1 / (x*x+1).
        return cr.apply_with_horner_factorials([&d0,&d1](size_t i) { return i ? d1.at(i-1)/i : d0; });
    }
}

template<typename RealType,size_t Order>
dimension<RealType,Order>
    fmod(const dimension<RealType,Order>& cr, const typename dimension<RealType,Order>::root_type& ca)
{
    using std::fmod;
    using root_type = typename dimension<RealType,Order>::root_type;
    return dimension<RealType,Order>(cr).set_root(0) += fmod(static_cast<root_type>(cr), ca);
}

template<typename RealType,size_t Order>
dimension<RealType,Order> round(const dimension<RealType,Order>& cr)
{
    using std::round;
    return dimension<RealType,Order>{round(static_cast<typename dimension<RealType,Order>::root_type>(cr))};
}

template<typename RealType,size_t Order>
int iround(const dimension<RealType,Order>& cr)
{
    using boost::math::iround;
    return iround(static_cast<typename dimension<RealType,Order>::root_type>(cr));
}

template<typename RealType,size_t Order>
dimension<RealType,Order> trunc(const dimension<RealType,Order>& cr)
{
    using std::trunc;
    return dimension<RealType,Order>{trunc(static_cast<typename dimension<RealType,Order>::root_type>(cr))};
}

template<typename RealType,size_t Order>
int itrunc(const dimension<RealType,Order>& cr)
{
    using boost::math::itrunc;
    return itrunc(static_cast<typename dimension<RealType,Order>::root_type>(cr));
}

template<typename RealType,size_t Order>
std::ostream& operator<<(std::ostream& out, const dimension<RealType,Order>& dim)
{
    const std::streamsize original_precision = out.precision();
    out.precision(std::numeric_limits<typename dimension<RealType,Order>::root_type>::digits10);
    out << "depth(" << dim.depth() << ')';
    for (size_t i=0 ; i<dim.v.size() ; ++i)
        out << (i?',':'(') << dim.v[i];
    out.precision(original_precision);
    return out << ')';
}

// Additional functions

template<typename RealType,size_t Order>
dimension<RealType,Order> acos(const dimension<RealType,Order>& cr)
{
    using std::acos;
    using root_type = typename dimension<RealType,Order>::root_type;
    constexpr size_t order = dimension<RealType,Order>::order_sum();
    const root_type d0 = acos(static_cast<root_type>(cr));
    if BOOST_AUTODIFF_IF_CONSTEXPR (order == 0)
        return dimension<RealType,0>(d0);
    else
    {
        auto d1 = dimension<root_type,order-1>(static_cast<root_type>(cr));
        d1 = -sqrt(1-(d1*=d1)).inverse(); // acos'(x) = -1 / sqrt(1-x*x).
        return cr.apply_with_horner_factorials([&d0,&d1](size_t i) { return i ? d1.at(i-1)/i : d0; });
    }
}

template<typename RealType,size_t Order>
dimension<RealType,Order> erfc(const dimension<RealType,Order>& cr)
{
    using std::erfc;
    using root_type = typename dimension<RealType,Order>::root_type;
    constexpr size_t order = dimension<RealType,Order>::order_sum();
    const root_type d0 = erfc(static_cast<root_type>(cr));
    if BOOST_AUTODIFF_IF_CONSTEXPR (order == 0)
        return dimension<RealType,0>(d0);
    else
    {
        auto d1 = dimension<root_type,order-1>(static_cast<root_type>(cr));
        d1 = -2*boost::math::constants::one_div_root_pi<root_type>()*exp(-(d1*=d1)); // erfc'(x)=-2/sqrt(pi)*exp(-x*x)
        return cr.apply_with_horner_factorials([&d0,&d1](size_t i) { return i ? d1.at(i-1)/i : d0; });
    }
}

template<typename RealType,size_t Order>
long lround(const dimension<RealType,Order>& cr)
{
    using std::lround;
    return lround(static_cast<typename dimension<RealType,Order>::root_type>(cr));
}

template<typename RealType,size_t Order>
long long llround(const dimension<RealType,Order>& cr)
{
    using std::llround;
    return llround(static_cast<typename dimension<RealType,Order>::root_type>(cr));
}

template<typename RealType,size_t Order>
long double truncl(const dimension<RealType,Order>& cr)
{
    using std::truncl;
    return truncl(static_cast<typename dimension<RealType,Order>::root_type>(cr));
}

} } } } } // namespace boost::math::differentiation::autodiff::v1

namespace std {
/// boost::math::tools::digits<RealType>() is handled by this std::numeric_limits<> specialization,
/// and similarly for max_value, min_value, log_max_value, log_min_value, and epsilon.
template <typename RealType,size_t Order>
class numeric_limits<boost::math::differentiation::autodiff::dimension<RealType,Order>>
    : public numeric_limits<typename boost::math::differentiation::autodiff::dimension<RealType,Order>::root_type>
{ };
} // namespace std

namespace boost { namespace math { namespace tools {

// See boost/math/tools/promotion.hpp
template <typename RealType0,size_t Order0,typename RealType1,size_t Order1>
struct promote_args_2<differentiation::autodiff::dimension<RealType0,Order0>,differentiation::autodiff::dimension<RealType1,Order1>>
{
    using type = differentiation::autodiff::dimension<typename promote_args_2<RealType0,RealType1>::type,
#ifndef BOOST_NO_CXX14_CONSTEXPR
        std::max(Order0,Order1)>;
#else
        Order0 < Order1 ? Order1 : Order0>;
#endif
};

template <typename RealType0,size_t Order0,typename RealType1>
struct promote_args_2<differentiation::autodiff::dimension<RealType0,Order0>,RealType1>
{
    using type = differentiation::autodiff::dimension<typename promote_args_2<RealType0,RealType1>::type,Order0>;
};

template <typename RealType0,typename RealType1,size_t Order1>
struct promote_args_2<RealType0,differentiation::autodiff::dimension<RealType1,Order1>>
{
    using type = differentiation::autodiff::dimension<typename promote_args_2<RealType0,RealType1>::type,Order1>;
};

} } } // namespace boost::math::tools

#ifdef BOOST_NO_CXX17_IF_CONSTEXPR
#include "autodiff_cpp11.hpp"
#endif

#endif // BOOST_MATH_AUTODIFF_HPP
