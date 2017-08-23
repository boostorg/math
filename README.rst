boost math
==========

boost math is a suite of tools useful for scientific computing and numerical analysis.
The tools provided include evaluation of special functions, statistical distributions,
regular and irregular interpolation, numerical integration, polynomial evaluation, and root-finding,
all at arbitrary precision.
Documentation for boost math can be found at boost.org_, as well as documentation for other libraries.


--------------------------------------
Desired features for upcoming releases
--------------------------------------

Despite being ~240,000 lines of source code, there is considerable scope for expansion of the library.
Current features which are missing include

- Monte-Carlo integration, Frank-Wolfe Bayesian quadrature
- Fast Fourier transforms (include discrete sine and cosine transforms, Hilbert transform)
- Integer factorization, including trial division, and the Pollard rho algorithm, sieves
- Continuous wavelet evaluation, discrete wavelet transforms
- Chebyshev polynomial interpolation
- Multidimensional interpolation
- Sparse grid quadrature

Assistance with these features would be greatly appreciated.

------------
Contributing
------------

boost math is a complex project, and contributing is also complex.
The following instructions will hopefully help the user get started.

- Navigate to https://github.com/boostorg/math, and click the "Fork" button.

- Clone boost:

    .. code:: bash

        $ git clone --recursive https://github.com/boostorg/boost.git

- Remove the math submodule from the boost superproject, and clone your fork into its location.

    .. code:: bash

        boost$ cd libs
        boost/libs$ rm -rf math; git clone https://github.com/MY-USERNAME/math.git
        boost/libs$ cd math;

- To keep your fork synced to the master, add an upstream:

    .. code:: bash

        boost/libs/math$ git remote add upstream https://github.com/boostorg/math.git

- Ensure boost builds successfully:

    .. code:: bash

        boost$ ./bootstrap.sh
        boost$ ./b2

- To build the documentation, you must install some dependencies:
    .. code:: bash

        boost/libs/math/doc$ sudo apt install xsltproc docbook-xsl docbook-xml
        boost/libs/math/doc$ ../../../b2 release
This will build the `boost::math` documentation into `doc/index.html` and hence can be viewed with any browser. In addition, `doc/math.pdf` should contain the same information.

- Write your feature. For the purposes of this contributing guide, we will test if a floating point number is an integer. In this case, place the file `is_float_int.hpp` in (say) the `include/boost/math/tools` directory.

    .. code:: cpp

        /*
         *  Copyright John Maddock, Nick Thompson, 2017
         *  Use, modification and distribution are subject to the
         *  Boost Software License, Version 1.0. (See accompanying file
         *  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
         */
        #ifndef BOOST_MATH_TOOLS_IS_FLOAT_INT_HPP
        #define BOOST_MATH_TOOLS_IS_FLOAT_INT_HPP

        namespace boost { namespace math {

        template<class Real>
        bool is_float_integer(Real x) {
            using std::floor;
            if (floor(x) == x) {
                return true;
            }
            return false;
        }
        }}
        #endif

A few things to note: All contributions to boost are released under the boost license, and as such you must include a copyright and license statement at the beginning of the file.
Next are the standard include guards, and namespaces.
Note that we write `using std::floor` and `floor(x) == x` instead of `std::floor(x) == x`.
This is because the latter form does not work with boost multiprecision types.
The former uses argument-dependent lookup to determine which version of `floor` (the one from `cmath` or the version from `boost::multiprecision`) should be used.

- Write a unit test. Place `is_float_int_test.cpp` in the `test` directory. This test should verify that your code works for many cases and many types.

    .. code:: cpp
        /*
         * Copyright Nick Thompson, 2017
         * Use, modification and distribution are subject to the
         * Boost Software License, Version 1.0. (See accompanying file
         * LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
         */
        #define BOOST_TEST_MODULE test_is_float_int

        #include <boost/test/included/unit_test.hpp>
        #include <boost/test/floating_point_comparison.hpp>
        #include <boost/math/concepts/real_concept.hpp>
        #include <boost/math/tools/is_float_int.hpp>
        #include <boost/multiprecision/cpp_bin_float.hpp>
        #include <boost/multiprecision/cpp_dec_float.hpp>

        using boost::math::is_float_int;
        using boost::multiprecision::cpp_bin_float_50;
        using boost::multiprecision::cpp_bin_float_100;

        template<class Real>
        void test_one()
        {
            std::cout << "Testing type " << boost::typeindex::type_id<Real>().pretty_name()  << " works with is_float_int\n";

            Real x = 1;
            BOOST_CHECK(is_float_int(x));
            x = 1.5;
            BOOST_CHECK(!is_float_int(x));
        }

        BOOST_AUTO_TEST_CASE(test_is_float_int)
        {
            // Standard floating-point types:
            test_one<float>();
            test_one<double>();
            test_one<long double>();
            // Ensure extended precision types work:
            test_one<boost::math::concepts::real_concept>();
            test_one<cpp_bin_float_50>();
            test_one<cpp_bin_float_100>();
        }

This test must be added to `test/Jamfile.v2` to start working with the CI system:

    .. code: bash

        [ run is_float_int_test.cpp ../../test/build//boost_unit_test_framework : : :  [ requires cxx11_smart_ptr cxx11_defaulted_functions ]  ]


- Write a compile test. A unit test may bring in other header files which define certain functions in the feature file.
But we need to know that the feature file compiles independently of any other file-this is what the compile test attempts to prove:

.. code: cpp

    //  Copyright John Maddock, Nick Thompson 2017.
    //  Use, modification and distribution are subject to the
    //  Boost Software License, Version 1.0. (See accompanying file
    //  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

    #include <boost/math/tools/is_float_int.hpp>
    #include "test_compile_result.hpp"

    void compile_and_link_test()
    {
       bool t = boost::math::tools::is_float_int<double>(1.4);
       check_result<bool>(t;
    }

Place this file in `test/compile_test/is_float_int_incl_test.cpp`.

- Write a concept test.

- Write an example. This is simply a `main` function that should a simple application of your code, preferably with explanatory comments.

.. code: cpp

    //  Copyright John Maddock, Nick Thompson 2017.
    //  Use, modification and distribution are subject to the
    //  Boost Software License, Version 1.0. (See accompanying file
    //  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

    //[is_float_int_example

    /*`This example demonstrates how to use the cubic b spline interpolator for regularly spaced data.
    */
    #include <boost/math/tools/is_float_int.hpp>

    int main()
    {
        double x = 1.2;
        std::cout << "Is 1.2 an integer? ";
        if (boost::math::is_float_int(x)) {
            std::cout << "Yes\n";
        }
        else {
            std::cout << "No\n";
        }

        return 0;
    }
    //] [/is_float_int_example]

Then add the following line to `example/Jamfile.v2`:

.. code:: bash

    [ run is_float_int_example.cpp : : : [ requires cxx11_smart_ptr cxx11_hdr_random cxx11_defaulted_functions ]  ]

- Write your documentation. For instance, in the file `doc/internals/is_float_int.qbk`, we could write

.. code:: bash

    [section:is_float_int Testing if a float is an integer]

    [h4 Synopsis]
    This tests if a function is an integer.

Then add your documentation to `doc/math.qbk` via `[include internals/is_float_int.qbk]`, and make sure your documentation builds via

.. code:: bash

    boost/libs/math/doc$ ../../../b2 release

and then examine the generated html to see that your documentation looks correct.


.. _boost.org: http://www.boost.org/doc/libs/1_65_0/libs/math/doc/html/index.html
