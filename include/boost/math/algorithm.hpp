//  (C) Copyright Jeremy William Murphy 2016.

//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_ALGORITHM
#define BOOST_MATH_ALGORITHM

#include <boost/math/algebraic_traits.hpp>
#include <boost/concept/assert.hpp>
#include <boost/concept_check.hpp>
#include <boost/bind.hpp>

#include <functional>
#include <numeric>
#include <utility>

namespace boost {
namespace math {
    namespace detail
    {
        struct not_equal_to_zero
        {
            typedef bool result_type;
            
            // Slightly awkward.
            template <template <typename> class Op, typename T>
            bool operator()(Op<T> const &, T const &x) const
            {
                return x != algebraic_traits<Op, T>::zero_element();
            }
        };
    }
    
    
    template<typename I, typename T, typename BinaryOperation, typename Predicate>
    inline std::pair<T, I>
    accumulate_if(I first, I last, T init, BinaryOperation op, Predicate p)
    {
        BOOST_CONCEPT_ASSERT((InputIterator<I>));
        
        for (; first != last && p(init); ++first )
            init = op ( init, *first );
        return std::make_pair(init, first);
    }
    

    template<typename I, typename T, typename BinaryOperation>
    inline std::pair<T, I>
    accumulate_if(I first, I last, T init, BinaryOperation op)
    {
        using boost::placeholders::_1;
        
        detail::not_equal_to_zero const pred;
        return accumulate_if(first, last, init, op, bind(pred, op, _1));
    }
} // math
} // boost

#endif

