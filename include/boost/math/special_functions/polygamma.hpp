
///////////////////////////////////////////////////////////////////////////////
//  Copyright 2013 Nikhar Agrawal
//  Copyright 2013 Christopher Kormanyos
//  Copyright 2013 John Maddock
//  Copyright 2013 Paul Bristow
//  Distributed under the Boost
//  Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef _BOOST_POLYGAMMA_2013_07_30_HPP_
  #define _BOOST_POLYGAMMA_2013_07_30_HPP_

  #include <boost/array.hpp>
  #include <boost/cstdint.hpp>
  #include <boost/math/special_functions/factorials.hpp>
  #include "detail/polygamma.hpp"

  namespace boost { namespace math {

  template<class T>	  
  struct promoteftod
  {
	  typedef T type;
  };  

  template<>
  struct promoteftod<float>
  {
	  typedef double type;
  };

  template<class T, class Policy>
  inline T polygamma(const int n, T x, const Policy &pol)
  {
	typedef typename promoteftod<T>::type result_type;
//	std::cout<<"~:"<<typeid(T).name()<<std::endl;
//	std::cout<<"~:"<<typeid(result_type).name()<<std::endl;
	result_type xx=result_type(x);
        result_type result= boost::math::detail::polygamma_imp(n,xx,pol);
	return T(result);
  }

  template<class T>
  inline T polygamma(const int n, T x)
  {
      return boost::math::polygamma(n,x,policies::policy<>());
  }

  template<class T, class Policy>
  inline T digamma(T x, const Policy &pol)
  {
      return boost::math::polygamma(0,x,pol);
  }

  template<class T>
  inline T digamma(T x)
  {
      return boost::math::digamma(x,policies::policy<>());
  }

  template<class T, class Policy>
  inline T trigamma(T x, const Policy &pol)
  {
      return boost::math::polygamma(1,x,pol);
  }

  template<class T>
  inline T trigamma(T x)
  {
      return boost::math::trigamma(x,policies::policy<>());
  }


} } // namespace boost::math

#endif // _BOOST_BERNOULLI_2013_05_30_HPP_

