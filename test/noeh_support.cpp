

#include <boost/config.hpp>
#include <iostream>
#include <iomanip>
#include <cstdlib>


#ifdef BOOST_NO_EXCEPTIONS

namespace boost {

   void throw_exception(const std::exception& e)
   {
      std::cout << e.what() << std::endl;
      std::abort();
   }


}

#endif
