
#include "test_rational.hpp"

#ifdef BOOST_HAS_LONG_LONG
template void do_test_spots<float, unsigned long long>(float, unsigned long long);
#endif
