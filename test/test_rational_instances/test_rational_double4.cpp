
#include "test_rational.hpp"

#ifdef BOOST_HAS_LONG_LONG
template void do_test_spots<double, unsigned long long>(double, unsigned long long);
#endif
