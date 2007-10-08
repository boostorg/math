
#include "test_rational.hpp"

#ifdef BOOST_HAS_LONG_LONG
template void do_test_spots<long double, unsigned long long>(long double, unsigned long long);
#endif
