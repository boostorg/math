
#include "test_rational.hpp"
#include <boost/math/concepts/real_concept.hpp>

#ifdef BOOST_HAS_LONG_LONG
template void do_test_spots<boost::math::concepts::real_concept, unsigned long long>(boost::math::concepts::real_concept, unsigned long long);
#endif
