import scipy.stats
print([shape_param.domain for shape_param in scipy.stats.beta._shape_info()])
# [[5e-324, inf], [5e-324, inf]]
print(scipy.stats.beta(a=1e-308, b=5).ppf(.2))
# 0.0
#print(scipy.stats.beta(a=1e-309, b=5).ppf(.2))
# OverflowError: Error in function boost::math::tgamma<d>(d): Overflow Error
# (raises error, but kernel doesn't crash)
print(scipy.stats.beta(a=1e-323, b=5).ppf(.2))
# KERNEL CRASHES with message: `Assertion failed: *p_derivative >= 0, file ..\..\scipy\_lib\boost/boost/math/special_functions/beta.hpp, line 739`
