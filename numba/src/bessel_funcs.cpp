#include <boost/math/special_functions/bessel.hpp>

/**
 * Bessel J_0 function
 * 
 * @param  n index
 * @param x argument
 * @return value
 */
extern "C"
double j0(double x) {
    return boost::math::cyl_bessel_j<int, double>(0, x);
}

/**
 * Bessel Y_0 function
 * 
 * @param  n index
 * @param x argument
 * @return value
 */
extern "C"
double y0(double x) {
    return boost::math::cyl_neumann<int, double>(0, x);
}

/**
 * Bessel J_1 function
 * 
 * @param  n index
 * @param x argument
 * @return value
 */
extern "C"
double j1(double x) {
    return boost::math::cyl_bessel_j<int, double>(1, x);
}

/**
 * Bessel Y_1 function
 * 
 * @param  n index
 * @param x argument
 * @return value
 */
extern "C"
double y1(double x) {
    return boost::math::cyl_neumann<int, double>(1, x);
}
