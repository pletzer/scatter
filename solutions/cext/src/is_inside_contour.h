#include <cmath>
#include <complex>

#ifndef IS_INSIDE_CONTOUR_H
#define IS_INSIDE_CONTOUR_H

/**
 * Check is a point is inside closed contour
 *
 * @param p point (2d array)
 * @param n number of points
 * @param xc array of x points, anticlockwise and must close
 * @param yc array of y points, anticlockwise and must close
 * @return 1 if p is inside, 0 otherwise
 */
extern "C" int
isInsideContour(const double p[], int n, const double xc[], const double yc[]);

#endif // IS_INSIDE_CONTOUR_H
