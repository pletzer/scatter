#include <cmath>
#include <complex>

#ifndef IS_INSIDE_CONTOUR_H
#define IS_INSIDE_CONTOUR_H

#define TWOPI 2. * M_PI

/**
 * Check is a point is inside closed contour by summing the 
 * the angles between point p, (xc[i], yc[i]) and (xc[i+1], yc[i+1]).
 * Point p id declared to be inside if the total angle amounts to 
 * 2*pi.
 *
 * @param p point (2d array)
 * @param n number of points
 * @param xc array of x points, anticlockwise and must close
 * @param yc array of y points, anticlockwise and must close
 * @param tol tolerance
 * @return 1 if p is inside, 0 otherwise
 */
extern "C" int
isInsideContour(const double p[], int n, const double xc[], const double yc[], double tol);

#endif // IS_INSIDE_CONTOUR_H
