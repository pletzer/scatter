#include "is_inside_contour.h"
#include <cmath>

extern "C" int
isInsideContour(const double p[], int n,
                const double xc[], const double yc[], double tol) {
    
    double a[2], b[2];
    double tot = 0.0;
    for (int i0 = 0; i0 < n - 1; ++i0) {
        int i1 = i0 + 1;
        a[0] = xc[i0] - p[0]; a[1] = yc[i0] - p[1];
        b[0] = xc[i1] - p[0]; b[1] = yc[i1] - p[1];
        tot += atan2(a[0]*b[1] - a[1]*b[0], a[0]*b[0] + a[1]*b[1]);
    }
    tot /= TWOPI;
    return std::abs(tot) > tol? 1: 0;
}