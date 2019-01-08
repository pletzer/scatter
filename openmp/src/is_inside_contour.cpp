#include "is_inside_contour.h"

extern "C" int
isInsideContour(const double p[], int n,
                const double xc[], const double yc[], const double tol) {
    
    double tot = 0.0;
    #pragma omp parallel for default(none) shared(n,xc,yc,p) reduction(+:tot)
    for (int i0 = 0; i0 < n - 1; ++i0) {
        int i1 = i0 + 1;
        double a[] = {xc[i0] - p[0], yc[i0] - p[1]};
        double b[] = {xc[i1] - p[0], yc[i1] - p[1]};
        tot += std::atan2(a[0]*b[1] - a[1]*b[0], a[0]*b[0] + a[1]*b[1]);
    }
    tot /= TWOPI;
    return std::abs(tot) > tol ? 1 : 0;
}
