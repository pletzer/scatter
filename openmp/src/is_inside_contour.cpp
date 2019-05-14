#include "is_inside_contour.h"

extern "C" int
isInsideContour(const double p[], int n,
                const double xc[], const double yc[]) {
    
    // find the minimum area between the point and the segment. If any
    // area is negative then the point is outside the contour
    double inside = 1.;
    #pragma omp parallel for default(none) shared(n,xc,yc,p) reduction(min:inside)
    for (int i0 = 0; i0 < n - 1; ++i0) {
        int i1 = i0 + 1;
        double a[] = {xc[i0] - p[0], yc[i0] - p[1]};
        double b[] = {xc[i1] - p[0], yc[i1] - p[1]};
        inside = std::min(inside, a[0]*b[1] - a[1]*b[0]);
    }
    return (inside > 1.e-10 ? 1 : 0);
}
