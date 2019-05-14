#include "is_inside_contour.h"

extern "C" int
isInsideContour(const double p[], int n,
                const double xc[], const double yc[]) {
    
    bool inside = true;
    for (int i0 = 0; i0 < n - 1; ++i0) {
        int i1 = i0 + 1;
        double a[] = {xc[i0] - p[0], yc[i0] - p[1]};
        double b[] = {xc[i1] - p[0], yc[i1] - p[1]};
        inside &= (a[0]*b[1] - a[1]*b[0] > 1.e-10);
    }
    return (inside? 1: 0);
}
