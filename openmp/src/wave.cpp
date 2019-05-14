#include "wave.h"
#include <boost/math/special_functions/bessel.hpp>

// Make sure that Boost is thread-safe - we call Boost functions
// from within a parallel (threaded) region
#ifndef BOOST_HAS_THREADS
#error This version of Boost does NOT support multithreading - please rebuild with multithreading enabled
#endif

/**
 * Hankel_n^{(1)} function
 * 
 * @param n index
 * @param x argument
 * @return value
 */
std::complex<double>
hankel1(const int n, const double x) {
    // We requested a Boost library with multithreading support
    double besselJ = boost::math::cyl_bessel_j<int, double>(n, x);
    double besselY = boost::math::cyl_neumann<int, double>(n, x);
    return besselJ + J1*besselY;
}

/**
 * Incident wave
 * 
 * @param kvec incident wave vector
 * @param point target point
 * @return amplitude
 */
std::complex<double> 
incident(const double kvec[], const double point[]) {
    // std::exp should be thread-safe
    return std::exp(J1*(kvec[0]*point[0] + kvec[1]*point[1]));
}

/**
 * Normal gradient of the incident wave, assumes incident wave is exp(1j * kvec.x)
 *
 * @param nvec normal vector pointing inwards
 * @param kvec incident wave vector
 * @param point (source) point
 * @return amplitude
 */
std::complex<double> 
gradIncident(const double nvec[], const double kvec[], 
             const double point[]) {
    return J1*(nvec[0]*kvec[0] + nvec[1]*kvec[1])*incident(kvec, point);
}

/**
 * Scattered wave contribution from a single segment
 *
 * @param kvec incident wave vector
 * @param p0 starting point of the segment
 * @param p1 end point of the segment
 * @param point observer point
 * @return wave contribution
 */
std::complex<double> 
computeScatteredWaveElement(const double kvec[], const double p0[], 
                            const double p1[], const double point[]) {

    // xdot is anticlockwise
    double xdot[] = {p1[0] - p0[0], p1[1] - p0[1]};

    // mid point of the segment
    double pmid[] = {0.5*(p0[0] + p1[0]), 0.5*(p0[1] + p1[1])};

    // segment length
    // std:sqrt should be thread-safe
    double dsdt = std::sqrt(xdot[0]*xdot[0] + xdot[1]*xdot[1]);

    // normal vector, pointintg inwards and normalised
    double nvec[] = {-xdot[1]/dsdt, xdot[0]/dsdt, 0.};

    // from segment mid-point to observer
    double rvec[] = {point[0] - pmid[0], point[1] - pmid[1]};
    double r = std::sqrt(rvec[0]*rvec[0] + rvec[1]*rvec[1]);

    double kmod = std::sqrt(kvec[0]*kvec[0] + kvec[1]*kvec[1]);
    double kr = kmod * r;

    // Green functions and normal derivatives
    std::complex<double> g = (J1/4.) * hankel1(0, kr);
    double nDotR = nvec[0]*rvec[0] + nvec[1]*rvec[1];
    std::complex<double> dgdn = (-J1/4.) * hankel1(1, kr) * kmod * nDotR / r;

    // contribution from the gradient of the incident wave on the surface
    // of the obstacle. The normal derivative of the scattered wave is
    // - normal derivative of the incident wave.
    std::complex<double> scattered_wave = - dsdt * g * gradIncident(nvec, kvec, pmid);

    // shadow side: total wave is nearly zero
    //              => scattered wave amplitude = -incident wave ampl.
    //
    // illuminated side:
    //              => scattered wave amplitude = +incident wave ampl.
    //
    double nDotK = nvec[0]*kvec[0] + nvec[1]*kvec[1];          
    double shadow = (nDotK > 0 ? 1.0 : -1.0);
        
    scattered_wave += shadow * dsdt * dgdn * incident(kvec, pmid);

    return scattered_wave;
}

extern "C" void
cincident (const double kvec[], const double point[], double* real_part, double* imag_part) {
    std::complex<double> res = incident(kvec, point);
    *real_part = res.real();
    *imag_part = res.imag();
}

extern "C" void computeScatteredWave(const double kvec[], int nc, const double xc[], const double yc[], 
                                     const double point[], double* real_part, double* imag_part) {
    // Use primitive data types for the parallel reduction operation, to avoid having to define a
    // custom OpenMP reduction for the "std::complex" data type
    double res_real=0., res_imag=0.;
    // OpenMP pragma defines the parallel region (where threads are spawned), data clauses and reduction
    // Variables that are defined inside the loop are automatically private
    #pragma omp parallel for default(none) shared(nc,xc,yc,kvec,point) reduction(+:res_real,res_imag)
    for (int i = 0; i < nc - 1; ++i) {
        double p0[] = {xc[i], yc[i]};
        double p1[] = {xc[i + 1], yc[i + 1]};
        // res is now just a helper variable and can be defined inside the loop
        // We also need to verify that "computeScatteredWaveElement" is thread-safe,
        // as it will now run concurrently in multiple threads
        std::complex<double> res = computeScatteredWaveElement(kvec, p0, p1, point);
        res_real += res.real();
        res_imag += res.imag();
    }
    *real_part = res_real;
    *imag_part = res_imag;
}

