#include <cmath>
#include <complex>

#ifndef WAVE_H
#define WAVE_H

#define TWOPI 2. * M_PI
#define FOURPI 2. * TWOPI
#define J1 std::complex<double>(0., 1.)


/**
 * Incident wave (C interface)
 * 
 * @param kvec incident wave vector
 * @param point target point
 * @return amplitude
 * @param real_part real part of the output
 * @param imag_part imaginary part of the output
 */
extern "C" void
cincident(const double kvec[], const double point[], 
          double* real_part, double* imag_part);

/**
 * Total scattered wave response, summing up contributions from each segment
 * (C interface)
 *
 * @param kvec incident wave vector
 * @param nc number of contour points
 * @param xc list of x coordinates representing the contour, must close
 * @param yc list of y coordinates representing the contour, must close
 * @param point observer point
 * @param real_part real part of the output
 * @param imag_part imaginary part of the output
 * @return wave response
*/
extern "C" void
computeScatteredWave(const double kvec[], int nc, const double xc[], const double yc[], 
                     const double point[], double* real_part, double* imag_part);

#endif // WAVE_H
