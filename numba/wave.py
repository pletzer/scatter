import numpy
import os
import ctypes

PI = numpy.pi
TWOPI = 2. * PI
FOURPI = 2. * TWOPI

# find the shared library under the build directory
libFile = ''
for root, dirs, files in os.walk('build/'):
    for file in files:
        if file[-3:] == '.so':
            libFile = os.path.join(root, file)

# open the shared library
lib = ctypes.CDLL(libFile)

# tell ctypes how to call the C++ implemented Bessel functions
c_j0 = lib.j0
c_j0.restype = ctypes.c_double
c_j0.argtypes = [ctypes.c_double]

c_y0 = lib.y0
c_y0.restype = ctypes.c_double
c_y0.argtypes = [ctypes.c_double]

c_j1 = lib.j1
c_j1.restype = ctypes.c_double
c_j1.argtypes = [ctypes.c_double]

c_y1 = lib.y1
c_y1.restype = ctypes.c_double
c_y1.argtypes = [ctypes.c_double]


def hankel1_0(x):
    return c_j0(x) + 1j* c_y0(x)


def hankel1_1(x):
    return c_j1(x) + 1j* c_y1(x)


def incident(kvec, point):
    """
    Incident wave

    @param kvec incident wave vector
    @param point target point
    @return complex number
    """
    return numpy.exp(1j*numpy.dot(kvec, point))

def gradIncident(nvec, kvec, point):
    """
    Normal gradient of the incident wave, assumes incident wave is exp(1j * kvec.x)

    @param nvec normal vector pointing inwards
    @param kvec incident wave vector
    @param point (source) point
    @return complex number
    """
    return 1j*numpy.dot(nvec, kvec)*incident(kvec, point)

def computeScatteredWaveElement(kvec, p0, p1, point):
    """
    Scattered wave contribution from a single segment
    @param kvec incident wave vector
    @param p0 starting point of the segment
    @param p1 end point of the segment
    @param point observer point
    @return complex value
    """

    # xdot is anticlockwise
    xdot = p1 - p0

    # mid point of the segment
    pmid = 0.5*(p0 + p1)

    # segment length
    dsdt = numpy.sqrt(numpy.dot(xdot, xdot))

    # normal vector, pointintg inwards and normalised
    nvec = numpy.array([-xdot[1], xdot[0],])
    nvec /= numpy.sqrt(numpy.dot(nvec, nvec))

    # from segment mid-point to observer
    rvec = point - pmid
    r = numpy.sqrt(numpy.dot(rvec, rvec))

    kmod = numpy.sqrt(numpy.dot(kvec, kvec))
    kr = kmod * r

    # Green functions and normal derivatives
    g = (1j/4.) * hankel1_0(kr)
    dgdn = (-1j/4.) * hankel1_1(kr) * kmod * numpy.dot(nvec, rvec) / r

    # contribution from the gradient of the incident wave on the surface
    # of the obstacle. The normal derivative of the scattered wave is 
    # - normal derivative of the incident wave.
    scattered_wave = - dsdt * g * gradIncident(nvec, kvec, pmid)

    # shadow side: total wave is nearly zero 
    #              => scattered wave amplitude = -incident wave ampl.
    #
    # illuminated side:
    #              => scattered wave amplitude = +incident wave ampl.
    shadow = 2*((numpy.dot(nvec, kvec) > 0.) - 0.5) # +1 on the shadow side, -1 on the illuminated side
    scattered_wave += shadow * dsdt * dgdn * incident(kvec, pmid)

    return scattered_wave

def computeScatteredWave(kvec, xc, yc, point):
    """
    Total scattered wave response, summing up 
    contributions from each segment

    @param kvec incident wave vector
    @param xc list of x coordinates representing the contour, must close
    @param yc list of y coordinates representing the contour, must close
    @param point observer point
    @return complex value
    """
    res = 0j
    n = len(xc)
    for i0 in range(n - 1):
        p0 = numpy.array([xc[i0], yc[i0],])
        i1 = i0 + 1
        p1 = numpy.array([xc[i1], yc[i1],])
        res += computeScatteredWaveElement(kvec, p0, p1, point)
    return res

