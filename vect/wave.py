import numpy
from scipy.special import hankel1

PI = numpy.pi
TWOPI = 2. * PI
FOURPI = 2. * TWOPI


def incident(kvec, points):
    """
    Incident wave

    @param kvec incident wave vector
    @param point target points, of size n * 3
    @return complex number
    """
    n = points.shape[0]
    return numpy.exp(1j*points.dot(kvec))


def gradIncident(kvec, nDotK, points):
    """
    Normal gradient of the incident wave, assumes incident wave is exp(1j * kvec.x)

    @param kvec incident wave vector
    @param nDotK nvec . kvec array of size n
    @param points (source) point, array of size n * 3
    @return complex number
    """
    return 1j*nDotK*incident(kvec, points)


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

    # |kvec|
    kmod = numpy.sqrt(kvec.dot(kvec))

    # number of points
    n = len(xc)

    # number of segments
    nm1 = n - 1

    # contour points as an n * 3 array
    pc = numpy.array([[xc[i], yc[i],] for i in range(n)])

    # vectors from one point to the next, nm1 * 3 array
    xdot = pc[1:, :] - pc[:-1, :]

    # mid-segment positions, array of size nm1 * 3
    pmid = 0.5*(pc[1:, :] + pc[:-1, :])

    # segment lengths, array of size nm1
    dsdt = numpy.sqrt(xdot[:, 0]*xdot[:, 0] + xdot[:, 1]*xdot[:, 1])

    # normal vectors at the mid segment locations, array of size nm1 * 3
    nvec = numpy.array([[-xdot[i, 1]/dsdt[i], xdot[i, 0]/dsdt[i]] for i in range(nm1)])

    # vector from the mid-segment positions to the observer, array of size nm1 * 3
    rvec = numpy.array([point - pmid[i, :] for i in range(nm1)])

    # distance from mid-segment positions to observer, array of size nm1
    r = numpy.sqrt(rvec[:, 0]*rvec[:, 0] + rvec[:, 1]*rvec[:, 1])

    kr = kmod * r
    nDotR = numpy.array([nvec[i, :].dot(rvec[i, :]) for i in range(nm1)])
    nDotK = numpy.array([nvec[i, :].dot(kvec) for i in range(nm1)])

    g =  0.25j * hankel1(0, kr)
    dgdn = -0.25j * kmod * nDotR * hankel1(1, kr) / r

    # shadow side: total wave is nearly zero 
    #              => scattered wave amplitude = -incident wave ampl.
    #
    # illuminated side:
    #              => scattered wave amplitude = +incident wave ampl.
    shadow = 2*((nDotK > 0.) - 0.5)

    # contribution from the gradient of the incident wave on the surface
    # of the obstacle. The normal derivative of the scattered wave is 
    # - normal derivative of the incident wave.
    return numpy.sum( dsdt * (-g * gradIncident(kvec, nDotK, pmid) + \
                              shadow * dgdn * incident(kvec, pmid)) )

