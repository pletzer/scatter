import numpy
from scipy.special import hankel1

ZHAT = numpy.array([0., 0., 1.])
PI = numpy.pi
TWOPI = 2. * PI
FOURPI = 2. * TWOPI


def incident(kvec, point):
	"""
	Incident wave

	@param kvec incident wave vector
	@param point target point
	@return complex number
	"""
	return numpy.exp(1j*kvec.dot(point))


def gradIncident(nvec, kvec, point):
	"""
	Normal gradient of the incident wave, assumes incident wave is exp(1j * kvec.x)

	@param nvec normal vector pointing outwards
	@param kvec incident wave vector
	@param point (source) point
	@return complex number
	"""
	return 1j*nvec.dot(kvec)*incident(kvec, point)


def computeScatteredWaveElement(kvec, p0, p1, point):
	"""
	Scattered wave contribution from a single segment
	@param kvec incident wave vector
	@param p0 starting point of the segment
	@param p1 end point of the segment
	@param point observer point
	@return complex value
	"""
	xdot = p1 - p0
	pmid = 0.5*(p0 + p1)

	# segment length
	dsdt = numpy.sqrt(xdot.dot(xdot))

	# normal vector
	nvec = numpy.cross(xdot, ZHAT)
	nvec /= numpy.sqrt(nvec.dot(nvec))

	kmod = numpy.sqrt(kvec.dot(kvec))

	rvec = point - pmid
	r = numpy.sqrt(rvec.dot(rvec))

	# CHECK SIGN!!!
	return +gradIncident(nvec, kvec, pmid) * (1j/4.) * dsdt * \
            hankel1(0, kmod * r)

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
		p0 = numpy.array([xc[i0], yc[i0], 0.])
		i1 = i0 + 1
		p1 = numpy.array([xc[i1], yc[i1], 0.])
		res += computeScatteredWaveElement(kvec, p0, p1, point)
	return res

