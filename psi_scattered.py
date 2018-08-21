import numpy
from scipy.special import hankel1

ZHAT = numpy.array([0., 0., 1.])
PI = numpy.pi
TWOPI = 2. * PI
FOURPI = 2. * TWOPI



def integrate(f, a, b):
	"""
	Definite integral using mid-point rule

	@param f function to integrate
	@param a lower bound
	@param b upper bound
	"""

	return f(0.5*(a + b)) * (b - a)

def gradPsiIncident(nvec, kvec, point):
	"""
	Normal gradient of the incident wave, assumes incident wave is exp(1j * kvec.x)

	@param nvec normal vector pointing outwards
	@param kvec incident wave vector
	@param point (source) point where gradPsiIncident is to be evaluated
	@return complex number
	"""
	return nvec.dot(1j*kvec*numpy.exp(1j*kvec.dot(point)))


def computePsiScatteredElement(kvec, p0, p1, point):
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
	ds = numpy.sqrt(xdot.dot(xdot))

	# normal vector
	nvec = numpy.cross(xdot, ZHAT)
	nvec /= numpy.sqrt(nvec.dot(nvec))

	kmod = numpy.sqrt(kvec.dot(kvec))

	def integrand(t):
		# vector from source to observer
		rvec = point - (p0 + xdot*t)
		r = numpy.sqrt(rvec.dot(rvec))
		return hankel1(0, kmod * r)

	return gradPsiIncident(nvec, kvec, pmid) * 1j * PI * ds * \
	       integrate( integrand, 0., 1.) / FOURPI


def computePsiScattered(kvec, xc, yc, point):
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
		print computePsiScatteredElement(kvec, p0, p1, point)
		res += computePsiScatteredElement(kvec, p0, p1, point)
	return res

