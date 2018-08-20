import argparse
import numpy
import scipy.special
from matplotlib import pylab

parser = argparse.ArgumentParser(description='Compute field scattered by an obstacle.')
parser.add_argument('-n', dest='n', ctype=int, help='number of segments')
parser.add_argument('-cx', dest='cx', ctype=str, help='x contour expression of 0 <= t <= 1')
parser.add_argument('-cy', dest='cy', ctype=str, help='y contour expression of 0 <= t <= 1')
parser.add_argument('-k', dest='k', ctype=float, help='x wavenumber')

args = parser.parse_args()

kvec = numpy.array([args.k + 0j, 0. + 0j])

def incidentPsi(p):
	return numpy.exp(1j* kvec.dot(p))

def incidentNormalGradPsi(p, nvec):
	return 1j * kvec.dot(nvec) * incidentPsi(p)

def normalVector(p0, p1):
	dx, dy = p1 - p0
	return numpy.array([dy, -dx])

# create grid
nx, ny = 10, 10
xmin, xmax = -10., 10.
ymin, ymax = -10., 10.
xg = numpy.linspace(xmin, xmax, nx + 1)
yg = numpy.linspace(ymin, ymax, ny + 1)

t = numpy.linspace(0., 1., args.n + 1)
xc = eval(args.cx)
yc = eval(args.cy)

# compute the field
dt = 1.0 / float(args.n)
psi_scattered = numpy.zeros((ny + 1, nx + 1), numpy.complex64)
psi_incident = numpy.zeros((ny + 1, nx + 1), numpy.complex64)
for j in range(ny + 1):
	y = yg[j]
	for i in range(nx + 1):
		x = xg[i]
		# need to check that x,y are outside contour
		# otherwise continue
		p = numpy.array([x, y])
		for k in range(args.n):
			t = k * dt
			p0 = numpy.array([eval(args.cx), eval(args.cy)])
			t += dt
			p1 = numpy.array([eval(args.cx), eval(args.cy)])
			pm = 0.5*(p0 + p1)
			d = p - pm
			r = numpy.sqrt(d.dot(d))
			nvec = normalVector(p0, p1)
			psi_incident[j, i] = incidentNormalGradPsi(pm, nvec)
			psi_scattered[j, i] += (-1j/4.) * scipy.special.hankel1(0, args.k * r) * psi_incident[j, i]

# plot
pylab.pcolot(xg, yg, psi_incident.real + psi_scattered.real)
pylab.show()