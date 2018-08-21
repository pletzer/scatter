import argparse
import numpy
import scipy.special
from matplotlib import pyplot
from numpy import cos, sin, pi

parser = argparse.ArgumentParser(description='Compute field scattered by an obstacle.')
parser.add_argument('-n', dest='n', type=int, default=100, help='number of segments')
parser.add_argument('-xc', dest='xContourExpr', type=str, default='cos(2*pi*t)', help='x contour expression of 0 <= t <= 1')
parser.add_argument('-yc', dest='yContourExpr', type=str, default='sin(2*pi*t)', help='y contour expression of 0 <= t <= 1')
parser.add_argument('-lambda', dest='lmbda', type=float, default=0.5, help='x wavenumber')
parser.add_argument('-plot', dest='plot', action='store_false', help='plot solution')

args = parser.parse_args()

# incident wavenumber
knum = 2 * numpy.pi / args.lmbda
kvec = numpy.array([knum, 0.], numpy.float64)
ikvec = numpy.array([knum*1j, 0.], numpy.complex64)

def incidentPsi(p):
	return numpy.exp(ikvec.dot(p))

def incidentNormalGradPsi(p, nvec):
	return ikvec.dot(nvec) * incidentPsi(p)

def normalVector(p0, p1):
	dx, dy = p1 - p0
	return numpy.array([dy, -dx], numpy.float64) / numpy.sqrt(dx*dx + dy*dy)

def greenFunction(pa, pb):
	# distance cannot be zero
	dp = pb - pa
	r = max(numpy.sqrt(numpy.dot(dp, dp)), 1.e-8)
	return 1j*numpy.pi * scipy.special.hankel1(0, knum * r)

def normalGradientGreenFunction(nvec, pa, pb):
	# distance cannot be zero
	dp = pb - pa
	r = max(numpy.sqrt(numpy.dot(dp, dp)), 1.e-8)
	gradR = dp / r
	return -1j*numpy.pi * scipy.special.hankel1(1, knum * r) * knum * nvec.dot(gradR)

# create grid
nx, ny = 100, 100
xmin, xmax = -20*args.lmbda, 10*args.lmbda
ymin, ymax = -10*args.lmbda, 20*args.lmbda
xg = numpy.linspace(xmin, xmax, nx + 1)
yg = numpy.linspace(ymin, ymax, ny + 1)

# contour points of the obstacle
t = numpy.linspace(0., 1., args.n + 1)
xc = eval(args.xContourExpr)
yc = eval(args.yContourExpr)

# compute the field
psi_scattered = numpy.zeros((ny + 1, nx + 1), numpy.complex64)
psi_incident = numpy.zeros((ny + 1, nx + 1), numpy.complex64)
fourPi = 4. * numpy.pi
for j in range(ny + 1):
	y = yg[j]
	for i in range(nx + 1):
		x = xg[i]
		# need to check that x,y are outside contour
		# otherwise continue
		p = numpy.array([x, y])

		# value of the incident wave on the observation point
		psi_incident[j, i] = incidentPsi(p)

		for el0 in range(args.n):
			el1 = el0 + 1

			# segment endpoints
			p0 = numpy.array([xc[el0], yc[el0]], numpy.float64)
			p1 = numpy.array([xc[el1], yc[el1]], numpy.float64)

			# mid point of the segment
			pm = 0.5*(p0 + p1)

			# normal vector, p0 -> p1 in anticlockwise direction 
			nvec = normalVector(p0, p1)

			# vector from mid segment to observation point
			rho = p - pm

			# segment length
			delta = p1 - p0
			ds = numpy.sqrt(delta.dot(delta))

			# normal gradient of the incident wave on the segment mid point
			normalGradIncidentPsiOnObstacle = incidentNormalGradPsi(pm, nvec)

			psi_scattered[j, i] += (ds/fourPi) * greenFunction(pm, p) * normalGradIncidentPsiOnObstacle
			if nvec.dot(kvec) > 0:
				# shadow
				# psi is approximately zero
				# psi_scattered ~ - psi_incident				
				psi_scattered[j, i] -= (ds/fourPi) * normalGradientGreenFunction(nvec, pm, p) * incidentPsi(pm)


if args.plot:
	# plot
	#pyplot.contour(xg, yg, psi_scattered.real, levels=numpy.linspace(-2., 2., 11), linestyles='dashed')
	#p = pyplot.contour(xg, yg, numpy.real(psi_incident + psi_scattered), levels=numpy.linspace(-2., 2., 11))
	p = pyplot.pcolor(xg, yg, numpy.real(psi_scattered + psi_incident), cmap='coolwarm', vmin=-2., vmax=2.)
	#p = pyplot.pcolor(xg, yg, numpy.real(psi_incident), cmap='coolwarm')
	pyplot.colorbar(p)

	# add the obstacle contour
	pyplot.plot(xc, yc, 'k-')
	pyplot.axes().set_aspect('equal')
	pyplot.show()