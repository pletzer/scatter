import argparse
import numpy
import scipy.special
from numpy import cos, sin, pi
import math
import saveVtk
import wave

parser = argparse.ArgumentParser(description='Compute field scattered by an obstacle.')
parser.add_argument('-lambda', dest='lmbda', type=float, default=0.5, help='x wavelength')
parser.add_argument('-nx', dest='nx', type=int, default=128, help='number of x cells')
parser.add_argument('-ny', dest='ny', type=int, default=128, help='number of y cells')
parser.add_argument('-nc', dest='nc', type=int, default=128, help='number of contour segments')
parser.add_argument('-xc', dest='xContourExpr', type=str, default='cos(2*pi*t + 0.5*sin(2*pi*t + 0.9))', help='x contour expression of 0 <= t <= 1')
parser.add_argument('-yc', dest='yContourExpr', type=str, default='sin(2*pi*t)', help='y contour expression of 0 <= t <= 1')
parser.add_argument('-save', dest='save', action='store_true', help='save time varying solution in VTK files')
parser.add_argument('-checksum', dest='checksum', action='store_true', help='compute and print a checksum of the scattered wave')


args = parser.parse_args()

twoPi = 2. * numpy.pi

# incident wavenumber
knum = 2 * numpy.pi / args.lmbda
kvec = numpy.array([knum, 0.,], numpy.float64)

def isInsideContour(p, xc, yc):
    """
    Check if a point is inside closed contour

    @param p point (2d array)
    @param xc array of x points, anticlockwise and must close
    @param yc array of y points, anticlockwise and must close
    @return True if p is inside, False otherwise
    """

    # number of segments
    numSeg = len(xc) - 1

    # vectors from point p to start point of segment
    a = numpy.array([xc[:-1], yc[:-1]])
    a[0, :] -= p[0]
    a[1, :] -= p[1]

    # vectors from point p to end point of segment
    b = numpy.array([xc[1:], yc[1:]])
    b[0, :] -= p[0]
    b[1, :] -= p[1]

    areas = a[0, :]*b[1, :] - a[1, :]*b[0, :]
   
    # return True if all the areas are positive
    return not numpy.any(areas < 1.e-10)

# contour points of the obstacle
t = numpy.linspace(0., 1., args.nc + 1)
xc = eval(args.xContourExpr)
yc = eval(args.yContourExpr)

# create grid 
nx, ny = args.nx, args.ny
xmin, xmax = xc.min() - 5*args.lmbda, xc.max() + 3*args.lmbda
ymin, ymax = yc.min() - 3*args.lmbda, yc.max() + 4*args.lmbda
xg = numpy.linspace(xmin, xmax, nx + 1)
yg = numpy.linspace(ymin, ymax, ny + 1)


# compute the field
scat = numpy.zeros((ny + 1, nx + 1), numpy.complex64)
inci = numpy.zeros((ny + 1, nx + 1), numpy.complex64)
for j in range(ny + 1):
    y = yg[j]
    for i in range(nx + 1):
        x = xg[i]

        # need to check that x,y are outside contour
        # otherwise continue
        p = numpy.array([x, y,])

        # skip if point is inside closed contour
        if isInsideContour(p, xc, yc):
            continue

        inci[j, i] = wave.incident(kvec, p.reshape([1, 2,]))
        scat[j, i] = wave.computeScatteredWave(kvec, xc, yc, p)

if args.checksum:
    print('Sum of scattered field |amplitudes|^2: {}'.format((scat*numpy.conj(scat)).sum().real))

if args.save:
    # number of time frames
    nanim = 20
    dOmegaTime = twoPi / float(nanim)
    for it in range(nanim):
        data = numpy.real(numpy.exp(-1j*it*dOmegaTime) * (scat + inci))
        saveVtk.saveData('scatter_{}.vtk'.format(it), xg, yg, data, 'total')
