import argparse
import numpy
import scipy.special
from numpy import cos, sin, pi
import math
import saveVtk
import wave
from numba import jit

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

@jit(nopython=True)
def isInsideContour(p, xc, yc):
    """
    Check if a point is inside closed contour

    @param p point (2d array)
    @param xc array of x points, anticlockwise and must close
    @param yc array of y points, anticlockwise and must close
    @return True if p is inside, False otherwise
    """
    inside = True
    for i0 in range(len(xc) - 1):
        i1 = i0 + 1
        a = numpy.array([xc[i0], yc[i0]]) - p[:2]
        b = numpy.array([xc[i1], yc[i1]]) - p[:2]
        # point is outside if any of the triangle extending from the point to the segment
        # has negative area (cross product < 0)
        # count a point on the contour as being outside (cross product == 0)
        inside &= (a[0]*b[1] - a[1]*b[0] > 1.e-10)
    return inside

# contour points of the obstacle
t = numpy.linspace(0., 1., args.nc + 1)
xc = eval(args.xContourExpr)
yc = eval(args.yContourExpr)

# create grid 
nx, ny = args.nx, args.ny
xmin, xmax = xc.min() - 5*args.lmbda, xc.max() + 3*args.lmbda
ymin, ymax = yc.min() - 3*args.lmbda, yc.max() + 4*args.lmbda
ny1, nx1 = ny + 1, nx + 1
xg = numpy.linspace(xmin, xmax, nx1)
yg = numpy.linspace(ymin, ymax, ny1)

@jit(nopython=True)
def computeField(flatIndx):
    """
    Compute the incident and scatted fields

    @param fltIndx flat index
    """

    # get the i j indices from the flat index
    j = flatIndx // nx1
    i = flatIndx % nx1

    # get the point coordinates
    x, y = xg[i], yg[j]

    # need to check that x,y are outside contour
    # otherwise continue
    p = numpy.array([x, y,])

    # skip if point is inside closed contour
    if isInsideContour(p, xc, yc):
        return (0j, 0j)
    else:
        inci_val = wave.incident(kvec, p)
        scat_val = wave.computeScatteredWave(kvec, xc, yc, p)
        return (inci_val, scat_val)

res = [computeField(k) for k in range(ny1 * nx1)]

# compute the field
inci = numpy.array([r[0] for r in res], numpy.complex64).reshape((ny1, nx1))
scat = numpy.array([r[1] for r in res], numpy.complex64).reshape((ny1, nx1))


if args.checksum:
    print('Sum of scattered field |amplitudes|^2: {}'.format((scat*numpy.conj(scat)).sum().real))

if args.save:
    # number of time frames
    nanim = 20
    dOmegaTime = twoPi / float(nanim)
    for it in range(nanim):
        data = numpy.real(numpy.exp(-1j*it*dOmegaTime) * (scat + inci))
        saveVtk.saveData('scatter_{}.vtk'.format(it), xg, yg, data, 'total')
