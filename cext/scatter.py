import argparse
import numpy
import scipy.special
from numpy import cos, sin, pi
import math
import saveVtk
import ctypes
import os
import wave

parser = argparse.ArgumentParser(description='Compute field scattered by an obstacle.')
parser.add_argument('-lambda', dest='lmbda', type=float, default=0.5, help='x wavelength')
parser.add_argument('-xmin', dest='xmin', type=float, default=-5.0, help='min x value of domain')
parser.add_argument('-xmax', dest='xmax', type=float, default=+3.0, help='max x value of domain')
parser.add_argument('-ymin', dest='ymin', type=float, default=-5.0, help='min y value of domain')
parser.add_argument('-ymax', dest='ymax', type=float, default=+3.0, help='min y value of domain')
parser.add_argument('-nx', dest='nx', type=int, default=128, help='number of x cells')
parser.add_argument('-ny', dest='ny', type=int, default=128, help='number of y cells')
parser.add_argument('-nc', dest='nc', type=int, default=128, help='number of contour segments')
parser.add_argument('-xc', dest='xContourExpr', type=str, default='cos(2*pi*t + 0.5*sin(2*pi*t + 0.9))', help='x contour expression of 0 <= t <= 1')
parser.add_argument('-yc', dest='yContourExpr', type=str, default='sin(2*pi*t)', help='y contour expression of 0 <= t <= 1')
parser.add_argument('-save', dest='save_numit', type=int, help='save time varying solution in numpy files')
parser.add_argument('-checksum', dest='checksum', action='store_true', help='compute and print a checksum of the scattered wave')

args = parser.parse_args()

twoPi = 2. * numpy.pi

# incident wavenumber
knum = 2 * numpy.pi / args.lmbda
kvec = numpy.array([knum, 0.,], numpy.float64)

# contour points of the obstacle
nc = args.nc
nc1 = nc + 1
t = numpy.linspace(0., 1., nc1)
xc = eval(args.xContourExpr)
yc = eval(args.yContourExpr)
xcPtr = xc.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
ycPtr = yc.ctypes.data_as(ctypes.POINTER(ctypes.c_double))

# create grid 
nx, ny = args.nx, args.ny
xg = numpy.linspace(args.xmin, args.xmax, nx + 1)
yg = numpy.linspace(args.ymin, args.ymax, ny + 1)

# find the library under the build directory
waveLibFile = ''
for root, dirs, files in os.walk('build/'):
    for file in files:
        if file[-3:] == '.so':
            waveLibFile = os.path.join(root, file)

# open the shared library 
wavelib = ctypes.CDLL(waveLibFile)

# create some types for calling C++
doubleStarType = ctypes.POINTER(ctypes.c_double) 

# containers to receive the output values of the C function
realVal, imagVal = ctypes.c_double(0.), ctypes.c_double(0.)

# returns void
wavelib.cincident.restype = None
# double*, double*, double*, double*
wavelib.cincident.argtypes = [doubleStarType, doubleStarType, 
                              doubleStarType, doubleStarType]

# returns void
wavelib.computeScatteredWave.restype = None
# double*, int, double*, double*, double*, double*, double* 
wavelib.computeScatteredWave.argtypes = [doubleStarType,
                                         ctypes.c_int, 
                                         doubleStarType,
                                         doubleStarType,
                                         doubleStarType,
                                         doubleStarType,
                                         doubleStarType]

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
        pPtr = p.ctypes.data_as(ctypes.POINTER(ctypes.c_double))

        # skip if point is inside closed contour
        if wavelib.isInsideContour(pPtr, nc1, xcPtr, ycPtr) == 1:
            continue

        # get the pointers from the numpy arrays
        kvecPtr = kvec.ctypes.data_as(ctypes.POINTER(ctypes.c_double))

        wavelib.cincident(kvecPtr, pPtr, ctypes.byref(realVal), ctypes.byref(imagVal))
        inci[j, i] = realVal.value + 1j*imagVal.value

        xcPtr = xc.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        ycPtr = yc.ctypes.data_as(ctypes.POINTER(ctypes.c_double))

        wavelib.computeScatteredWave(kvecPtr, nc1, xcPtr, ycPtr, pPtr, 
                                     ctypes.byref(realVal), ctypes.byref(imagVal))
        scat[j, i] = realVal.value + 1j*imagVal.value

if args.checksum:
    print('Sum of scattered field |amplitudes|^2: {}'.format((scat*numpy.conj(scat)).sum().real))

if args.save_numit:
    # number of time frames
    dOmegaTime = twoPi / float(args.save_numit)
    for it in range(args.save_numit):
        data = numpy.real(numpy.exp(-1j*it*dOmegaTime) * (scat + inci))
        numpy.save(f'scatter_{it:05d}.npy', data)
