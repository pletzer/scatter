import argparse
import numpy
import scipy.special
from numpy import cos, sin, pi
import math
import saveVtk
import ctypes
import glob
import os

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

# contour points of the obstacle
nc = args.nc
nc1 = nc + 1
t = numpy.linspace(0., 1., nc1)
xc = eval(args.xContourExpr)
yc = eval(args.yContourExpr)

# create grid 
nx, ny = args.nx, args.ny
xmin, xmax = xc.min() - 5*args.lmbda, xc.max() + 3*args.lmbda
ymin, ymax = yc.min() - 3*args.lmbda, yc.max() + 4*args.lmbda
xg = numpy.linspace(xmin, xmax, nx + 1)
yg = numpy.linspace(ymin, ymax, ny + 1)

# find the shared library under the build directory next to this file, so that the
# script works from any directory
_buildDir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'build')
_libFiles = sorted(glob.glob(os.path.join(_buildDir, '**', '*.so'), recursive=True))
if not _libFiles:
    raise FileNotFoundError('no shared library found under {}: '
                            'build it first with "python setup.py build"'.format(_buildDir))
waveLibFile = _libFiles[-1]

# open the shared library 
wavelib = ctypes.CDLL(waveLibFile)

# create some types for calling C++
doubleArr = numpy.ctypeslib.ndpointer(dtype=numpy.float64, flags='C_CONTIGUOUS')

# returns void
wavelib.computeField.restype = None
# double[], int, double[], double[], int, int, double[], double[], double[], double[]
wavelib.computeField.argtypes = [doubleArr, ctypes.c_int, doubleArr, doubleArr,
                                 ctypes.c_int, ctypes.c_int, doubleArr, doubleArr,
                                 doubleArr, doubleArr]

# compute the field at all the grid points in one call, the C++ code loops over
# the observer points in parallel (OpenMP). The (complex) results are returned as
# interleaved real and imaginary parts
nx1, ny1 = nx + 1, ny + 1
inciC = numpy.zeros((ny1, nx1), numpy.complex128)
scatC = numpy.zeros((ny1, nx1), numpy.complex128)
wavelib.computeField(kvec, nc1, xc, yc, nx1, ny1, xg, yg, inciC.view(numpy.float64), scatC.view(numpy.float64))
inci = inciC.astype(numpy.complex64)
scat = scatC.astype(numpy.complex64)

if args.checksum:
    print('Sum of scattered field |amplitudes|^2: {}'.format((scat*numpy.conj(scat)).sum().real))

if args.save:
    # number of time frames
    nanim = 20
    dOmegaTime = twoPi / float(nanim)
    for it in range(nanim):
        data = numpy.real(numpy.exp(-1j*it*dOmegaTime) * (scat + inci))
        saveVtk.saveData('scatter_{}.vtk'.format(it), xg, yg, data, 'total')
