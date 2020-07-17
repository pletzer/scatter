import random
import numpy
import pandas
from numpy import cos, sin, pi
import math
import os
import ctypes
import wave

random.seed(123)

# number of cases/files
n = 10

# wavelength
lmbda = 0.2234

# radius
amin, amax = 0.1, 1.2

# elongation
kmin, kmax = 0.1, 2.3

# triangularity
dmin, dmax = 0., 0.9

# phase
pmin, pmax = -numpy.pi, numpy.pi

# number of contour points 
nc = 128
nc1 = nc + 1
t = numpy.linspace(0., 1., nc1)

# grid
nx = 127
ny = 127
xmin, xmax = -5., +3.
ymin, ymax = -5., 5.
xg = numpy.linspace(xmin, xmax, nx + 1)
yg = numpy.linspace(ymin, ymax, ny + 1)

dic_data = {
    'a': [],
    'k': [],
    'd': [],
    'phase':[],
    'id': []
}

twoPi = 2. * numpy.pi
# incident wavenumber
knum = 2 * numpy.pi / lmbda
kvec = numpy.array([knum, 0.,], numpy.float64)

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

# field
scat = numpy.zeros((ny + 1, nx + 1), numpy.complex64)
inci = numpy.zeros((ny + 1, nx + 1), numpy.complex64)

for it in range(n):

    a = amin + (amax - amin)*random.random()
    k = kmin + (kmax - kmin)*random.random()
    d = dmin + (dmax - dmin)*random.random()
    phase = pmin + (pmax - pmin)*random.random()

    print(f'iter = {it:06d} a = {a:6.4f} k = {k:6.4f} d = {d:6.4f} phase = {phase:6.3f}')
    dic_data['a'].append(a)
    dic_data['k'].append(k)
    dic_data['d'].append(d)
    dic_data['phase'].append(phase)
    dic_data['id'].append(it)

    xc = eval(f'{a}*cos(2*pi*(t - {phase}))')
    yc = eval(f'{a}*{k}*sin(2*pi*(t - {phase}) - {d}*sin(2*pi*(t - {phase})))')
    xcPtr = xc.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    ycPtr = yc.ctypes.data_as(ctypes.POINTER(ctypes.c_double))


    # compute the field
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


    # random wave phase
    OmegaTime = 0. + (2*pi - 0.)*random.random()
    data = numpy.real(numpy.exp(-1j*OmegaTime) * (scat + inci))
    numpy.save(f'scatter_{it:05d}.npy', data)

# create dataframe
for key in 'a', 'k', 'd', 'phase':
    dic_data[key] = numpy.array(dic_data[key], numpy.float32)
dic_data['id'] = numpy.array(dic_data['id'], numpy.int32)
df = pandas.DataFrame(dic_data)
df.to_csv('scatter.csv')
