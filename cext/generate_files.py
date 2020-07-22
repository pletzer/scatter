import random
import numpy
import pandas
from numpy import cos, sin, pi
import math
import os
import ctypes
import wave


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


def computeScatteredField(args):
    """
    Compute the scattered field

    @param xg 1-d array of x grid points
    @param yg 1-d array of y grid points
    @param wavelib ctypes CDLL handle to the C library
    @param kvec wave vector
    @param inci incident wave, array of size ny1, nx1
    @param xc x values of the object's contour
    @param yc y values of the object's contour
    @return array of size ny1, nx1
    """

    # unpack
    xg, yg, kvec, inci, xcfun, ycfun = args


    ny1, nx1 = inci.shape
    p = numpy.array([0., 0.,])
    pPtr = p.ctypes.data_as(ctypes.POINTER(ctypes.c_double))

    # number of contour points 
    nc = 128
    nc1 = nc + 1
    t = numpy.linspace(0., 2.*pi, nc1)

    # contour points
    xc = eval(xcfun)
    yc = eval(ycfun)
    xcPtr = xc.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    ycPtr = yc.ctypes.data_as(ctypes.POINTER(ctypes.c_double))

    # get the pointers from the numpy arrays
    kvecPtr = kvec.ctypes.data_as(ctypes.POINTER(ctypes.c_double))

    # containers to receive the output values of the C function
    realVal, imagVal = ctypes.c_double(0.), ctypes.c_double(0.)


    nc1 = len(xc)

    scat = numpy.zeros((ny1, nx1), numpy.complex64)

    for j in range(ny1):

        y = yg[j]
        for i in range(nx1):

            x = xg[i]

            p[:] = x, y

            # always outside, no need to check...
            # skip if point is inside closed contour
            #if wavelib.isInsideContour(pPtr, nc1, xcPtr, ycPtr) == 1:
            #    continue

            wavelib.cincident(kvecPtr, pPtr, ctypes.byref(realVal), ctypes.byref(imagVal))
            inci[j, i] = realVal.value + 1j*imagVal.value

            wavelib.computeScatteredWave(kvecPtr, nc1, xcPtr, ycPtr, pPtr, 
                                         ctypes.byref(realVal), ctypes.byref(imagVal))

            scat[j, i] = realVal.value + 1j*imagVal.value

    return scat


random.seed(123)

# number of cases/files
ncases = 1000

# wavelength
lmbda = 0.43456

# centre of object
xcmin, xcmax = 1.0, 3.0
ycmin, ycmax = -1.0, 1.4

# radius
amin, amax = 0.1, 0.8

# elongation
kmin, kmax = 0.13, 2.2

# triangularity
dmin, dmax = 0., 0.9

# phase
pmin, pmax = -numpy.pi, numpy.pi


# grid
nx = 127
ny = 127
xmin, xmax = -10., 0.
ymin, ymax = -5., 5.
xg = numpy.linspace(xmin, xmax, nx + 1)
yg = numpy.linspace(ymin, ymax, ny + 1)

twoPi = 2. * numpy.pi
# incident wavenumber
knum = 2 * numpy.pi / lmbda
kvec = numpy.array([knum, 0.,], numpy.float64)


# field
inci = numpy.zeros((ny + 1, nx + 1), numpy.complex64)

dic_data = {
    'xcentre': [xcmin + (xcmax - xcmin)*random.random() for i in range(ncases)],
    'ycentre': [ycmin + (ycmax - ycmin)*random.random() for i in range(ncases)],
    'a': [amin + (amax - amin)*random.random() for i in range(ncases)],
    'k': [kmin + (kmax - kmin)*random.random() for i in range(ncases)],
    'd': [dmin + (dmax - dmin)*random.random() for i in range(ncases)],
    'phase':[pmin + (pmax - pmin)*random.random() for i in range(ncases)],
    'id': [i for i in range(ncases)]
}

# random time
omegaTimes = [0. + (2*pi - 0.)*random.random() for i in range(ncases)]

input_values = [(xg, yg, kvec, inci, 
               f'{dic_data["a"][i]} * cos(t) + {dic_data["xcentre"][i]}',
               f'{dic_data["a"][i]} * {dic_data["k"][i]} * sin(t - {dic_data["d"][i]}*sin(t - {dic_data["phase"][i]})) + {dic_data["ycentre"][i]}',)
                 for i in range(ncases)]


for i in range(ncases):
    print(i)
    # compute the scatted field
    scat = computeScatteredField(input_values[i])
    # add to the incident field, evaluate at a random time and take the real part
    numpy.save(f'scatter_{i:05d}.npy', numpy.real( numpy.exp(-1j*omegaTimes[i]) * (scat + inci) ) )



# create dataframe
for key in 'xcentre', 'ycentre', 'a', 'k', 'd', 'phase':
    dic_data[key] = numpy.array(dic_data[key], numpy.float32)
dic_data['id'] = numpy.array(dic_data['id'], numpy.int32)
df = pandas.DataFrame(dic_data)
df.to_csv('scatter.csv')
