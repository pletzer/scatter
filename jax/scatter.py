import argparse
import numpy
import jax
import jax.numpy as jnp
from jax.numpy import cos, sin, pi
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
kvec = jnp.array([knum, 0.,], jnp.float64)

def isInsideContour(p, xc, yc):
    """
    Check if a point is inside closed contour

    @param p point (2d array)
    @param xc array of x points, anticlockwise and must close
    @param yc array of y points, anticlockwise and must close
    @return True if p is inside, False otherwise
    """
    # areas of the triangles extending from the point to each segment, a point
    # is outside if any of them is negative (cross product < 0). A point on the
    # contour is counted as outside (cross product == 0)
    ax, ay = xc[:-1] - p[0], yc[:-1] - p[1]
    bx, by = xc[1:] - p[0], yc[1:] - p[1]
    return jnp.all(ax*by - ay*bx > 1.e-10)


@jax.jit
def computeField(points, xc, yc):
    """
    Compute the incident and scattered fields at all observer points

    @param points observer points, array of size n * 2
    @return incident and scattered fields, complex arrays of size n
    """
    def one_point(p):
        inside = isInsideContour(p, xc, yc)
        # the field is left at zero inside the obstacle
        inci = jnp.where(inside, 0j, wave.incident(kvec, p[None, :])[0])
        scat = jnp.where(inside, 0j, wave.computeScatteredWave(kvec, xc, yc, p))
        return inci, scat
    return jax.vmap(one_point)(points)

# contour points of the obstacle
t = jnp.linspace(0., 1., args.nc + 1)
xc = jnp.asarray(eval(args.xContourExpr))
yc = jnp.asarray(eval(args.yContourExpr))

# create grid 
nx, ny = args.nx, args.ny
xmin, xmax = xc.min() - 5*args.lmbda, xc.max() + 3*args.lmbda
ymin, ymax = yc.min() - 3*args.lmbda, yc.max() + 4*args.lmbda
xg = numpy.linspace(xmin, xmax, nx + 1)
yg = numpy.linspace(ymin, ymax, ny + 1)


# compute the field at all the grid points at once (this also compiles the code)
X, Y = numpy.meshgrid(xg, yg)
points = jnp.array(numpy.stack([X.ravel(), Y.ravel()], axis=1))
inci, scat = computeField(points, jnp.array(xc), jnp.array(yc))
inci = numpy.asarray(inci).astype(numpy.complex64).reshape((ny + 1, nx + 1))
scat = numpy.asarray(scat).astype(numpy.complex64).reshape((ny + 1, nx + 1))

if args.checksum:
    print('Sum of scattered field |amplitudes|^2: {}'.format((scat*numpy.conj(scat)).sum().real))

if args.save:
    # number of time frames
    nanim = 20
    dOmegaTime = twoPi / float(nanim)
    for it in range(nanim):
        data = numpy.real(numpy.exp(-1j*it*dOmegaTime) * (scat + inci))
        saveVtk.saveData('scatter_{}.vtk'.format(it), xg, yg, data, 'total')
