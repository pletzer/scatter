import jax
import jax.numpy as jnp

# the reference results are computed in double precision
jax.config.update('jax_enable_x64', True)

PI = jnp.pi
TWOPI = 2. * PI
FOURPI = 2. * TWOPI

# JAX has no Bessel functions of the second kind (and no J_0, J_1 either), so
# we use the rational/trigonometric approximations of Abramowitz & Stegun
# (as given in Numerical Recipes), accurate to ~1e-8 for real x > 0. This is
# more than enough for the single precision results of this code.

_TWO_OVER_PI = 0.636619772


def _poly(y, coeffs):
    """Horner evaluation of c[0] + y*(c[1] + y*(c[2] + ...))"""
    res = coeffs[-1]
    for c in coeffs[-2::-1]:
        res = c + y*res
    return res


def _j0(x):
    y = x*x
    small = _poly(y, [57568490574.0, -13362590354.0, 651619640.7, -11214424.18,
                      77392.33017, -184.9052456]) / \
            _poly(y, [57568490411.0, 1029532985.0, 9494680.718, 59272.64853,
                      267.8532712, 1.0])
    z = 8. / x
    y = z*z
    xx = x - 0.785398164
    a1 = _poly(y, [1.0, -0.1098628627e-2, 0.2734510407e-4, -0.2073370639e-5,
                   0.2093887211e-6])
    a2 = _poly(y, [-0.1562499995e-1, 0.1430488765e-3, -0.6911147651e-5,
                   0.7621095161e-6, -0.934935152e-7])
    large = jnp.sqrt(_TWO_OVER_PI/x)*(jnp.cos(xx)*a1 - z*jnp.sin(xx)*a2)
    return jnp.where(x < 8., small, large)


def _y0(x):
    y = x*x
    small = _poly(y, [-2957821389.0, 7062834065.0, -512359803.6, 10879881.29,
                      -86327.92757, 228.4622733]) / \
            _poly(y, [40076544269.0, 745249964.8, 7189466.438, 47447.26470,
                      226.1030244, 1.0]) + _TWO_OVER_PI*_j0(x)*jnp.log(x)
    z = 8. / x
    y = z*z
    xx = x - 0.785398164
    a1 = _poly(y, [1.0, -0.1098628627e-2, 0.2734510407e-4, -0.2073370639e-5,
                   0.2093887211e-6])
    a2 = _poly(y, [-0.1562499995e-1, 0.1430488765e-3, -0.6911147651e-5,
                   0.7621095161e-6, -0.934935152e-7])
    large = jnp.sqrt(_TWO_OVER_PI/x)*(jnp.sin(xx)*a1 + z*jnp.cos(xx)*a2)
    return jnp.where(x < 8., small, large)


def _j1(x):
    y = x*x
    small = x*_poly(y, [72362614232.0, -7895059235.0, 242396853.1, -2972611.439,
                        15704.48260, -30.16036606]) / \
            _poly(y, [144725228442.0, 2300535178.0, 18583304.74, 99447.43394,
                      376.9991397, 1.0])
    z = 8. / x
    y = z*z
    xx = x - 2.356194491
    a1 = _poly(y, [1.0, 0.183105e-2, -0.3516396496e-4, 0.2457520174e-5,
                   -0.240337019e-6])
    a2 = _poly(y, [0.04687499995, -0.2002690873e-3, 0.8449199096e-5,
                   -0.88228987e-6, 0.105787412e-6])
    large = jnp.sqrt(_TWO_OVER_PI/x)*(jnp.cos(xx)*a1 - z*jnp.sin(xx)*a2)
    return jnp.where(x < 8., small, large)


def _y1(x):
    y = x*x
    small = x*_poly(y, [-0.4900604943e13, 0.1275274390e13, -0.5153438139e11,
                        0.7349264551e9, -0.4237922726e7, 0.8511937935e4]) / \
            _poly(y, [0.2499580570e14, 0.4244419664e12, 0.3733650367e10,
                      0.2245904002e8, 0.1020426050e6, 0.3549632885e3, 1.0]) \
            + _TWO_OVER_PI*(_j1(x)*jnp.log(x) - 1.0/x)
    z = 8. / x
    y = z*z
    xx = x - 2.356194491
    a1 = _poly(y, [1.0, 0.183105e-2, -0.3516396496e-4, 0.2457520174e-5,
                   -0.240337019e-6])
    a2 = _poly(y, [0.04687499995, -0.2002690873e-3, 0.8449199096e-5,
                   -0.88228987e-6, 0.105787412e-6])
    large = jnp.sqrt(_TWO_OVER_PI/x)*(jnp.sin(xx)*a1 + z*jnp.cos(xx)*a2)
    return jnp.where(x < 8., small, large)


def _safe(x):
    # both branches of jnp.where are evaluated, keep x away from 0 so that
    # neither produces inf/nan (this matters for gradients, not for values)
    return jnp.maximum(x, 1.e-12)


def hankel1_0(x):
    x = _safe(x)
    return _j0(x) + 1j*_y0(x)


def hankel1_1(x):
    x = _safe(x)
    return _j1(x) + 1j*_y1(x)


def incident(kvec, points):
    """
    Incident wave

    @param kvec incident wave vector, size 2
    @param points target points, array of size n * 2
    @return complex array of size n
    """
    return jnp.exp(1j*points.dot(kvec))


def gradIncident(kvec, nDotK, points):
    """
    Normal gradient of the incident wave, assumes incident wave is exp(1j * kvec.x)

    @param kvec incident wave vector
    @param nDotK nvec . kvec array of size n
    @param points (source) points, array of size n * 2
    @return complex array of size n
    """
    return 1j*nDotK*incident(kvec, points)


def computeScatteredWave(kvec, xc, yc, point):
    """
    Total scattered wave response at a single observer point, summing up
    contributions from each segment. Written for one point, use jax.vmap to
    apply to many points.

    @param kvec incident wave vector
    @param xc x coordinates representing the contour, must close
    @param yc y coordinates representing the contour, must close
    @param point observer point, size 2
    @return complex value
    """
    kmod = jnp.sqrt(kvec.dot(kvec))

    # contour points, n * 2 array
    pc = jnp.stack([xc, yc], axis=1)

    # vectors from one point to the next, mid-segment points and lengths
    xdot = pc[1:, :] - pc[:-1, :]
    pmid = 0.5*(pc[1:, :] + pc[:-1, :])
    dsdt = jnp.sqrt((xdot*xdot).sum(axis=1))

    # normal vectors at the mid segment locations
    nvec = jnp.stack([-xdot[:, 1]/dsdt, xdot[:, 0]/dsdt], axis=1)

    # vectors from the mid-segment positions to the observer and distances
    rvec = point - pmid
    r = jnp.sqrt((rvec*rvec).sum(axis=1))
    kr = kmod * r

    nDotR = (nvec*rvec).sum(axis=1)
    nDotK = nvec.dot(kvec)

    g = 0.25j * hankel1_0(kr)
    dgdn = -0.25j * kmod * nDotR * hankel1_1(kr) / r

    # shadow side: total wave is nearly zero
    #              => scattered wave amplitude = -incident wave ampl.
    #
    # illuminated side:
    #              => scattered wave amplitude = +incident wave ampl.
    shadow = jnp.where(nDotK > 0., 1., -1.)

    return jnp.sum(dsdt * (-g * gradIncident(kvec, nDotK, pmid) +
                           shadow * dgdn * incident(kvec, pmid)))
