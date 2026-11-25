This version uses [JAX](https://jax.readthedocs.io): the whole field is computed
in one jit-compiled call, with `jax.vmap` applying the single-point code over all the observer points
(the segments are vectorised inside `wave.computeScatteredWave`, as in the `vect` version).

```
pip install jax
python scatter.py -checksum
```

Notes
 * JAX has no Bessel functions of the second kind, so `wave.py` implements J0, J1, Y0 and Y1
   with the Abramowitz & Stegun approximations (relative error ~1e-8, ample for the single
   precision results of this code). Double precision (`jax_enable_x64`) is switched on.
 * The grid is processed in a single batch, so memory scales like `nx * ny * nc`; for much
   larger problems process the points in chunks or use `jax.lax.map`.
 * The code runs unchanged on a GPU if a GPU build of JAX is installed.
