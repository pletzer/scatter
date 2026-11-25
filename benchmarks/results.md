Problem: nx=256 ny=256 nc=256 (66049 observer points x 256 segments).  
Machine: Apple M4, 10 cores, Darwin 25.6.0, Python 3.13.16.  
Wall-clock time of the whole run: best of 3 for runs under 20 s, a single run otherwise.

| variant | workers | time (s) | speedup vs original |
|---|---:|---:|---:|
| original | 1 | 125.06 | 1.0x |
| vect | 1 | 30.98 | 4.0x |
| numba | 1 | 6.74 | 19x |
| cext | 1 | 1.65 | 76x |
| multiproc | 10 | 25.19 | 5.0x |
| mpi | 10 | 28.52 | 4.4x |
| openmp | 10 | 0.74 | 169x |
| jax | auto | 1.70 | 74x |
