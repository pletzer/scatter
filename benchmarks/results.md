Problem: nx=256 ny=256 nc=256 (66049 observer points x 256 segments).  
Machine: Apple M4, 10 cores, Darwin 25.6.0, Python 3.13.16.  
Wall-clock time of the whole run, best of 1.

| variant | workers | time (s) | speedup vs original |
|---|---:|---:|---:|
| original | 1 | 129.82 | 1.0x |
| vect | 1 | 32.18 | 4.0x |
| numba | 1 | 7.25 | 18x |
| cext | 1 | 1.71 | 76x |
| multiproc | 10 | 30.46 | 4.3x |
| mpi | 10 | 28.89 | 4.5x |
| openmp | 10 | 7.73 | 17x |
| jax | auto | 1.99 | 65x |
