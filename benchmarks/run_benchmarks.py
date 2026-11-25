"""
Time every variant of scatter on the same problem and write a table and a bar plot.

    python benchmarks/run_benchmarks.py [-nx 256 -ny 256 -nc 256] [--repeat 1]

Results (results.csv, results.md, timings.png) are written next to this script.
The times are wall-clock times of the whole `python scatter.py` process, i.e. they
include interpreter start-up, imports and JIT compilation. The solutions are timed
(the exercises are the serial starting points), plus `original` and `jax`.
"""
import argparse
import csv
import glob
import os
import platform
import re
import shutil
import subprocess
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
HERE = os.path.dirname(os.path.abspath(__file__))
NCPU = os.cpu_count()

# name, directory, number of parallel workers (processes/threads), how many repeats cap
VARIANTS = [
    ('original',          'original',            1),
    ('vect',              'solutions/vect',      1),
    ('numba',             'solutions/numba',     1),
    ('cext',              'solutions/cext',      1),
    ('multiproc',         'solutions/multiproc', NCPU),
    ('mpi',               'solutions/mpi',       NCPU),
    ('openmp',            'solutions/openmp',    NCPU),
    ('jax',               'jax',                 None),   # XLA manages its own threads
]
COMPILED = ('numba', 'cext', 'openmp')


def build(directory):
    if not glob.glob(os.path.join(directory, 'build', '**', '*.so'), recursive=True):
        subprocess.run([sys.executable, 'setup.py', 'build'], cwd=directory,
                       check=True, capture_output=True)


def command(name, workers):
    cmd = [sys.executable, 'scatter.py']
    if name == 'mpi':
        mpirun = shutil.which('mpirun')
        cmd = [mpirun, '--oversubscribe', '-n', str(workers)] + cmd
    return cmd


def run_once(name, directory, workers, args):
    env = dict(os.environ)
    # one thread unless the variant is meant to use several
    env['OMP_NUM_THREADS'] = str(workers if name in ('multiproc', 'openmp') else 1)
    env['OMPI_MCA_rmaps_base_oversubscribe'] = '1'
    cmd = command(name, workers) + ['-checksum', '-nx', str(args.nx), '-ny', str(args.ny),
                                    '-nc', str(args.nc)]
    t0 = time.perf_counter()
    res = subprocess.run(cmd, cwd=directory, env=env, capture_output=True, text=True, check=True)
    dt = time.perf_counter() - t0
    checksum = float(re.search(r'amplitudes\|\^2: ([-+0-9.eE]+)', res.stdout).group(1))
    return dt, checksum


def machine():
    brand = platform.processor()
    if sys.platform == 'darwin':
        brand = subprocess.run(['sysctl', '-n', 'machdep.cpu.brand_string'],
                               capture_output=True, text=True).stdout.strip() or brand
    return '{}, {} cores, {} {}, Python {}'.format(brand, NCPU, platform.system(),
                                                 platform.release(), platform.python_version())


def main():
    p = argparse.ArgumentParser()
    p.add_argument('-nx', type=int, default=256)
    p.add_argument('-ny', type=int, default=256)
    p.add_argument('-nc', type=int, default=256)
    p.add_argument('--repeat', type=int, default=1, help='runs per variant, best time is kept')
    p.add_argument('--only', nargs='*', help='only time these variants')
    p.add_argument('--from-csv', action='store_true',
                   help='do not run anything, regenerate results.md and timings.png from results.csv')
    args = p.parse_args()

    rows = []
    if args.from_csv:
        with open(os.path.join(HERE, 'results.csv')) as f:
            for r in csv.DictReader(f):
                rows.append({'variant': r['variant'],
                             'workers': int(r['workers']) if r['workers'].isdigit() else r['workers'],
                             'seconds': float(r['seconds']), 'checksum': float(r['checksum'])})
    for name, rel, workers in ([] if args.from_csv else VARIANTS):
        if args.only and name not in args.only:
            continue
        directory = os.path.join(ROOT, rel)
        if name in COMPILED:
            build(directory)
        runs = [run_once(name, directory, workers, args) for _ in range(args.repeat)]
        best = min(t for t, _ in runs)
        checksum = runs[0][1]
        rows.append({'variant': name, 'workers': workers if workers else 'auto',
                     'seconds': round(best, 2), 'checksum': checksum})
        print('{:10s} {:8.2f} s  checksum {}'.format(name, best, checksum), flush=True)

    ref = next((r for r in rows if r['variant'] == 'original'), rows[0])
    for r in rows:
        r['speedup'] = round(ref['seconds'] / r['seconds'], 2)
        # all variants must compute the same thing (float32 sums: ~1e-4 relative)
        assert abs(r['checksum'] - ref['checksum']) < 1.e-4 * abs(ref['checksum']), r

    with open(os.path.join(HERE, 'results.csv'), 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=['variant', 'workers', 'seconds', 'speedup', 'checksum'])
        w.writeheader()
        w.writerows(rows)

    with open(os.path.join(HERE, 'results.md'), 'w') as f:
        f.write('Problem: nx={} ny={} nc={} ({} observer points x {} segments).  \n'.format(
            args.nx, args.ny, args.nc, (args.nx + 1)*(args.ny + 1), args.nc))
        f.write('Machine: {}.  \n'.format(machine()))
        f.write('Wall-clock time of the whole run, best of {}.\n\n'.format(args.repeat))
        f.write('| variant | workers | time (s) | speedup vs original |\n|---|---:|---:|---:|\n')
        for r in rows:
            f.write('| {variant} | {workers} | {seconds:.2f} | {speedup:.2f}x |\n'.format(**r))

    plot(rows, args)


def plot(rows, args):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    ink, muted, grid, bar = '#1f2328', '#59636e', '#e1e4e8', '#2a6fdb'
    rows = sorted(rows, key=lambda r: -r['seconds'])
    names = ['{} ({})'.format(r['variant'], r['workers']) if r['workers'] not in (1,) else r['variant']
             for r in rows]
    secs = [r['seconds'] for r in rows]

    fig, ax = plt.subplots(figsize=(8, 4.2), dpi=160)
    y = range(len(rows))
    ax.barh(y, secs, height=0.6, color=bar)
    ax.set_yticks(list(y))
    ax.set_yticklabels(names, color=ink, fontsize=10)
    ax.invert_yaxis()
    top = max(secs)
    for yi, r in zip(y, rows):
        ax.text(r['seconds'] + top*0.01, yi, '{:.2f} s  ({:.2f}x)'.format(r['seconds'], r['speedup']),
                va='center', color=ink, fontsize=9)
    ax.set_xlim(0, top*1.22)
    ax.set_xlabel('wall-clock time (s), lower is better', color=muted, fontsize=9)
    ax.xaxis.grid(True, color=grid, linewidth=0.8)
    ax.set_axisbelow(True)
    for s in ('top', 'right', 'left'):
        ax.spines[s].set_visible(False)
    ax.spines['bottom'].set_color(grid)
    ax.tick_params(axis='x', colors=muted, labelsize=9)
    ax.tick_params(axis='y', length=0)
    fig.suptitle('scatter: execution time by implementation', x=0.02, ha='left',
                 color=ink, fontsize=12, fontweight='bold')
    ax.set_title('{}x{} grid, {} segments; {}; speedup vs original in brackets'.format(
        args.nx + 1, args.ny + 1, args.nc, machine().split(',')[0] + ', ' + str(NCPU) + ' cores'),
        loc='left', color=muted, fontsize=8.5, pad=8)
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, 'timings.png'))


if __name__ == '__main__':
    main()
