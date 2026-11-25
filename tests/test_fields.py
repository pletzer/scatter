"""
Compare the field computed by every variant, point by point, with the
reference `original` implementation. The checksum used in ci/check_all.sh
is a single number and could hide e.g. a sign error on part of the contour.

Run from the repository root:  pytest
"""
import glob
import os
import shutil
import subprocess
import sys

import numpy
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ARGS = ['-nx', '12', '-ny', '10', '-nc', '24', '-save']
NFRAMES = 20
COMPILED = ('cext', 'numba', 'openmp')
MPIRUN = shutil.which('mpirun')


def build(directory):
    if not glob.glob(os.path.join(directory, 'build', '**', '*.so'), recursive=True):
        subprocess.run([sys.executable, 'setup.py', 'build'], cwd=directory,
                       check=True, capture_output=True)


def read_vtk(filename):
    """Return the 2d array stored in a file written by saveVtk.saveData"""
    with open(filename) as f:
        lines = f.read().splitlines()
    nx1, ny1 = (int(s) for s in lines[4].split()[1:3])
    start = next(i for i, l in enumerate(lines) if l.startswith('total ')) + 1
    values = numpy.array(' '.join(lines[start:]).split(), dtype=numpy.float64)
    return values.reshape((ny1, nx1))


def run(directory, command=(), tmp_path=None):
    """Run scatter.py in `directory` with the current directory set to tmp_path,
    return the list of fields (one per time frame)"""
    subprocess.run(list(command) + [sys.executable, os.path.join(directory, 'scatter.py')] + ARGS,
                   cwd=tmp_path, check=True, capture_output=True)
    return [read_vtk(tmp_path / 'scatter_{}.vtk'.format(it)) for it in range(NFRAMES)]


@pytest.fixture(scope='module')
def reference(tmp_path_factory):
    return run(os.path.join(ROOT, 'original'), tmp_path=tmp_path_factory.mktemp('ref'))


def variants():
    for top in ('.', 'solutions'):
        for name in ('original', 'vect', 'numba', 'cext', 'multiproc', 'mpi', 'openmp', 'jax'):
            path = os.path.join(top, name)
            if os.path.isdir(os.path.join(ROOT, path)):
                yield pytest.param(path, id=path.replace('./', '').replace(os.sep, '/'))


@pytest.mark.parametrize('path', variants())
def test_field_matches_original(path, reference, tmp_path):
    directory = os.path.join(ROOT, path)
    command = ()
    if os.path.basename(path) in COMPILED:
        build(directory)
    if os.path.basename(path) == 'mpi':
        if MPIRUN is None:
            pytest.skip('mpirun not available')
        if not path.startswith('solutions'):
            pytest.skip('the mpi exercise does not gather the field yet')
        command = (MPIRUN, '--oversubscribe', '-n', '2')
    fields = run(directory, command, tmp_path)
    for it in range(NFRAMES):
        numpy.testing.assert_allclose(fields[it], reference[it], rtol=1.e-4, atol=1.e-4,
                                      err_msg='time frame {}'.format(it))


@pytest.mark.parametrize('name', COMPILED)
def test_helpful_error_if_not_built(name, tmp_path):
    """Copy a compiled variant without its build directory, from another cwd"""
    src = os.path.join(ROOT, name)
    dst = tmp_path / name
    shutil.copytree(src, dst, ignore=shutil.ignore_patterns('build', '__pycache__'))
    result = subprocess.run([sys.executable, str(dst / 'scatter.py'), '-nx', '4', '-ny', '4'],
                            cwd=tmp_path, capture_output=True, text=True)
    assert result.returncode != 0
    assert 'python setup.py build' in result.stderr
