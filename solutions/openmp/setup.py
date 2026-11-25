import os
import sys
from setuptools import setup, Extension

# directory where BOOST is installed
boost_dir = '/opt/homebrew' if os.path.isdir('/opt/homebrew/include/boost') else '/usr/local'
if os.environ.get('BOOST_DIR'):
    boost_dir = os.environ.get('BOOST_DIR')
elif os.environ.get('EBROOTBOOST'):
    boost_dir = os.environ.get('EBROOTBOOST')

# OpenMP flags. Apple clang needs the libomp package (brew install libomp)
include_dirs = [boost_dir + '/include']
compile_args = ['-fopenmp']
link_args = ['-fopenmp']
if sys.platform == 'darwin':
    omp_dir = os.environ.get('OMP_DIR', '/opt/homebrew/opt/libomp')
    include_dirs.append(omp_dir + '/include')
    compile_args = ['-Xpreprocessor', '-fopenmp']
    link_args = ['-L' + omp_dir + '/lib', '-lomp', '-Wl,-rpath,' + omp_dir + '/lib']

setup(
    name='scatter',
    description='A python code that computes the wave scattered by an obstacle',
    license='BSD',
    scripts=['scatter.py'],
    ext_modules=[Extension('wave', ['src/wave.cpp', 'src/is_inside_contour.cpp'],
                                   include_dirs=include_dirs,
                       extra_compile_args=compile_args,
                       extra_link_args=link_args,
                       ),
                     ],
)
