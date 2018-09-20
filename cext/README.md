 To compile the code type

```
BOOST_DIR=<path-to-boost> python setup.py build
```
You might need in addition to set either
```
export LD_LIBRARY_PATH=<path-to-boost>:$LD_LIBRARY_PATH
```
on Linux, or
```
export DYLD_LIBRARY_PATH=<path-to-boost>:$DYLD_LIBRARY_PATH
```
on Mac OS X.
