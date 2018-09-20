# scatter

Scatter is a simple program that computes the scattering of an incident wave on a 2D
object. The total wave (incident + scattered) satisfies homomegenous Neunmann boundary
conditions on the obstacle.

The numerical method integrates Green's functions along the boundary and 
requires that the normal derivative of the scattered wave cancels the normal 
derivative of the incident wave there. Using the assumption that the wavelength is much 
smaller than the size of the obstacle, the total wave amplitude is approximately set to zero on
the shadow side. On the illuminated side, the scattered amplitude set to the incident wave 
amplitude as described in Morse and Feshbach, Methods of Theoretical Physics, pp 1551-1552.

The numerical approach involves looping over observer points, nodes on a uniform grid, with the 
scattered wave amplitude computed by summing all the reflection contributions from each boundary 
segment.

This program shows how a Python code can be modified to run faster using various techniques. 

## Requirements

 * Python 2.7.x (tested 2.7.14) or 3.x (tested 3.6.3)
 * Scipy (tested 1.0.0, 1.0.1)
 * Numpy (tested 1.13.1, 1.13.3)
 * Boost (1.6.8, 1.6.1)
 * C++ compiler

## Directory layout

The directories below 
```
|/'
   |-original
   |-vect
   |-cext
```
contain different versions of the scatter code, all producing the same result but each having 
a different implementation. For instance `vect` has a vectorised implementation and `cext` 
has some functions coded in C++ for additional speed. 

## How to run the code

Regardless of the version of the code, use
```
python scatter.py
```
If you change the code and need to verify that the result has not changed, 
```
python scatter.py -checksum
```
as this will print out a check sum. Additional options can be invoked to adjust the problem size
and other parameters (type `python scatter.py --help` for a full list). The result can be saved in a VTK file using
```
python scatter.py -save
``` 
VisIt, Paraview and other tools can then be used to visualise the field.

![Total wave amplitude at time t = 0](https://raw.githubusercontent.com/pletzer/scatter/master/scatter_result.png)








