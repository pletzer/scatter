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

## Requirements

 * Python 2.7.x (tested 2.7.14) or 3.x (tested 3.6.3)
 * Scipy (tested 1.0.0, 1.0.1)
 * Numpy (tested 1.13.1, 1.13.3)

## How to run the code

Type 
```
python scatter.py --help
```
for a list of options. For a wavelength of 0.3, a circular obstacle using 64 segments on a 100 x 200 grid type
```
python scatter.py -nx 100 -ny 200 -nc 64 -xc "cos(2*pi*t)" -yc "sin(2*pi*t)" -lambda 0.3 -save -checksum
```

## Example results

The above command will generate multiple time frames, saved in VTK files and which can be visualized using 
VisIt, Paraview or other tools. Below is the first time frame.

![Total wave amplitude at time t = 0](https://raw.githubusercontent.com/pletzer/scatter/master/scatter_result.png)








