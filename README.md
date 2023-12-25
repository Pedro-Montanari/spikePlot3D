# spikePlot3D
Generate spike-like 3D graphs in a structured grid from sample
points in an unstructured grid

This function takes sample points (defined by "x" and "y" coordinates
and a value "v") on an unstructured cartesian domain, distributes the 
value "v" of each sample point to a polar domain of radius "r0" 
according to a decaying exponential function, and sums the 
contributions of all polar domains to a structured cartesian domain.

The exponential function used in this algorithm is given by 
$$V = v*\left(1- \frac{r}{r_0}\right)\left(\frac{\exp(r_0-r)}{\exp(r_0)}\right)$$,
where "V" is the calculated value in the current point, "v" is the 
function value at the sample point, and "r" is the distance from the 
sample point to the current point. At the sample point the value of the
function is "v", and at the radius "r0" the value of the function is
zero.

Please read the license carefully and see the example file for a basic usage 
guide and the function file for documentation. 

The code uses the function "linmap" inspired by machinelearning1's post found on
https://machinelearning1.wordpress.com/2014/07/13/linear-vector-mapping-scaling-matlab/
