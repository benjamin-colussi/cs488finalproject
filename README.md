# READ ME BRO

## Done:
* set up github repo
* added smallpathtracer.cpp and got it working
* dissected smallpathtracer.cpp
* refactored existing renderer from base code
* output single .ppm file
* set up command line loading
* switched from clion to vscode
* got parallelization working with OpenMP
* installed gcc and using compiler optimizations
* we r rly flyin now




## Friday:
* implement Malley for cosine weighted sampling
* fix direct light sampling
* implement mirror reflections
* cosine weighted importance sampling for diffuse surfaces
* cosine weighted importance sampling for light surfaces
* multiple importance sampling
* filter camera rays




## To do:
Work through the LTE
Figure out importance sampling and MIS
Set up BRDF class
Set up the material class w/ proper BSDF, PDF, sampling
Consult paper on drawing better orthonormal bases
Set up geometry class with material, intersection, etc.




## To do:
* implement area lighting as spheres
* fix geometric normal calculations in triangle intersect
* unidirectional path tracing
* set up static scene with atmosphere and moons
* atmospheric scattering
* create test object files in Blender

## Later:
* better OOP - multiple files for classes, separate header/implementation, etc.
* make sure not dividing by 0 in ray tracing - as per feedback on A1
* get rid of jaggies, aliasing in simple shapes, can see lines in shaded surfaces

## Extensions:
* photon mapping



## Questions:
* for path tracing, when we write the LTE as an integral over paths, do we still generate random direction rays? or should we be sampling from all the objects in scene?
* why divide by 2 pi for uniform hemisphere? steradians? help visualizing?