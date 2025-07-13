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

## I believe the following are working now:
* uniform and cosine-weighted hemisphere sampling
* along with PDFs ...
* actually we still have to implement the light PDFs



## Next:
* i kinda wanna try to fix the VSCode conda shell bug again ... just try for an hour or so ...
* read about multiple importance sampling and then implement
* fix all boolean checks to not use "== true" to avoid that one bug again




* i think we need to change NEE code to use area formulation
* we should try to learn about MIS next, because we have sampling methods ready to go

## Currently working on:
* implemented malley for cos weighted importance sampling with the duff othonormal basis
* now let's work on sampling from lights better ...
* this will involve what? ...
* 

* fix direct light sampling
* implement mirror reflections
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
* generating cos weighted hemisphere direction in small path tracer vs pbrt?

* is it possible to sample the same path twice? should we account for this?

* if we sample from a BRDF for a diffuse surface, we r getting a random direction - so we use the solid angle integral?
* if we sample from the surface of a light, we get a point on the surface and test visibility - so we use the area integral?
* can we mix the two approaches? do we have to be careful about anything?
* is this where multiple importance sampling comes in?

* how to model the light itself accurately? we use wattage but im confused ...
* should it fall off at distance increases?

* a little stumped on perfect reflection ...

* what is the reason for passing by reference? is it to not overflow?

* is the method in the lecture slides better than disk method and malley in pbr for drawing cosine weighted samples?

* when drawing "samples" are we referring to drawing "paths" or "next vertices or directions" ???
* because when we are "sampling" to estimate the integral, we are technically drawing samples from surface area or solid angle
* in order to estimate the integrand, which 
* we are sampling these paths, and technically as we sample an additional vertex, we are sampling an entirely new path, which is our method of incrementally building the paths
* im trying to understand how the sampling methods align with the mathematics fo the integral

* what the eff is f