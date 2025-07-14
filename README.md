###### ###### ###### ###### ###### ###### ###### ###### ###### ######
###### READ ME BRO   ###### ###### ###### ###### ###### ###### ######
###### ###### ###### ###### ###### ###### ###### ###### ###### ######



## Done:
* set up github repo
* refactored existing renderer from base code
* output single .ppm file
* set up command line loading
* switched from clion to vscode
* got parallelization working with OpenMP
* installed gcc and using compiler optimizations - we r rly flyin now !!!
* i uninstalled anaconda so ill probs have to fix this later ... womp


## I believe the following are working now:
* spherical area lights
* uniform and cosine-weighted disk malley duff hemisphere sampling - along with PDFs


## Next:
* read about multiple importance sampling and then implement
* perfect reflection and direct light hits
* filter camera rays using tent or something
* make rays constructor normalize direction automatically

## Next next:
fix all boolean checks to not use "== true" to avoid that one bug again
rename spherical light source and just have one set of light sources
Set up BRDF class
Set up the material class w/ proper BSDF, PDF, sampling
add materials to lights, improve class structure
geometry class allows surface sample
material class allows direction sample
make sphere a class, and then can add material, which can sample the spectrum or power
geometry class with material, intersection, etc.



###### ###### ###### ###### ###### ###### ###### ###### ###### ######
###### ###### ###### ###### ###### ###### ###### ###### ###### ######
###### ###### ###### ###### ###### ###### ###### ###### ###### ######



## To do:
* fix geometric normal calculations in triangle intersect
* unidirectional path tracing
* set up static scene with atmosphere and moons
* atmospheric scattering
* create test object files in Blender

## Later:
* better OOP - multiple files for classes, separate header/implementation, etc.
* get rid of jaggies, aliasing in simple shapes, can see lines in shaded surfaces



###### ###### ###### ###### ###### ###### ###### ###### ###### ######
###### ###### ###### ###### ###### ###### ###### ###### ###### ######
###### ###### ###### ###### ###### ###### ###### ###### ###### ######



## Questions:
* is it possible to sample the same path twice? should we account for this?

* if we sample from a BRDF for a diffuse surface, we r getting a random direction - so we use the solid angle integral?
* if we sample from the surface of a light, we get a point on the surface and test visibility - so we use the area integral?
* can we mix the two approaches? do we have to be careful about anything?
* is this where multiple importance sampling comes in?

* how to model the light itself accurately? we use wattage but im confused ...
* should it fall off at distance increases? or is this accounted for by geometry term? then what is L_e?

* what is the reason for passing by reference? is it to not overflow?
* multiplying by 1 / something

* is the method in the lecture slides better than disk method and malley in pbr for drawing cosine weighted samples?

* when drawing "samples" are we referring to drawing "paths" or "next vertices or directions" ? because if we are drawing samples, then we need to divide by the total product of pdfs?
wouldnt this need to be adjusted by some correlation? the paths are no longer independent ...
* we are sampling these paths, and technically as we sample an additional vertex, we are sampling an entirely new path, which is our method of incrementally building the paths
* im trying to understand how the sampling methods align with the mathematics fo the integral

* what is f and L and the BSDF

* how to get realistic lighting in cornell box? always too bright ... calculating strength of light wrong? or maybe has to do with my MIS weighting

* if i am sampling a point on the surface of the sphere, am i using surface area formulation or solid angle? can i use the same method to draw solid angle as drawing points on sphere?
should i be using conic sampling to get a direction towards sphere area lights?

* NEE: should i be keeping track of total paths and dividing by the whole thing at the end for each pixel?



###### ###### ###### ###### ###### ###### ###### ###### ###### ######
###### ###### ###### ###### ###### ###### ###### ###### ###### ######
###### ###### ###### ###### ###### ###### ###### ###### ###### ######