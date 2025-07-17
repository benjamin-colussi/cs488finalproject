# READ ME BRO

## Completed objectives:
1. Refactor
2. Path tracing
3. Spherical area lights

## Remaining objectives:
4. BxDF models
5. Volumetric scattering
6. Atmospheric scattering
7. Build scene



## Done:
* set up github repo
* refactored existing renderer from base code
* output single .ppm file
* set up command line loading
* switched from clion to vscode
* got parallelization working with OpenMP
* installed gcc and using compiler optimizations - we r rly flyin now !!!
* i uninstalled anaconda so ill probs have to fix this later ... womp

## I believe the following are working:
* spherical area lights
* uniform and cosine-weighted disk malley duff hemisphere sampling - along with PDFs
* MIS - might have an issue - not sure if weights are summing to 1 - doesnt rly seem to reduce variance - is this because of only diffuse surfaces?

## Next:
* just check if MIS works better with conic sampling - it seems to work the same? so i guess ive implemented it properly?
* also is my RR cutting at the correct point?
* perfect reflection
* perfect refraction
* filter camera rays using tent or something
* make RR more robust for possible infinite reflections/refractions
* BRDF models



## Tidying up:
fix geometric normal calculations in triangle intersect
fix all boolean checks to not use "== true" to avoid that one bug again
Set up BRDF class
Set up the material class w/ proper BSDF, PDF, sampling
add materials to lights, improve class structure
geometry class allows surface sample
material class allows direction sample
make sphere a class, and then can add material, which can sample the spectrum or power
geometry class with material, intersection, etc.
try TEV for viewing output

* better OOP - multiple files for classes, separate header/implementation, etc.
* get rid of jaggies, aliasing in simple shapes, can see lines in shaded surfaces



## Questions:
* is it possible to sample the same path twice? should we account for this?
* i think my MIS doesnt make sense - weights sum to 1 if the same path is sampled
* weights dont sum to 1 if 


* what is the reason for passing by reference? is it to not overflow?
* multiplying by 1 / something


* if i am sampling a point on the surface of the sphere, am i using surface area formulation or solid angle? can i use the same method to draw solid angle as drawing points on sphere?
should i be using conic sampling to get a direction towards sphere area lights?