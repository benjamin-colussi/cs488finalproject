# READ ME BRO

## Completed objectives:
1. Refactor? I hardly know her ...
2. Path tracing w/ RR, NEE & MIS :^D
3. Spherical area lights :^O

## Remaining objectives:
4. BxDF models
    (i) specular reflection - done
    (ii) specular refraction w/ Fresnel - do later
    (iii) microfacet model for rough metal - do later
5. Volumetric scattering - ray marching
6. Atmospheric scattering - Beer, Rayleigh, Mie
7. Build scene - do later



## Done:
* set up github repo
* refactored existing renderer from base code
* output single .ppm file
* set up command line loading
* switched from clion to vscode
* got parallelization working with OpenMP
* installed gcc and using compiler optimizations - we r rly flyin now !!!
* i uninstalled anaconda so ill probs have to fix this later ... womp

## Next:
* atmospheric scattering
* specular refraction & fresnel
* filter camera rays
* BRDF models
* SAH-BVH



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
* fix SAH-BVH
* fix image class with my own gamma correction and tone mapping
