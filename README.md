# READ ME BRO

## Completed objectives:
1. Refactor? I hardly know her ...
2. Path tracing w/ RR, NEE & MIS :^D
3. Spherical area lights :^O
4. Materials
    (i) specular reflection
    (ii) specular refraction w/ Fresnel
    (iii) microfacet model
5. Volumetric scattering - Beer - multiple scattering, free flight sampling, NEE, MIS
7. Built the scene

## Remaining objectives:
6. Atmospheric scattering - Rayleigh, Mie
8. WRITE THE REPORT !!!
9. Model some test scenes and images.

## Tidying up:
1. Refactor
    (i) Make sure I understand OpenMP and how I'm using it.
    (ii) Clean up command line interface and main function.
    (iii) Add gamma correction and clamping to Image class.
    (iv) Tidy up code.
2. Path Tracing
    (i) Filter camera rays using tent to combat aliasing.
6. Atmospheric scattering
    (i) try to play around with rayleigh scattering by changing rgb values of scattering coefficients
    (ii) implement a better phase function with more forward scattering
    (iii) implement sampling from the phase function





## Tidying up:
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
* make sure my code doesnt use the small path tracer code ...
