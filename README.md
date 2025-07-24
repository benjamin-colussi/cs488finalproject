# READ ME BRO

## Completed objectives:
1. Refactor? I hardly know her ...
2. Path tracing w/ RR, NEE & MIS :^D
3. Spherical area lights :^O
4. Materials
    (i) specular reflection
5. Volumetric scattering - Beer

## Remaining objectives:
4. Materials
    (ii) specular refraction w/ Fresnel
    (iii) microfacet model for rough metal and glossy surfaces
6. Atmospheric scattering - Rayleigh, Mie
7. Build scene



## Next:
* going to have to fix my light sampling using method of projection
  because if we sample a direction towards light, and we dont have occlusion, then we should definitely have a light hit - BIAS
* atmospheric scattering - implement rayleigh later
* specular refraction & fresnel
* filter camera rays - i think im getting some aliasing for some reason ...
* material models
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
* make sure my code doesnt use the small path tracer code ...
