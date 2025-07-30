# READ ME BRO



## Completed objectives:
1. Refactor? I hardly know her ...
2. Path tracing w/ RR, NEE & MIS :^D
3. Spherical area lights :^O
4. Materials
    (i) specular reflection
    (ii) specular refraction w/ Fresnel
5. Volumetric scattering - Beer



## Remaining objectives:
4. Materials
    (iii) microfacet model for rough metal and glossy surfaces
6. Atmospheric scattering - Rayleigh, Mie
7. Build scene
8. WRITE THE REPORT !!!



## Tidying up:
1. Refactor
    (i) Make sure I understand OpenMP and how I'm using it.
    (ii) Clean up command line interface.
2. Path tracing
    (i) Make light sampling more robust.
    (ii) Filter camera rays.
3. Light sources
    (i) Experiment with colours.
4. Materials
    (i) I think I may have to change my compuation of the Fresnel term when refracting to be based off the refracted direction and not the reflected direction.
5. Volumetric scattering
    (i) Maybe try that other sampling method to reduce variance.
6. Atmospheric scattering
    (i) Try to implement Rayleigh or Mie scattering to get some realistic atmospheric effects.
7. Model the scene
    (i) Try to build one simple hexagon in Blender.



## Next:
* going to have to fix my light sampling using method of projection
  because if we sample a direction towards light, and we dont have occlusion, then we should definitely have a light hit - BIAS
* atmospheric scattering - implement rayleigh later
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



## Questions:
* what are the units of the BxDF/phase functions?
* MIS w/ fog - how do the coefficients affect output? should they be less than 1
* MIS w/ fog - just confirm my logic is correct
* sampling specular materials - how do the BxDF work with importance sampling?

* sampling solid angle via spherical coordinates - do we have to deal with change of variables - or can we directly change spherical to cartesian / solid angle ???
