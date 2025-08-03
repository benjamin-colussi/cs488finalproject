# CS 488 - Final Project

## Remaining objectives:
1. WRITE REPORT !!!
2. WRITE README !!!
3. MAKE VIDEO !!!
4. RENDER FINAL SCENES !!!

## Tidying up:
* clean up the main function and command line interface
* add gamma correction to image class
* use the std c++ clamping
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






## Acknowledgements:
The code in this project has been implemented as an extension of the CS 488/688
base code written by Toshiya Hachisuka.

The book "Physically Based Rendering: From Theory To Implementation" by Matt Pharr,
Wenzel Jakob, and Greg Humphreys was heavily consulted during development.

The program "smallpt.cpp" by Kevin Beason was referenced when making certain implementation
decisions.



## What does my program do?
This project implements a path tracing algorithm to render physically-based global
illumination of 3D-modelled computer geometry. Just specify a wavefront file and sample
count on the command line after building, and it outputs a .ppm lossless image file
depicting the scene from a pre-specified camera angle.



## Compilation:
Must have the following files to build:
* main.cpp
* cs488.h
* linalg.h

Build on the command line using GCC with the following command:
* g++-15 -O3 -fopenmp main.cpp -o finalproject

Run on the command line with the following command:
* ./finalproject ./wavefrontfile.obj spp
* ex. ./finalproject ./media/cornellbox.obj 4

The first argument is the wavefront file you would like to render. You may
include .mtl files to use with your .obj files for materials. The second argument
is square root of the number of samples the program will use in the Monte Carlo
estimation in the implemented path tracing algorithm.

I have included a small number of .obj and accompanying .mtl files I used when
testing my program and to build my final scene in the ./media folder. The program
will output a single .ppm lossless image file in the same directory as the input.
You will need to use a .ppm file viewer to see the program output. I hear "tev" is
a good one built by someone in the industry (https://github.com/Tom94/tev),
but I did not use this during the development of my final project. I used a simple
extension in Visual Studio Code called "PBM/PPM/PGM Viewer for Visual Studio Code"
(https://marketplace.visualstudio.com/items?itemName=ngtystr.ppm-pgm-viewer-for-vscode).



## Specification:
In the first section of cs488.h, where global constants and variables are defined,
there are two boolean "switches" that are important for building and running the program
differently. The first is "ATMOSPHERIC_SCATTERING" which builds the program with
volume scattering (set == true to turn on). The second is "FINAL_SCENE" which changes the
path tracing algorithm very slightly. Light sources, when hit directly from a camera ray,
will appear less bright, so as to simulate a glowing orb (set == true to turn on).



## Implementation:
Describes some software considerations, where appropriate, about:
* Algorithms, data structures, and complexities.
* Caveats, bugs, cautions, assumptions.

SAH-BVH
overall pt alg from camera ray jitters to shading to dividing by spp to writing image
surface vs. volume scattering algorithms
complexity relative to spp ???
describe my software engineering choices i guess ...
mention openmp
mention the code in main and how to build differently with lights





## Objective 1: Refactor existing renderer
Completed.

I built my final project on top of the provided CS 488 base code. I didn't need all of
the functionality, so I stripped it of everything that was not required.

There are six main sections of code within cs488.h:
1. Global constants and variables
2. Image class (stripped down from base code)
3. Material class (where I implemented BRDF models)
4. Triangles (slightly stripped down and modified from base code)
5. BVH implementation (same as base code - includes my SAH-BVH implementation)
6. Sphere and Scene classes and shaders (where the path tracing algorithm is implemented)

The main.cpp file handles command line argument passing, which specifies the .obj file
to load as a triangle mesh, and the number of samples to use in the Monte Carlo estimation.
Sphere objects are also manually added to the scene in the main function. This is where
light sources are added, as well as any spheres of other materials.

I took some inspiration from smallpt.cpp to implement parallelization with OpenMP,
as well as the simple gamma correction and writing to a .ppm file.

The compiler command "#pragma omp parallel for schedule(dynamic, 1) private(shade)"
specifies to parallelize the following for loop on a dynamic schedule,
because each pixel value can be computed independently of the others, but each
will require an amount of work that is unknown at compile time.
The "private" clause specifies that the variable "shade" will have its own copy of
the variable, so that pixel values do not interfere with each other.



## Objective 2: Path tracing
Completed.

I have implemented a unidirectional path tracing algorithm that traces rays from the
camera (or backwards) into the scene in order to construct paths as samples for
Monte Carlo estimation of the light transport equation.

I have implemented BRDF sampling for path continuation (more on this in the materials
section). Sampling methods, along with BRDF evaluation and PDF computation are used
in the estimator of the integral equation.

I have imlemented Russian roulette in order to decrease runtime of my program without
introducing bias. Russian roulette, in the very best case, will not increase variance,
and in the the worst case, increase variance, but the computational savings make up
for this by allowing us to estimate using more samples, which reduces variance even
more than the Russian roulette increases it.

I have implemented next event estimation in order to sample paths that contribute
more radiance to the estimator, which in turn reduces variance and helps our
estimator converge faster.

I have implemented multiple importance sampling of the BRDF and NEE sampling strategies.
The balance heuristic is used in order choose optimal weighting of the PDFs for each of
the two sampling strategies. MIS with the balance heuristic reduces variance
by weighting samples highly if they have higher overlapping values in the PDFs
of all the sampling strategies.



## Objective 3: Light sources
Completed.

I have implemented coloured spherical area lights, which is required for my scene.
I started by building the "Sphere" class, which is also used for other non-emitting
objects. Ray-sphere intersection detection is used, and detects when inside or
outside of a sphere. You can specify a material for each sphere object, as with
triangle meshes.



## Objective 4: Materials
Completed.

I first implemented specular reflection, which was relatively simple. By calculating
the reflected direction and continuing the path in that direction, we are perfectly
importance sampling the Dirac delta distribution and we simply continue the path without
next event estimation.

I then implemented Fresnel reflection/refraction for glass using Schlick's approximation,
because it seems like the best choice when concerned with performance. We importance
sample the double Dirac delta distribution by choosing whether to reflect or refract
based on the calculated Fresnel reflectance term which is always between 0 and 1, so we
can simply sample a uniformly distributed random variable and continue our path by
calculating the reflected or refracted direction.

Lastly, I implemented a microfacet model 



## Objective 5: Volumetric scattering
Completed.

I initially intended to implement a ray marching algorithm to simulate single scattering,
but I instead chose to implement free-flight sampling on top of my unidirectional
path tracing algorithm. This method is much more robust and much less biased,
because it estimates the actual volume rendering equation by sampling distance
as well as solid angle. This allows us to estimate multiple scattering and not just
single scattering, by allowing light to travel in all directions and eventually
end uop along a ray and scatters towards the camera.

We first sample a direction, and then sample free-flight distance along this direction.
We then check for a ray-surface interaction, and check if the sampled distance is less
than the distance of the surface interaction. If the sampled distance is less,
we have a volume interaction. If the sampled distance is greater, we have a surface
interaction. In both cases we perform next event estimation, and then sample solid angle
to continue the path. In the case of a volume interaction, we sample the phase function
and compute the accumulated radiance based on the volume rendering equation. In the case
of a surface interaction, we accumulate radiance in the usual unidirectional path tracing
manner.

We also implement multiple importance sampling of the sampling methods used in volume
scattering. The only thing that changes from the usual path tracing MIS is that we are now
sampling distance at each step to organically continue the path. This means that our
BRDF/phase function sampling method must now be multiplied by the probability of sampling
the distance at each step. We weight the estimate of the accumulated radiance using the
balance heuristic as before, but now our two sampling strategies are BRDF/phase function
multiplied by distance sampling, and sampling solid angle towards a light source during NEE.

This method allows us to Monte Carlo estimate the double integral of the volume rendering
equation, which involves and integral over distance of an integral over solid angle.
When we have a volume interaction, we multiply by the probability of reaching the
sampled distance. When we have a surface interaction, we multiply by the probability
of reaching or exceeding the distance to the surface, because that is the probability
of sampling a distance at or beyond the surface, if the surface wasn't there.



## Objective 6: Atmospheric scattering
Completed.

I initially intended to implement various atmospheric scattering models but I decided
this was not fully necessary for the scene I was trying to create, considering it is set
in a fictional world, so the atmospheric effects could be much different than Earth’s.
I also ran out of time ...

In the end, I implemented the Beer-Lambert law to simulate the attenuation of light
through a homogeneous atmosphere, in order to achieve a hazy atmospheric effect,
and create a glow around the coloured light sources in my scene.

The Beer-Lambert law states that the radiance of light decays exponentially as it travels
through the medium according to the absorption and scattering coefficients of the medium.
The absorption coefficient is the probability density that the light is absorbed per unit
distance travelled in the medium. Similarly, the scattering coefficient is the probability
density that there is a scattering event per unit distance travelled in the medium.

We use free-flight sampling to Monte Carlo estimate the volume rendering equation.
The Beer-Lambert law describes the survival of radiance, or the probability of a light
ray not interacting up to a certain distance.

We must sample from the exponential distribution in order to importance sample the
transmittance term that arises from the Beer-Lambert law. I have implemented this
sampling by using the inversion method of the CDF to transform our uniform random
numbers into exponential distribution samples.

I am using the simplest phase function (1 / 4 * PI), so I sample a direction on the
unit sphere, which can also be accomplished by sampling spherical coordinates using
the inversion method of the CDF. This allows us to importance sample the phase function.

There is a lot of cancellation in the volume rendering equation when using the
Beer-Lambert law for volume scattering in homogeneous media. When volume scattering,
The probability of sampling a free-flight distance cancels with the transmittance term
that arises from the Beer-Lambert law, leaving only (1 / attenuation coefficient).
We also importance sample the phase function exactly by sampling the unit sphere,
so this cancels perfectly. All that is left is the scattering coefficient from the
in-scattering term. So we only multiply by the scattering coefficient over the
attenuation coefficient.

In the case of surface scattering, the probability of reaching the surface or exceeding
it is (1 - CDF(s)) which is equivalent to the transmittance of the previous vertex to
the surface, so it cancel perfectly and we multiply throughput by 1.



## Objective 7: Build scene
Completed.

Used Blender to create a simple hexagonal shape. Added bevelled edges to create a rocky
effect. Copied the object and placed them in a circular pattern. Sphere objects which have
been implemented in my program could then be defined within main.cpp to add to the scene.
There are some comments within main on how to choose which spheres to use for which
.obj files.





## Revised objectives:

1. Refactor existing renderer:
    * Create an offline rendering process that writes to a .ppm lossless image file.
    * Implement parallelization using OpenMP.
    * Load object files from the command line.

2. Unidirectional (backward) Monte Carlo path tracing:
    * Monte Carlo integration: The Monte Carlo estimator is used to approximate the integral in the light transport (rendering) equation.
    * Subpixel sampling using stratified sampling (“jittered” samples): Pixels are divided into equally-spaced strata and a randomly-directed ray is traced through each to reduce variance and combat aliasing.
    * Russian roulette: Unbiased stopping condition used to improve efficiency without increasing variance.
    * Implement the inversion method of sampling from PDFs.
    * Multiple importance sampling: Optimal weighting of importance sampling when sampling from multiple functions such as reflection functions, scattering functions, incident radiance from light sources, etc.

3. Light source models:
    * Implement area light sources.
    * Implement coloured light sources.
    * Implement spherical objects, ray-sphere intersections, and an efficient sampling method for generating a point on the visible surface of the sphere. Might be better to implement as circles instead of spheres.

4. Reflection and transmission models:
    * Implement improved BSDF models using microfacet theory, as described in PBR 9.6, ie. microfacet distribution, masking/shadowing function, only sampling from visible surface, etc.
    * Implement Fresnel equations.

5. Volumetric scattering:
    * Implement ray marching.

6. Atmospheric scattering:
    * Use ray marching to accurately simulate the scattering of light in an atmosphere using physical equations such as Beer-Lambert, Rayleigh, Mie, etc.

7. Modelling the scene:
    * Use Blender to create a simple, metal, hexagonal landscape with a central, glass orb.
    * Use Blender to create simple test objects.
    * Position static light sources and atmosphere in the scene.
