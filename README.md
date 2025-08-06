# CS 488 - Final Project



## Acknowledgements:
The code in this project has been implemented as an extension of the CS 488/688
base code written by Toshiya Hachisuka.

The book "Physically Based Rendering: From Theory To Implementation" by Matt Pharr,
Wenzel Jakob, and Greg Humphreys was heavily consulted during development.

The program "smallpt.cpp" by Kevin Beason was referenced during development.
Some small amounts of code were used directly.



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
there are boolean "switches" that are important for building and running the program
differently.

The first is "ATMOSPHERIC_SCATTERING" which builds the program with
volume scattering (set == true to turn on).

There are accompanying sphere objects defined in main.cpp to be used with the .obj files.
You may choose a scene to render by setting the corresponding switch == true.

To render cornellbox.obj, cornellbox-metal.obj, cornellbox-glass.obj, set
RENDER_CORNELL_BOX == true.

To render diorama.obj, set RENDER_GLASS_DIORAMA or RENDER_MICROFACET_DIORAMA == true.

To render mirrodin.obj, set RENDER_FINAL_SCENE == true. This switch changes the
path tracing algorithm very slightly. Light sources, when hit directly from a camera ray,
will appear less bright, so as to simulate a glowing orb.
This is not entirely physically realistic.



## Implementation:
I chose to implement an iterative path tracing algorithm, so I am accumulating radiance
and throughput as we iterate. We use OpenMP to parallelize the main for loop iterating
over all pixels. Camera rays are stratified and randomly jittered within their strata.
I did not implement subpixel filtering because it didn't seem necessary during
development and testing, but I intend to add this in the future.

My program only considers geometric normals and not shading normals. Shading could
easily be implemented because I already have interpolated shading normals implemented
for triangle meshes from the assignments.

I implemented SAH_BVH as an extra for assignent 1, but it seems to seg fault for only
one of the test scenes from the base code. I intended to fix this, but it seems to work
with all the scenes I am including in my final project, and increased performance.



## Objectives:

1. Refactor existing renderer:
    * Create an offline rendering process that writes to a .ppm lossless image file.
    * Implement parallelization using OpenMP.
    * Load object files from the command line.

2. Unidirectional (backward) Monte Carlo path tracing:
    * Monte Carlo integration: The Monte Carlo estimator is used to approximate the
    integral in the light transport (rendering) equation.
    * Subpixel sampling using stratified sampling (“jittered” samples):
    Pixels are divided into equally-spaced strata and a randomly-directed ray is
    traced through each to reduce variance and combat aliasing.
    * Russian roulette: Unbiased stopping condition used to improve efficiency
    without increasing variance.
    * Implement the inversion method of sampling from PDFs.
    * Multiple importance sampling: Optimal weighting of importance sampling when
    sampling from multiple functions such as reflection functions, scattering functions,
    incident radiance from light sources, etc.

3. Light source models:
    * Implement area light sources.
    * Implement coloured light sources.
    * Implement spherical objects, ray-sphere intersections, and an efficient sampling
    method for generating a point on the visible surface of the sphere.
    Might be better to implement as circles instead of spheres.

4. Reflection and transmission models:
    * Implement improved BSDF models using microfacet theory, as described in PBR 9.6,
    ie. microfacet distribution, masking/shadowing function, only sampling from
    visible surface, etc.
    * Implement Fresnel equations.

5. Volumetric scattering:
    * Implement ray marching.

6. Atmospheric scattering:
    * Use ray marching to accurately simulate the scattering of light in an atmosphere
    using physical equations such as Beer-Lambert, Rayleigh, Mie, etc.

7. Modelling the scene:
    * Use Blender to create a simple, metal, hexagonal landscape with a central, glass orb.
    * Use Blender to create simple test objects.
    * Position static light sources and atmosphere in the scene.



## Tidying up / future work:
* implement that awesome distance sampling method - novak
* add gamma correction to image class, use the std c++ clamping
* fix all boolean checks to not use "== true" to avoid that one bug again
* improve class structure
* geometry class allows surface sample - geometry class with material, intersection, etc.
* material class allows direction sample
* better OOP - multiple files for classes, separate header/implementation, etc.
* get rid of jaggies, aliasing in simple shapes, can see lines in shaded surfaces
* use shading normals
* subpixel filtering
* fix SAH-BVH to not seg fault on that one scene ...
