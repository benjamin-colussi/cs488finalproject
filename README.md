# CS 488 - Final Project



## Remaining objectives:
1. WRITE REPORT !!!
2. UPDATE README !!!



## Tidying up:
1. Refactor
    (i) Make sure I understand OpenMP and how I'm using it.
    (ii) Clean up command line interface and main function.
    (iii) Add gamma correction and clamping to Image class.

2. Tidy up code.
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



## Completion:
You should also note whether each objective is completed in your README.
If you wrote code for objectives that you did not get completed,
you might request partial credit there.



### Objective 1: Refactor existing renderer
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



### Objective 2: Path tracing
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



### Objective 3: Light sources
Completed.

I have implemented coloured spherical area lights, which is required for my scene.
I started by building the "Sphere" class, which is also used for other non-emitting
objects. Ray-sphere intersection detection is used, and detects when inside or
outside of a sphere. You can specify a material for each sphere object, as with
triangle meshes.



### Objective 4: Materials

### Objective 5: Volumetric scattering

### Objective 6: Atmospheric scattering

### Objective 7: Build scene
Completed.

Used Blender to create a simple hexagonal shape. Added bevelled edges to create a rocky
effect. Copied the object and placed them in a circular pattern. Sphere objects which have
been implemented in my program could then be defined within main.cpp to add to the scene.
There are some comments within main on how to choose which spheres to use for which
.obj files.
