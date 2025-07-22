///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////                                                                              ////////////////////////////////
/////////   CS 488 - SPRING 2025                                                       ////////////////////////////////
/////////   FINAL PROJECT - GLIMMERVOID                                                ////////////////////////////////
/////////   by Benjamin Colussi                                                        ////////////////////////////////
/////////                                                                              ////////////////////////////////
/////////   compile:  g++-15 -std=c++14 -Wall -O3 -fopenmp main.cpp -o finalproject    ////////////////////////////////
/////////   run:      time ./finalproject ./media/cornellbox.obj 4                     ////////////////////////////////
/////////                                                                              ////////////////////////////////
/////////   16x16=256, 32x32=1024, 64x64=4096, 128x128=16384                           ////////////////////////////////
/////////   256x256=65536, 512x512=262144, 1024x1024=1048576                           ////////////////////////////////
/////////                                                                              ////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "cs488.h"

inline float clamp(const float x) { return x < 0 ? 0 : x > 1 ? 1 : x; }
inline int toInt(const float x) { return int(pow(clamp(x), 1 / 2.2) * 255 + .5); }



// render
int main(const int argc, const char* argv[]) {

    // .ppm file path
    std::string filename = "no-object.ppm";

    // load .obj on command line
    TriangleMesh mesh;
    if (argc > 1) {
        bool objLoadSucceed = mesh.load(argv[1]);
        if (!objLoadSucceed) {
            printf("Invalid .obj file.\n");
            printf("Making a single triangle instead.\n");
            mesh.createSingleTriangle();
        }
        else {
            const std::string s = argv[1];
            filename = s.substr(0, s.length() - 4) + ".ppm";
        }
    }
    else {
        printf("Specify .obj file in the command line arguments. Example: CS488.exe cornellbox.obj\n");
        printf("Making a single triangle instead.\n");
        mesh.createSingleTriangle();
    }
    globalScene.addObject(&mesh);


    
    ///////////////////////////////////////////////////////////////////////////////////////////////
    // lighting - camera is at (0, 0, 1.5) - looking towards -z
    float3 centre(-0.25f, 0.25f, -0.25f);
    float radius(0.025f);
    Material material;
    material.type = LIGHT;
    material.emission = float3(100.0f);
    Sphere light(centre, radius, material);
    globalScene.addLight(&light);
    ///////////////////////////////////////////////////////////////////////////////////////////////



    // path trace
    globalScene.preCalc();
    globalScene.pathTrace(argc == 3 ? std::stoi(argv[2]) : 2);

    // write to .ppm file
    FILE *f = fopen(filename.c_str(), "w");
    fprintf(f, "P3\n%d %d\n%d\n", globalWidth, globalHeight, 255);
    for (int j = globalHeight - 1; j >= 0; --j) {
        for (int i = 0; i < globalWidth; ++i) {
            fprintf(f, "%d %d %d ", toInt(globalImage.pixel(i, j)[0]), toInt(globalImage.pixel(i, j)[1]), toInt(globalImage.pixel(i, j)[2]));
        }
    }

}
