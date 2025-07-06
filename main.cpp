///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////   CS 488 - SPRING 2025          ////////////////////////////////////////////////////////////////////////////
//////////   FINAL PROJECT - GLIMMERVOID   ////////////////////////////////////////////////////////////////////////////
//////////   by Benjamin Colussi           ////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// incloods
#include "cs488.h"

// compile using gcc as follows for best performance:
// g++-15 -std=c++17 -O3 -fopenmp main.cpp -o finalproject
// run as follows:
// time ./finalproject ./media/cornellbox.obj 4

/* notes:
$ g++ -std=c++17 -Wall main.cpp -o finalproject
$ type g++
>> g++ is hashed (/usr/bin/g++) // uses clang
$ type g++-15
>> g++-15 is hashed (/usr/local/bin/g++-15) // uses gcc
*/

// utilities
inline float clamp(const float x) { return x < 0 ? 0 : x > 1 ? 1 : x; }
inline int toInt(const float x) { return static_cast<int>(clamp(x) * 255); }



// main renderer
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

    // set up lighting
    // PointLightSource pointLightSource;
    // light.position = float3(3.0f, 3.0f, 3.0f);
    // light.wattage = float3(1000.0f, 1000.0f, 1000.0f);
    SphericalLightSource sphericalLightSource;
    sphericalLightSource.centre = float3(3.0f);
    sphericalLightSource.radius = 0.5f;
    sphericalLightSource.emission = float3(1000.0f);
    globalScene.addSphericalLightSource(&sphericalLightSource);

    // scene calculations
    globalScene.preCalc();

    // ray trace
    // globalScene.rayTrace();

    // path trace
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
