///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////   CS 488 - SPRING 2025          ////////////////////////////////////////////////////////////////////////////
//////////   FINAL PROJECT - GLIMMERVOID   ////////////////////////////////////////////////////////////////////////////
//////////   by Benjamin Colussi           ////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

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

// export PATH="/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin"
// getconf PATH

// export PATH=/sbin:/bin:/usr/sbin:/usr/bin >>> got rid of the ld's
// which -a ld

// default no command line argument is 2x2=4
// 4x4=16, 8x8=64, 16x16=256, 32x32=1024, 64x64=4096
// 128x128=16384, 256x256=65536, 512x512=262144, 1024x1024=1048576



// utilities
inline float clamp(const float x) { return x < 0 ? 0 : x > 1 ? 1 : x; }
// inline int toInt(const float x) { return static_cast<int>(clamp(x) * 255); } // no gamma correction
inline int toInt(const float x) { return int(pow(clamp(x), 1 / 2.2) * 255 + .5); } // with gamma correction



// renderer
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
    // set up lighting - camera is at (0, 0, 1.5) looking negative z
    // float3 centre(0.0f, 0.0f, 5.0f); float radius(3.0f); // initial testing
    // float3 centre(-0.5f, 1.0f, 2.0f); // secondary testing

    float3 centre(-0.25f, 0.25f, -0.25f); // in box top left back corner
    // float3 centre(-0.25f, 0.25f, 0.25f); // in box top left front corner
    
    float radius(0.025f);
    // float radius(0.05f);
    // float radius(0.5f);
    
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
