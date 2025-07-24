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
inline int toInt(const float x) { return int(pow(clamp(x), 1 / 2.2) * 255 + 0.5); }

// render
int main(const int argc, const char* argv[]) {

    // parse command line
    std::string filename = "./media/balls.ppm";
    TriangleMesh mesh;
    bool objLoadSucceed = mesh.load(argv[1]);

    // object file
    if (objLoadSucceed) {

        // set up object
        const std::string s = argv[1];
        filename = s.substr(0, s.length() - 4) + ".ppm";
        globalScene.addObject(&mesh);

        // set up lighting
        Material light;
        light.type = LIGHT;
        light.emission = float3(100.0f);
        // Sphere lightBall = Sphere(float3(-0.25f, 0.25f, -0.25f), 0.025f, light); // initial fog testing
        Sphere lightBall = Sphere(float3(0.0f, 0.4f, 0.0f), 0.05f, light); // secondary fog testing
        globalScene.addBall(&lightBall);

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

    // no object
    else {

        // message
        printf("Failed to load .obj file. Making balls instead.\n");

        // light
        Material light;
        light.type = LIGHT;
        light.emission = float3(100.0f);
        Sphere lightBall = Sphere(float3(0.0f, 1.0f, -3.0f), 0.5f, light);
        globalScene.addBall(&lightBall);

        // green
        Material green;
        green.type = LAMBERTIAN;
        green.Ka = float3(0.08700f, 0.23040f, 0.04680f);
        green.Kd = float3(0.11600f, 0.30720f, 0.06240f);
        green.Ks = float3(0.00000f, 0.00000f, 0.00000f);
        green.Ns = 5.00000f;

        // white
        Material white;
        white.type = LAMBERTIAN;
        white.Ka = float3(0.58800f, 0.51060f, 0.32220f);
        white.Kd = float3(0.78400f, 0.68080f, 0.42960f);
        white.Ks = float3(0.00000f, 0.00000f, 0.00000f);
        white.Ns = 5.00000f;

        // red
        Material red;
        red.type = LAMBERTIAN;
        red.Ka = float3(0.40020f, 0.04920f, 0.00240f);
        red.Kd = float3(0.53360f, 0.06560f, 0.00320f);
        red.Ks = float3(0.00000f, 0.00000f, 0.00000f);
        red.Ns = 5.00000f;

        // metal
        Material metal;
        metal.type = METAL;
        metal.Ks = float3(0.64f);

        // balls
        Sphere leftBall = Sphere(float3(-0.5f, -0.5f, -0.75f), 0.5f, green);
        globalScene.addBall(&leftBall);
        Sphere middleBall = Sphere(float3(0.0f, -2.25f, -1.0f), 2.0f, white);
        globalScene.addBall(&middleBall);
        Sphere rightBall = Sphere(float3(0.5f, -0.5f, -1.25f), 0.5f, red);
        globalScene.addBall(&rightBall);

        // balls
        // Material metal;
        // metal.type = METAL;
        // metal.Ks = float3(0.64f);
        // Sphere leftBall = Sphere(float3(-50.0f, -100.0f, -150.0f), 75.0f, metal);
        // globalScene.addBall(&leftBall);
        // Sphere middleBall = Sphere(float3(0.0f, -1045.0f, -175.0f), 1000.0f, metal);
        // globalScene.addBall(&middleBall);
        // Sphere rightBall = Sphere(float3(50.0f, -100.0f, -200.0f), 75.0f, metal);
        // globalScene.addBall(&rightBall);

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

}
