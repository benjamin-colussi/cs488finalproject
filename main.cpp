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





        // red
        Material redLight;
        redLight.type = LIGHT;
        redLight.emission = float3(53.36f, 6.56f, 0.32f);
        redLight.Kd = float3(0.53360f, 0.06560f, 0.00320f);
        Sphere redMoon = Sphere(float3(-5.5f, 5.0f, -13.0f), 1.0f, redLight);
        globalScene.addLight(&redMoon);

        // green
        Material greenLight;
        greenLight.type = LIGHT;
        greenLight.emission = float3(11.6f, 30.72f, 6.24f);
        greenLight.Kd = float3(0.11600f, 0.30720f, 0.06240f);
        Sphere greenMoon = Sphere(float3(-5.5f, 5.0f, -16.0f), 1.0f, greenLight);
        globalScene.addLight(&greenMoon);

        // white
        Material whiteLight;
        whiteLight.type = LIGHT;
        whiteLight.emission = float3(78.4f, 68.08f, 42.96f);
        whiteLight.Kd = float3(0.78400f, 0.68080f, 0.42960f);
        Sphere whiteMoon = Sphere(float3(4.0f, 4.5f, -18.0f), 1.0f, whiteLight);
        globalScene.addLight(&whiteMoon);

        // blue
        Material blueLight;
        blueLight.type = LIGHT;
        blueLight.emission = float3(11.6f, 6.24f, 30.72f);
        blueLight.Kd = float3(0.11600f, 0.06240f, 0.30720f);
        Sphere blueMoon = Sphere(float3(6.0f, 4.0f, -12.0f), 1.0f, blueLight);
        globalScene.addLight(&blueMoon);

        // black
        Material blackLight;
        blackLight.type = MICROFACET;
        blackLight.Kd = float3(0.1f);
        Sphere blackMoon = Sphere(float3(1.5f, 3.5f, -10.0f), 1.0f, blackLight);
        globalScene.addBall(&blackMoon);

        // glass
        Material glass;
        glass.type = GLASS;
        Sphere glassOrb = Sphere(float3(-0.5f, 1.0f, -14.5f), 4.0f, glass);
        globalScene.addBall(&glassOrb);





        /*

        Material light;
        light.type = LIGHT;
        light.emission = float3(100.0);
        Sphere bigLight = Sphere(float3(0.0f, 2.5f, -5.0f), 0.25f, light);
        globalScene.addLight(&bigLight);

        Sphere backTopCorner = Sphere(float3(-0.25f, 0.25f, -0.25f), 0.025f, light);
        globalScene.addLight(&backTopCorner);
        Sphere backBottomCorner = Sphere(float3(0.25f, -0.25f, -0.25f), 0.025f, light);
        globalScene.addLight(&backBottomCorner);
        Sphere centreCeiling = Sphere(float3(0.0f, 0.4f, 0.0f), 0.05f, light);
        globalScene.addLight(&centreCeiling);
        Sphere centreFloor = Sphere(float3(0.0f, -0.4f, 0.0f), 0.025f, light);
        globalScene.addLight(&centreFloor);
        Sphere centreBack = Sphere(float3(0.0f, 0.0f, -0.4f), 0.025f, light);
        globalScene.addLight(&centreBack);
        Sphere leftWall = Sphere(float3(-0.4f, 0.0f, 0.0f), 0.025f, light);
        globalScene.addLight(&leftWall);
        Sphere rightWall = Sphere(float3(0.4f, 0.0f, 0.0f), 0.025f, light);
        globalScene.addLight(&rightWall);

        // specular metal
        Material metal;
        metal.type = METAL;
        metal.Ks = float3(0.64f);
        Sphere metalBall = Sphere(float3(0.15f, 0.175f, -0.1f), 0.15f, metal);
        globalScene.addBall(&metalBall);

        // specular glass
        Material glass;
        glass.type = GLASS;
        Sphere glassBall = Sphere(float3(-0.15f, 0.0f, 0.1f), 0.15f, glass);
        globalScene.addBall(&glassBall);

        // microfacet
        Material microfacet;
        microfacet.type = MICROFACET;
        microfacet.Kd = float3(0.64f); // steel
        microfacet.Kd = float3(0.75f, 0.5f, 0.5f); // copper
        microfacet.Kd = float3(1.022f, 0.782f, 0.344f); // gold
        microfacet.Kd = float3(0.1f, 0.1f, 0.99f); // blue
        Sphere microfacetBall = Sphere(float3(0.15f, -0.175f, 0.0f), 0.15f, microfacet);
        globalScene.addBall(&microfacetBall);

        */



        

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
        Sphere lightBall = Sphere(float3(0.0f, 0.5f, -3.0f), 0.15f, light);
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
        Sphere leftBall = Sphere(float3(-0.75f, -0.75f, -1.25f), 0.5f, green);
        globalScene.addBall(&leftBall);
        Sphere middleBall = Sphere(float3(0.0f, -1.25f, -1.5f), 1.0f, white);
        globalScene.addBall(&middleBall);
        Sphere rightBall = Sphere(float3(0.75f, -0.75f, -1.75f), 0.5f, red);
        globalScene.addBall(&rightBall);

        // moons
        Sphere greenMoon = Sphere(float3(0.27f, 0.25f, -0.5f), 0.25f, green);
        globalScene.addBall(&greenMoon);
        Sphere redMoon = Sphere(float3(-0.27f, 0.25f, -0.25f), 0.25f, red);
        globalScene.addBall(&redMoon);

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
