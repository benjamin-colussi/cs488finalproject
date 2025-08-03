///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////                                                                              ////////////////////////////////
/////////   CS 488 - SPRING 2025                                                       ////////////////////////////////
/////////   FINAL PROJECT - MIRRODIN                                                   ////////////////////////////////
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

        // light material
        Material light;
        light.type = LIGHT;
        light.emission = float3(100.0);

        // glass material
        Material glass; glass.type = GLASS;

        // metal material
        Material metal; metal.type = METAL; metal.Ks = float3(0.64f);

        // microfacet material
        const float3 microfacetColour(0.64f); // grey
        // const float3 microfacetColour(0.955f, 0.697f, 0.652f); // copper
        // const float3 microfacetColour(0.972f, 0.960f, 0.915f); // silver
        // const float3 microfacetColour(1.022f, 0.782f, 0.344f); // gold
        Material microfacet; microfacet.type = MICROFACET; microfacet.Kd = microfacetColour; microfacet.setRoughness(0.5f);
        Material microfacet1; microfacet1.type = MICROFACET; microfacet1.Kd = microfacetColour; microfacet1.setRoughness(0.0f);
        Material microfacet2; microfacet2.type = MICROFACET; microfacet2.Kd = microfacetColour; microfacet2.setRoughness(0.05f);
        Material microfacet3; microfacet3.type = MICROFACET; microfacet3.Kd = microfacetColour; microfacet3.setRoughness(0.15f);
        Material microfacet4; microfacet4.type = MICROFACET; microfacet4.Kd = microfacetColour; microfacet4.setRoughness(0.25f);
        Material microfacet5; microfacet5.type = MICROFACET; microfacet5.Kd = microfacetColour; microfacet5.setRoughness(0.35f);
        Material microfacet6; microfacet6.type = MICROFACET; microfacet6.Kd = microfacetColour; microfacet6.setRoughness(0.45f);
        Material microfacet7; microfacet7.type = MICROFACET; microfacet7.Kd = microfacetColour; microfacet7.setRoughness(0.55f);
        Material microfacet8; microfacet8.type = MICROFACET; microfacet8.Kd = microfacetColour; microfacet8.setRoughness(0.65f);
        Material microfacet9; microfacet9.type = MICROFACET; microfacet9.Kd = microfacetColour; microfacet9.setRoughness(0.75f);
        Material microfacet10; microfacet10.type = MICROFACET; microfacet10.Kd = microfacetColour; microfacet10.setRoughness(0.85f);
        Material microfacet11; microfacet11.type = MICROFACET; microfacet11.Kd = microfacetColour; microfacet11.setRoughness(0.95f);
        Material microfacet12; microfacet12.type = MICROFACET; microfacet12.Kd = microfacetColour; microfacet12.setRoughness(1.0f);

        

        // red
        Material redLight;
        redLight.type = LIGHT;
        redLight.emission = float3(533.6f, 65.6f, 3.2f) / 2;
        redLight.Kd = float3(0.53360f, 0.06560f, 0.00320f);
        // Sphere redMoon = Sphere(float3(-5.5f, 5.0f, -13.0f), 1.0f, redLight);
        // Sphere redMoon = Sphere(float3(-5.5f, 4.0f, -11.0f), 1.5f, redLight);
        Sphere redMoon = Sphere(float3(-8.0f, 7.0f, -30.0f), 5.0f, redLight);
        globalScene.addLight(&redMoon);

        // green
        Material greenLight;
        greenLight.type = LIGHT;
        greenLight.emission = float3(116.0f, 307.2f, 62.4f) / 2;
        greenLight.Kd = float3(0.11600f, 0.30720f, 0.06240f);
        // Sphere greenMoon = Sphere(float3(-5.5f, 5.0f, -16.0f), 1.0f, greenLight);
        // Sphere greenMoon = Sphere(float3(-5.0f, 5.0f, -16.0f), 1.5f, greenLight);
        Sphere greenMoon = Sphere(float3(3.0f, 10.0f, -35.0f), 5.0f, greenLight);
        globalScene.addLight(&greenMoon);

        // white
        Material whiteLight;
        whiteLight.type = LIGHT;
        whiteLight.emission = float3(884.0f, 780.8f, 29.6f) / 2;
        whiteLight.Kd = float3(0.88400f, 0.78080f, 0.02960f);
        // Sphere whiteMoon = Sphere(float3(4.0f, 4.5f, -18.0f), 1.0f, whiteLight);
        // Sphere whiteMoon = Sphere(float3(4.0f, 4.5f, -18.0f), 1.5f, whiteLight);
        // Sphere whiteMoon = Sphere(float3(5.0f, 4.5f, -30.0f), 5.0f, whiteLight);
        // globalScene.addLight(&whiteMoon);

        // blue
        Material blueLight;
        blueLight.type = LIGHT;
        blueLight.emission = float3(116.0f, 62.4f, 307.2f) / 2;
        blueLight.Kd = float3(0.11600f, 0.06240f, 0.30720f);
        // Sphere blueMoon = Sphere(float3(6.0f, 4.0f, -12.0f), 1.0f, blueLight);
        // Sphere blueMoon = Sphere(float3(7.0f, 4.0f, -14.0f), 1.5f, blueLight);
        Sphere blueMoon = Sphere(float3(10.0f, 4.0f, -30.0f), 5.0f, blueLight);
        globalScene.addLight(&blueMoon);

        // black
        Material blackLight;
        blackLight.type = MICROFACET;
        blackLight.Kd = float3(0.1f);
        blackLight.setRoughness(0.5f);
        // Sphere blackMoon = Sphere(float3(1.5f, 3.5f, -10.0f), 1.0f, blackLight);
        // Sphere blackMoon = Sphere(float3(1.0f, 3.5f, -10.0f), 1.5f, blackLight);
        // Sphere blackMoon = Sphere(float3(2.5f, 3.5f, -10.0f), 1.0f, blackLight);
        // globalScene.addBall(&blackMoon);

        // orb
        Sphere orb = Sphere(float3(-0.5f, 1.0f, -14.5f), 4.0f, glass);
        globalScene.addBall(&orb);



        /*
        
        // cornell box
        // Sphere backTopCorner = Sphere(float3(-0.25f, 0.25f, -0.25f), 0.025f, light); globalScene.addLight(&backTopCorner);
        // Sphere centreCeiling = Sphere(float3(0.0f, 0.4f, 0.0f), 0.05f, light); globalScene.addLight(&centreCeiling); // for fog

        // diorama lights
        // Sphere dioramaCeiling = Sphere(float3(0.0f, 0.4f, 0.0f), 0.05f, light); globalScene.addLight(&dioramaCeiling);
        // Sphere dioramaFloor = Sphere(float3(0.0f, -0.4f, 0.0f), 0.05f, light); globalScene.addLight(&dioramaFloor);
        // Sphere dioramaBackWall = Sphere(float3(0.0f, 0.0f, -0.4f), 0.05f, light); globalScene.addLight(&dioramaBackWall);
        // Sphere dioramaLeftWall = Sphere(float3(-0.55f, 0.0f, 0.0f), 0.05f, light); globalScene.addLight(&dioramaLeftWall);
        // Sphere dioramaRightWall = Sphere(float3(0.55f, 0.0f, 0.0f), 0.05f, light); globalScene.addLight(&dioramaRightWall);
        // Sphere dioramaLeftBehind = Sphere(float3(-0.5f, 0.0f, 1.0f), 0.075f, light); globalScene.addLight(&dioramaLeftBehind);
        // Sphere dioramaRightBehind = Sphere(float3(0.5f, 0.0f, 1.0f), 0.075f, light); globalScene.addLight(&dioramaRightBehind);

        // glass diorama
        // Sphere glassBallBig1l = Sphere(float3(-0.3f, 0.2f, -0.1), 0.18f, glass); globalScene.addBall(&glassBallBig1l);
        // Sphere glassBallBig2l = Sphere(float3(-0.3f, -0.2f, -0.1), 0.18f, glass); globalScene.addBall(&glassBallBig2l);
        // Sphere glassBallBig1r = Sphere(float3(0.3f, 0.2f, 0.1), 0.18f, glass); globalScene.addBall(&glassBallBig1r);
        // Sphere glassBallBig2r = Sphere(float3(0.3f, -0.2f, 0.1), 0.18f, glass); globalScene.addBall(&glassBallBig2r);
        // Sphere glassBall1r = Sphere(float3(0.0f, 0.3f, -0.2), 0.09f, glass); globalScene.addBall(&glassBall1r);
        // Sphere glassBall2r = Sphere(float3(0.2f, 0.3f, -0.2), 0.09f, glass); globalScene.addBall(&glassBall2r);
        // Sphere glassBall3r = Sphere(float3(0.4f, 0.3f, -0.2), 0.09f, glass); globalScene.addBall(&glassBall3r);
        // Sphere glassBall4r = Sphere(float3(0.0f, 0.1f, -0.2), 0.09f, glass); globalScene.addBall(&glassBall4r);
        // Sphere glassBall5r = Sphere(float3(0.2f, 0.1f, -0.2), 0.09f, glass); globalScene.addBall(&glassBall5r);
        // Sphere glassBall6r = Sphere(float3(0.4f, 0.1f, -0.2), 0.09f, glass); globalScene.addBall(&glassBall6r);
        // Sphere glassBall7r = Sphere(float3(0.0f, -0.1f, -0.2), 0.09f, glass); globalScene.addBall(&glassBall7r);
        // Sphere glassBall8r = Sphere(float3(0.2f, -0.1f, -0.2), 0.09f, glass); globalScene.addBall(&glassBall8r);
        // Sphere glassBall9r = Sphere(float3(0.4f, -0.1f, -0.2), 0.09f, glass); globalScene.addBall(&glassBall9r);
        // Sphere glassBall10r = Sphere(float3(0.0f, -0.3f, -0.2), 0.09f, glass); globalScene.addBall(&glassBall10r);
        // Sphere glassBall11r = Sphere(float3(0.2f, -0.3f, -0.2), 0.09f, glass); globalScene.addBall(&glassBall11r);
        // Sphere glassBall12r = Sphere(float3(0.4f, -0.3f, -0.2), 0.09f, glass); globalScene.addBall(&glassBall12r);
        // Sphere glassBall1l = Sphere(float3(0.0f, 0.3f, 0.2), 0.09f, glass); globalScene.addBall(&glassBall1l);
        // Sphere glassBall2l = Sphere(float3(-0.2f, 0.3f, 0.2), 0.09f, glass); globalScene.addBall(&glassBall2l);
        // Sphere glassBall3l = Sphere(float3(-0.4f, 0.3f, 0.2), 0.09f, glass); globalScene.addBall(&glassBall3l);
        // Sphere glassBall4l = Sphere(float3(0.0f, 0.1f, 0.2), 0.09f, glass); globalScene.addBall(&glassBall4l);
        // Sphere glassBall5l = Sphere(float3(-0.2f, 0.1f, 0.2), 0.09f, glass); globalScene.addBall(&glassBall5l);
        // Sphere glassBall6l = Sphere(float3(-0.4f, 0.1f, 0.2), 0.09f, glass); globalScene.addBall(&glassBall6l);
        // Sphere glassBall7l = Sphere(float3(0.0f, -0.1f, 0.2), 0.09f, glass); globalScene.addBall(&glassBall7l);
        // Sphere glassBall8l = Sphere(float3(-0.2f, -0.1f, 0.2), 0.09f, glass); globalScene.addBall(&glassBall8l);
        // Sphere glassBall9l = Sphere(float3(-0.4f, -0.1f, 0.2), 0.09f, glass); globalScene.addBall(&glassBall9l);
        // Sphere glassBall10l = Sphere(float3(0.0f, -0.3f, 0.2), 0.09f, glass); globalScene.addBall(&glassBall10l);
        // Sphere glassBall11l = Sphere(float3(-0.2f, -0.3f, 0.2), 0.09f, glass); globalScene.addBall(&glassBall11l);
        // Sphere glassBall12l = Sphere(float3(-0.4f, -0.3f, 0.2), 0.09f, glass); globalScene.addBall(&glassBall12l);

        // microfacet diorama
        // Sphere metalBall = Sphere(float3(-0.3f, 0.2f, -0.1), 0.18f, metal); globalScene.addBall(&metalBall);
        // Sphere microfacetBall = Sphere(float3(-0.3f, -0.2f, -0.1), 0.18f, microfacet); globalScene.addBall(&microfacetBall);
        // Sphere microfacetBall1 = Sphere(float3(0.0f, 0.3f, -0.1), 0.09f, microfacet1); globalScene.addBall(&microfacetBall1);
        // Sphere microfacetBall2 = Sphere(float3(0.2f, 0.3f, -0.1), 0.09f, microfacet2); globalScene.addBall(&microfacetBall2);
        // Sphere microfacetBall3 = Sphere(float3(0.4f, 0.3f, -0.1), 0.09f, microfacet3); globalScene.addBall(&microfacetBall3);
        // Sphere microfacetBall4 = Sphere(float3(0.0f, 0.1f, -0.1), 0.09f, microfacet4); globalScene.addBall(&microfacetBall4);
        // Sphere microfacetBall5 = Sphere(float3(0.2f, 0.1f, -0.1), 0.09f, microfacet5); globalScene.addBall(&microfacetBall5);
        // Sphere microfacetBall6 = Sphere(float3(0.4f, 0.1f, -0.1), 0.09f, microfacet6); globalScene.addBall(&microfacetBall6);
        // Sphere microfacetBall7 = Sphere(float3(0.0f, -0.1f, -0.1), 0.09f, microfacet7); globalScene.addBall(&microfacetBall7);
        // Sphere microfacetBall8 = Sphere(float3(0.2f, -0.1f, -0.1), 0.09f, microfacet8); globalScene.addBall(&microfacetBall8);
        // Sphere microfacetBall9 = Sphere(float3(0.4f, -0.1f, -0.1), 0.09f, microfacet9); globalScene.addBall(&microfacetBall9);
        // Sphere microfacetBall10 = Sphere(float3(0.0f, -0.3f, -0.1), 0.09f, microfacet10); globalScene.addBall(&microfacetBall10);
        // Sphere microfacetBall11 = Sphere(float3(0.2f, -0.3f, -0.1), 0.09f, microfacet11); globalScene.addBall(&microfacetBall11);
        // Sphere microfacetBall12 = Sphere(float3(0.4f, -0.3f, -0.1), 0.09f, microfacet12); globalScene.addBall(&microfacetBall12);

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
