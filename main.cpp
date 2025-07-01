//==========================
// CS 488 - FINAL PROJECT
// (by Benjamin Colussi)
//==========================



// utilities
#include "cs488.h"
inline float clamp(const float x) {
    return x < 0 ? 0 : x > 1 ? 1 : x;
}
inline int toInt(const float x) {
    // return int(pow(clamp(x), 1 / 2.2) * 255 + .5); // gamma correction
    return int(clamp(x) * 255 + .5); // no gamma correction
}



// global variables
static TriangleMesh mesh;
static PointLightSource light;



// main renderer
int main(int argc, const char* argv[]) {

    // .ppm file path
    std::string filename = "no-object.ppm";

    // load .obj on command line
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
    light.position = float3(3.0f, 3.0f, 3.0f);
    light.wattage = float3(1000.0f, 1000.0f, 1000.0f);
    globalScene.addLight(&light);

    // scene calculations
    globalScene.preCalc();

    // ray trace
    globalScene.Raytrace();
    if (MEASURE_RAYS_PER_SECOND) {
        const float rps = RAY_ACCUMULATION / TIME_ACCUMULATION.count();
        std::cout << rps << " rays / second" << std::endl;
        RAY_ACCUMULATION = 0;
        TIME_ACCUMULATION = std::chrono::duration<double>::zero();
    }

    // write to .ppm file
    FILE *f = fopen(filename.c_str(), "w");
    fprintf(f, "P3\n%d %d\n%d\n", globalWidth, globalHeight, 255);
    for (int j = globalHeight - 1; j >= 0; --j) {
        for (int i = 0; i < globalWidth; ++i) {
            fprintf(f, "%d %d %d ", toInt(FrameBuffer.pixel(i, j)[0]), toInt(FrameBuffer.pixel(i, j)[1]), toInt(FrameBuffer.pixel(i, j)[2]));
        }
    }

}
