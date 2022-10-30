#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#include <chrono>
#include <iostream>
#include "raytracer.hpp"

int main(int argc, char* argv[])
{
    // Sample usage for reading an XML scene file
    parser::Scene scene;

    scene.loadFromXml(argv[1]);

    DorkTracer::Raytracer renderer(scene);

    auto start = std::chrono::steady_clock::now();

    int i = 0;
    // TODO loop over all cameras.
    int camIndex = 0;

    DorkTracer::Camera& cam = scene.cameras[camIndex];
    int width = cam.imageWidth, height = cam.imageHeight;
    unsigned char* image = new unsigned char [width * height * 3];

    for (int y = 0; y < height; ++y)
    {
        for (int x = 0; x < width; x++)
        {
            Vec3i color = renderer.RenderPixel(x, y, camIndex);
            image[i++] = color.x;
            image[i++] = color.y;
            image[i++] = color.z;
        }
    }
    stbi_write_png(cam.imageName.c_str(), width, height, 3, image, width * 3);

    
    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;
    std::cout << "Rendering took: " << elapsed_seconds.count() << "s\n";
}

