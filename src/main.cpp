#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#include <chrono>
#include <iostream>
#include <thread>
#include "raytracer.hpp"

#define THREAD_COUNT 8

struct RenderThreadArgs{
    DorkTracer::Raytracer* renderer;
    DorkTracer::Camera* cam;
    int threadIndex;
    unsigned char *imagePtr;
};

void renderThreadMain(RenderThreadArgs args)
{
    int threadIndex = args.threadIndex;
    DorkTracer::Camera* cam = args.cam;
    DorkTracer::Raytracer* renderer = args.renderer;
    unsigned char *image = args.imagePtr;
    
    // calculate the starting and ending indices for width & height using these parameters.
    int width = cam->imageWidth, height = cam->imageHeight;

    int startingHeight = threadIndex * (height / THREAD_COUNT);
    int endingHeight = startingHeight + (height / THREAD_COUNT);
    // Iterate over the part of the image plane that belongs to this thread.
    for (int y = startingHeight; y < endingHeight; y++)
    {
        for (int x = 0; x < width; x++)
        {
            Vec3i color = renderer->RenderPixel(x, y, *cam);

            int imgIdx = 3 * (x + y * width);
            image[imgIdx] = color.x;
            image[imgIdx+1] = color.y;
            image[imgIdx+2] = color.z;
        }
    }
}

int main(int argc, char* argv[])
{
    // Sample usage for reading an XML scene file
    DorkTracer::Scene scene;

    scene.loadFromXml(argv[1]);

    DorkTracer::Raytracer renderer(scene);

    auto start = std::chrono::steady_clock::now();

    for(int i = 0; i < scene.cameras.size(); i++)
    {
        DorkTracer::Camera& cam = scene.cameras[i];
        int width = cam.imageWidth, height = cam.imageHeight;
        unsigned char* image = new unsigned char [width * height * 3];

        std::cout<<"Resolution: "<<width <<"x"<<height << std::endl;
        std::vector<std::thread> renderThreads;
        renderThreads.resize(THREAD_COUNT);

        std::vector<RenderThreadArgs> args;
        args.resize(THREAD_COUNT);

        for(int t = 0; t < THREAD_COUNT; t++)
        {
            args[t].cam = &cam;
            args[t].renderer = &renderer;
            args[t].threadIndex = t;
            args[t].imagePtr = image;
            renderThreads[t] = std::thread(renderThreadMain, args[t]);
        }

        for(int i = 0; i < THREAD_COUNT; i++)
        {
            renderThreads[i].join();
        }

        stbi_write_png(cam.imageName.c_str(), width, height, 3, image, width * 3);
    }

    
    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;
    std::cout << "Rendering took: " << elapsed_seconds.count() << "s\n";
}
