#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#define STB_IMAGE_IMPLEMENTATION

#define TINYEXR_IMPLEMENTATION

#include <chrono>
#include <iostream>
#include <thread>
#include <random>
#include "raytracer.hpp"
#include "helperMath.h"
#include "gaussian.h"

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

    int samplesPerPixel = cam->samplesPerPixel;
    bool isMultisampling = samplesPerPixel > 1;
    int nRows = std::sqrt(samplesPerPixel);
    int nCols = nRows;

    std::vector<Vec2f> samples(samplesPerPixel);

    std::mt19937 randomGen(rand());
    std::uniform_real_distribution<> uniform01(0.0f, 1.0f);

    float pixelWidth = 1.0f;
    Gaussian2D gaussian2d(pixelWidth / 6.0f);

    for (int y = startingHeight; y < endingHeight; y++)
    {
        for (int x = 0; x < width; x++)
        {
            Vec3f color{0.0f,0.0f,0.0f};
            if(isMultisampling)
            {
                // send samplerPerPixel rays thru each pixel, determine sample locations using:
                // Stratified Random Sampling.

                int i = 0;
                for(int row = 0; row < nRows; row++)
                {
                    for(int col = 0; col < nCols; col++)
                    {
                        float psi1 = uniform01(randomGen);
                        float psi2 = uniform01(randomGen);
                        samples[i].x = (col+psi1) / nCols;
                        samples[i].y = (row+psi2) / nRows;
                        i++;
                    }
                }

                
                // apply 2D gaussian filter.
                float sumOfWeights = 0.0f;
                for(i = 0; i < samplesPerPixel; i++)
                {
                    Vec3f col = renderer->RenderPixel(samples[i].x + x, samples[i].y + y, *cam);

                    float xDistFromCenter = samples[i].x - 0.5f;
                    float yDistFromCenter = samples[i].y - 0.5f;

                    assert(xDistFromCenter > -0.50001f && xDistFromCenter < 0.5001f);
                    assert(yDistFromCenter > -0.50001f && yDistFromCenter < 0.5001f);

                    // get gaussian weight for this one.
                    float gWeight = gaussian2d.GetWeight(xDistFromCenter,yDistFromCenter);
                    color.x += col.x * gWeight;
                    color.y += col.y * gWeight;
                    color.z += col.z * gWeight;
                    sumOfWeights += gWeight;
                }
                color.x = color.x / sumOfWeights;
                color.y = color.y / sumOfWeights;
                color.z = color.z / sumOfWeights;
            }
            else
            {
                color = renderer->RenderPixel(x, y, *cam);
            }

            // // TODO: apply a tonemapping operator instead of simple clamping.
            Vec3i finalColor;
            if(cam->hasTonemapper)
            {
                finalColor = cam->GetTonemappedColor(color);
            }
            else finalColor = clamp(color);

            uint32_t imgIdx = 3 * (x + y * width);
            image[imgIdx] = finalColor.x;
            image[imgIdx+1] = finalColor.y;
            image[imgIdx+2] = finalColor.z;
        }
    }

}

int main(int argc, char* argv[])
{
    // Sample usage for reading an XML scene file
    DorkTracer::Scene scene;
    
    scene.loadFromXml(argv[1]);
    auto start = std::chrono::steady_clock::now();

    DorkTracer::Raytracer renderer(scene);
    
    for(int i = 0; i < scene.cameras.size(); i++)
    {
        DorkTracer::Camera& cam = scene.cameras[i];
        int width = cam.imageWidth, height = cam.imageHeight;
        unsigned char* image = new unsigned char [width * height * 3];

        std::cout<<"Resolution: "<<width <<"x"<<height << std::endl;
        std::cout<<"Running on: "<<THREAD_COUNT <<" threads." << std::endl;

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

        stbi_write_png((cam.imageName + ".png").c_str(), width, height, 3, image, width * 3);
    }

    
    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;
    std::cout << "Rendering took: " << elapsed_seconds.count() << "s\n";
}
