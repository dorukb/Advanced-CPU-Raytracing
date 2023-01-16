#ifndef __DORKTRACER_ENVLIGH__
#define __DORKTRACER_ENVLIGH__
#include "image.h"
#include <random>

namespace DorkTracer
{
    class SphericalEnvironmentLight
    {

    public:
        int id;
        Image* image;
        SphericalEnvironmentLight(int id, Image* img)
        {
            this->id = id;
            this->image = img;
                    
            this->gen = std::mt19937(rand());
            this->uniformDistro = std::uniform_real_distribution<>(-1.0f, 1.0f);
        }
        Vec3f GetSample(Vec3f& dir)
        {
            // Assumes the image is a spherical environment map.
            float u = (1 + (std::atan2(dir.x, -dir.z) / M_PI)) / 2.0f;
            float v = std::acos(dir.y) / M_PI;

            int w = image->width;
            int h = image->height;

            int i = w * u;
            int j = h * v;
            Vec3f luminance = image->GetSample(i,j) * 2 * M_PI;
            return luminance;
        }

        Vec3f GetDirection(Vec3f& surfaceNormal)
        {
            Vec3f n = makeUnit(surfaceNormal);

            // do Random Rejection Sampling.
            // return a direction in the upper hemisphere.
            int failSafe = 0;

            Vec3f candidate;

            while(true)
            {
                candidate.x = uniformDistro(gen);
                candidate.y = uniformDistro(gen);
                candidate.z = uniformDistro(gen);
                float length = len(candidate);

                if(length <= 1.0f && dot(n, candidate) > 0.0f)
                {
                    candidate / length; // normalize.
                    break;
                }

                if(failSafe++ > 1000) {
                    std::cout<<"[ERROR] Rejection sampling is not terminating." << std::endl;
                }
            }
            return candidate;
        }

    private:    
        std::mt19937 gen;
        std::uniform_real_distribution<> uniformDistro;
    };
}

#endif