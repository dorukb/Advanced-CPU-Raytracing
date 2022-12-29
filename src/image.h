#ifndef __DORKTRACER_IMAGE__
#define __DORKTRACER_IMAGE__

#include "stb_image.h"
#include <string>
#include "helperMath.h"

namespace DorkTracer
{
    class Image
    {

    public:
        int id;
        int width;
        int height;

        Image(std::string filename, int id)
        {
            this->id = id;
            LoadImage(filename);
        }

        Vec3i GetSample(int i, int j)
        {
            uint32_t imgIdx = channels * (i + j * width);
            
            Vec3i color;
            color.x = image[imgIdx];
            color.y = image[imgIdx+1];
            color.z = image[imgIdx+2];

            return color;
        }
        int GetSampleGreyscale(int i, int j)
        {
            uint32_t imgIdx = i + j * width;
            return image[imgIdx];
        }

    private:
        int channels;
        unsigned char* image;

        void LoadImage(std::string filename)
        {
            image = stbi_load(filename.c_str(), &width, &height, &channels, 0);
            if(image == NULL) {
                printf("Error in loading the image\n");
            }
            printf("Loaded image with a width of %dpx, a height of %dpx and %d channels\n", width, height, channels);
        }

    };
}
#endif