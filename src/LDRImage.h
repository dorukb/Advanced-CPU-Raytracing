#include <string>
#include "stb_image.h"
#include "image.h"

namespace DorkTracer
{
    class LDRImage : public Image
    {
    public:
        LDRImage(std::string filename, int id)
        {
            this->id = id;
            LoadImage(filename);
        }

        Vec3f GetSample(int i, int j) override
        {
            uint32_t imgIdx = channels * (i + j * width);
            
            Vec3f color;
            color.x = image[imgIdx];
            color.y = image[imgIdx+1];
            color.z = image[imgIdx+2];        

            return color;
        }
        int GetSampleGreyscale(int i, int j) override
        {
            uint32_t imgIdx = i + j * width;
            return image[imgIdx];
        }
        void* getData(){
            return image;
        }
    protected:

        unsigned char* image;   
        void LoadImage(std::string filename) override
        {
            image = stbi_load(filename.c_str(), &width, &height, &channels, 0);
            if(image == NULL) {
                printf("Error in loading the image: %s\n", filename.c_str());
            }
            printf("Loaded image with a width of %dpx, a height of %dpx and %d channels\n", width, height, channels);
        } 

    };
}