#ifndef __DORKTRACER_HDRIMAGE__
#define __DORKTRACER_HDRIMAGE__
#include <string>
#include <iostream>
#include <vector>

#include "image.h"

namespace DorkTracer
{
    class HDRImage : public Image
    {

    public:
        std::vector<float> src;
        // EXRImage exrImage;

        HDRImage(std::string filename, int id)
        {
            this->id = id;
            LoadImage(filename);
        }

        Vec3f GetSample(int i, int j) override
        {
            uint32_t imgIdx = 3 * (i + j * width);
            
            Vec3f color;
            color.x = src[imgIdx];
            color.y = src[imgIdx+1];
            color.z = src[imgIdx+2];
            return color;
        }
       
        int GetSampleGreyscale(int i, int j) override
        {
            uint32_t imgIdx = i + j * width;
            return src[imgIdx];
        }



    protected:

        void LoadImage(std::string filename) override
        {  
            this->channels = 3;

            float *rgba = nullptr;
            const char *err = nullptr;
            int ret = LoadEXR(&rgba, &width, &height, filename.c_str(), &err);
            if (TINYEXR_SUCCESS != ret) {
                std::cerr << "Failed to load EXR file [" << filename << "] code = " << ret << std::endl;
                if (err) {
                    std::cerr << err << std::endl;
                    FreeEXRErrorMessage(err);
                }
            }

            std::cout << "loaded EXR. width x height = " << width << "x" << height << std::endl;
            src.resize(size_t(width * height * 3));

            // ignore alpha for now
            for (size_t i = 0; i < size_t(width * height); i++) 
            {
                src[3 * i + 0] = rgba[4 * i + 0];
                src[3 * i + 1] = rgba[4 * i + 1];
                src[3 * i + 2] = rgba[4 * i + 2];
            }

            free(rgba);
        }

    };
}

#endif