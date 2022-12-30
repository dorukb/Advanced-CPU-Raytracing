#ifndef __DORKTRACER_HDRIMAGE__
#define __DORKTRACER_HDRIMAGE__
#include <string>
#include <iostream>

#include "stb_image.h"
#include "image.h"
#include <vector>
#include <algorithm>

namespace DorkTracer
{
    class HDRImage : public Image
    {

    public:
        HDRImage(std::string filename, int id)
        {
            this->id = id;
            LoadImage(filename);
        }

        Vec3f GetSample(int i, int j) override
        {
            uint32_t imgIdx = 3 * (i + j * width);
            
            Vec3f color;
            float R = src[imgIdx];
            float G = src[imgIdx+1];
            float B = src[imgIdx+2];
            
            // first compute luminance from channels.
            float y_i = 0.2126 * R + 0.7152 * G + 0.0722 * B;

            //TODO: here apply tonemapping to luminance
            float y_o = Tonemap(y_i);

            float sat = 1.0f; // TODO: get saturation from Camera!!
            // compute RGB
            float r_o = clip(y_o * std::pow((R/y_i), sat), 0.0f, 1.0f);
            float g_o = clip(y_o * std::pow((G/y_i), sat), 0.0f, 1.0f);
            float b_o = clip(y_o * std::pow((B/y_i), sat), 0.0f, 1.0f);

            color.x = std::min(255.0f, 255 * std::pow(r_o, (1/2.2f)));
            color.y = std::min(255.0f, 255 * std::pow(g_o, (1/2.2f)));
            color.z = std::min(255.0f, 255 * std::pow(b_o, (1/2.2f)));
            
            return color;
        }
        float clip(float n, float lower, float upper) {
            return std::max(lower, std::min(n, upper));
        }
        int GetSampleGreyscale(int i, int j) override
        {
            // uint32_t imgIdx = i + j * width;
            // return image[imgIdx];
        }

        // Y_i: input luminance
        float Tonemap(float Y_i)
        {
            const float referenceKey = 0.18f; // get from camera TMOOptions!
            const float burnOutPercentile = 0.99f; // get from camera TMOOptions!

            // Eqn 1, find the average key in the scene, done during parsing for performance.

            // Eqn 2. Map the avg lum to corresponding photographic zone.
            float Lxy= referenceKey * Y_i / this->avgLuminance;

            int lastIdx = sortedPixelLuminances.size() - 1;
            int burnThresholdLuminanceIdx = std::min(lastIdx, (int) (burnOutPercentile * lastIdx));
            
            float thresholdLuminance = sortedPixelLuminances[burnThresholdLuminanceIdx];
            float LwhiteSqr = thresholdLuminance * thresholdLuminance; // to keep the lingo same with the paper.

            float tonemappedLuminanceLd = (Lxy * (1 + (Lxy / LwhiteSqr))) / (1 + Lxy);
            // float tonemappedLuminanceLd = Lxy  / (1 + Lxy);

            return tonemappedLuminanceLd;
        }

    protected:
        // unsigned char* image;
        std::vector<float> src;
        std::vector<float> sortedPixelLuminances;
        EXRImage exrImage;

        const float delta = 0.001;
        float avgLuminance = 0.0f;


        void LoadImage(std::string filename) override
        {  
            std::cout << "NumChannels: " << exrImage.num_channels << std::endl;
            this->width = exrImage.width;
            this->height = exrImage.height;
            this->channels = 3;
            
            {
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
                sortedPixelLuminances.resize(size_t(width * height * 3));

                // ignore alpha for now
                
                double R,G,B;
                double logLuminancesSum = 0.0f;
                

                for (size_t i = 0; i < size_t(width * height); i++) {
                    R = src[3 * i + 0] = rgba[4 * i + 0];
                    G = src[3 * i + 1] = rgba[4 * i + 1];
                    B = src[3 * i + 2] = rgba[4 * i + 2];

                    sortedPixelLuminances[3 * i + 0] = rgba[4 * i + 0];
                    sortedPixelLuminances[3 * i + 1] = rgba[4 * i + 1];
                    sortedPixelLuminances[3 * i + 2] = rgba[4 * i + 2];

                    
                    // Compute average luminance here for performance improvements.
                    // first compute luminance from channels.
                    double luminance = 0.2126 * R + 0.7152 * G + 0.0722 * B;
                    // Eqn 1, find the average key in the scene
                    logLuminancesSum += std::log(delta + luminance);
                }

                long pixelCount = width * height;
                this->avgLuminance = std::exp(logLuminancesSum / (float) pixelCount);
                std::cout <<"Avg luminance: " << this->avgLuminance <<" sum: " << logLuminancesSum << std::endl;

                std::sort(sortedPixelLuminances.begin(), sortedPixelLuminances.end());
                free(rgba);
            }
        }

    };
}

#endif