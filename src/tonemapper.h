#include <string>
#include <algorithm>

#include "HDRImage.h"

namespace DorkTracer
{
    class Tonemapper
    {
    public:
        std::string operatorType;
        float keyValue;
        float burnPerct;
        float saturation;
        float gamma;
        float delta = 0.001f;

        Tonemapper(std::string opType, float keyVal, float burnPerct, float saturation, float gamma)
        {
            this->operatorType = opType;
            this->keyValue = keyVal;
            this->burnPerct = burnPerct;
            this->saturation = saturation;
            this->gamma = gamma;
        }
                
        // Tonemaps the whole hdr image, writes to the ldrDest buffer.
        void Tonemap(int width, int height, float* hdrSrc, unsigned char* ldrDest)
        { 
            // compute avg luminances
            // create sorted luminance list for burnout.
            
            sortedPixelLuminances.resize(size_t(width * height * 3));
            double logLuminancesSum = 0.0f;
            for (size_t i = 0; i < size_t(width * height); i++) 
            {
                double r = sortedPixelLuminances[3 * i + 0] = hdrSrc[3 * i + 0];
                double g = sortedPixelLuminances[3 * i + 1] = hdrSrc[3 * i + 1];
                double b = sortedPixelLuminances[3 * i + 2] = hdrSrc[3 * i + 2];

                double luminance = 0.2126 * r + 0.7152 * g + 0.0722 * b;
                // Eqn 1, find the average key in the scene
                logLuminancesSum += std::log(delta + luminance);
            }

            long pixelCount = width * height;
            this->avgLuminance = std::exp(logLuminancesSum / (double) pixelCount);
            std::cout <<"Avg luminance: " << this->avgLuminance <<" sum: " << logLuminancesSum << std::endl;

            std::sort(sortedPixelLuminances.begin(), sortedPixelLuminances.end());


            // std::cout << "Key: " << keyValue << " burn perct: " << burnPerct << std::endl;
            // Tonemap each pixel, writes the tonemapped color to the corresponding index of ldrDest.
            for (size_t i = 0; i < size_t(width * height); i++) 
            {
                TonemapPixel(i, hdrSrc, ldrDest);
            }
            
        }

    private:
        double avgLuminance;
        std::vector<float> sortedPixelLuminances;


        void TonemapPixel(int i, float* hdrSrc, unsigned char* ldrDest)
        {
            double R = hdrSrc[3 * i + 0];
            double G = hdrSrc[3 * i + 1];
            double B = hdrSrc[3 * i + 2];
            // first compute luminance from channels.
            double y_i = 0.2126 * R + 0.7152 * G + 0.0722 * B;

            double y_o = Reinhard(y_i);

            // compute new RGB
            double r_o = clip(y_o * std::pow((R/y_i), saturation), 0.0f, 1.0f);
            double g_o = clip(y_o * std::pow((G/y_i), saturation), 0.0f, 1.0f);
            double b_o = clip(y_o * std::pow((B/y_i), saturation), 0.0f, 1.0f);

            double gammaInv = 1.0f / gamma;
            Vec3i color;
            color.x = std::floor(std::min(255.0, 255 * std::pow(r_o, gammaInv)));
            color.y = std::floor(std::min(255.0, 255 * std::pow(g_o, gammaInv)));
            color.z = std::floor(std::min(255.0, 255 * std::pow(b_o, gammaInv)));
            
            ldrDest[3*i+0] = color.x;
            ldrDest[3*i+1] = color.y;
            ldrDest[3*i+2] = color.z;
        }
        // Y_i: input luminance
        float Reinhard(double Y_i)
        {
            // Eqn 1, find the average key in the scene, done beforehand for performance.

            // Eqn 2. Map the avg lum to corresponding photographic zone.
            double Lxy= (keyValue * Y_i) / this->avgLuminance;

            if(burnPerct > 0.01)
            {
                // i.e, if burnPerct = 3, then the highest 3% of the luminances should map to pure white.
                float thresholdPerct = (100.0f- burnPerct) / 100;
                int lastIdx = sortedPixelLuminances.size() - 1;

                int burnThresholdLuminanceIdx = std::min(lastIdx, (int) (thresholdPerct * lastIdx));
                double thresholdLuminance = sortedPixelLuminances[burnThresholdLuminanceIdx];
                thresholdLuminance = thresholdLuminance * keyValue / this->avgLuminance;

                double LwhiteSqr = thresholdLuminance * thresholdLuminance; // to keep the lingo same with the paper.
                double res = (Lxy * (1 + (Lxy / LwhiteSqr))) / (1.0f + Lxy);
                return res;
            }
            else
            {
                return  Lxy / (1 + Lxy);
            }
        }

        float clip(float n, float lower, float upper) 
        {
            return std::max(lower, std::min(n, upper));
        }

    };
}
