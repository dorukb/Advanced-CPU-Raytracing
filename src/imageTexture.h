#ifndef __DORKTRACER_IMAGETEXTURE__
#define __DORKTRACER_IMAGETEXTURE__

#include "texture.h"
#include "image.h"

namespace DorkTracer
{
    class ImageTexture : public Texture
    {

    public:
        float normalizer;
        float sampleMultiplier;

        ImageTexture(int id, float normalizer, float sampleMultiplier, std::string& opMode, std::string& interpolationMode, Image* image)
        : Texture(id, opMode)
        {
            this->normalizer = normalizer;
            this->sampleMultiplier = sampleMultiplier;
            this->img = image;
            this-> interpolationMode = InterpolationMode::Bilinear;
            if(interpolationMode == "nearest"){
                this-> interpolationMode = InterpolationMode::Nearest;
            }
        }
        Vec3f GetDirectSample(int i, int j)
        {
            return img->GetSample(i,j);
        }
        Vec3f GetRGBSample(float u, float v)
        {
            int i,j;
            if(interpolationMode == InterpolationMode::Nearest)
            {
                i = (int) (u * img->width);
                j = (int) (v * img->height);
                
                i = std::min(img->width-1, i);
                j = std::min(img->height-1, j);
                return img->GetSample(i, j);
            }
            else return interpolateBilinear(u, v);
        }
        float GetNormalizer(){
            return this->normalizer;
        }
        float GetSampleMultiplier(){
            return this->sampleMultiplier;
        }
        float GetHeight()
        {
            return img->height;
        }
        float GetWidth()
        {
            return img->width;
        }

    private:
        Image* img;
        enum InterpolationMode
        {
            Nearest,
            Bilinear,
        }interpolationMode;
        
        float clip(float n, float lower, float upper) {
            return std::max(lower, std::min(n, upper));
        }

        Vec3f interpolateBilinear(float u, float v) 
        {
            // Clamp the coordinates (x, y) to the bounds of the texture map
            float i = clip(u*img->width, 0.0f, (float)(img->width - 1));
            float j = clip(v*img->height, 0.0f, (float)(img->height - 1));

            float p = std::floor(i);
            float q = std::floor(j);
            float dx = i - p;
            float dy = j - q;


            // Calculate the weights for each of the four nearest points
            float w1 = (1 - dx) * (1 - dy);
            float w2 = dx * (1 - dy);
            float w3 = (1 - dx) * dy;
            float w4 = dx * dy;

            // Calculate the interpolated value of the function at (x, y)
            Vec3f finalColor =  img->GetSample(p,q) * w1 + img->GetSample(p+1, q) * w2
                                + img->GetSample(p, q+1) * w3 + img->GetSample(p+1, q+1) * w4;
            return finalColor;
        }       
    };
}
#endif