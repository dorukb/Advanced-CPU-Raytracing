#ifndef __DORKTRACER_IMAGETEXTURE__
#define __DORKTRACER_IMAGETEXTURE__

#include "texture.h"
#include "image.h"

namespace DorkTracer
{
    class ImageTexture : public Texture
    {

    public:
        ImageTexture(int id,std::string& opMode, std::string& interpolationMode, Image* image)
        : Texture(id, opMode)
        {
            this->img = image;
            this-> interpolationMode = InterpolationMode::Bilinear;
            if(interpolationMode == "nearest"){
                this-> interpolationMode = InterpolationMode::Nearest;
            }
        }

        Vec3f GetSample(float u, float v)
        {
            int i,j;
            if(interpolationMode == InterpolationMode::Nearest){
                i = (int) (u * img->width);
                j = (int) (v * img->height);
            }
            else{
                // TODO: Implement bilinear interpolation.
                i = (int) (u * img->width);
                j = (int) (v * img->height);
            }

            return img->GetSample(i, j);
        }

    private:
        Image* img;
        enum InterpolationMode
        {
            Nearest,
            Bilinear,
        }interpolationMode;
    };
}
#endif