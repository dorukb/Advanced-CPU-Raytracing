#ifndef __DORKTRACER_PERLINTEXTURE__
#define __DORKTRACER_PERLINTEXTURE__

#include "texture.h"
#include <string>

namespace DorkTracer
{
    class PerlinTexture : public Texture
    {
        
    public:

        PerlinTexture(int id, std::string& opMode, float noiseScale, std::string& noiseConversion)
        : Texture(id, opMode)
        {
            this->scale = noiseScale;
            if(noiseConversion == "absval"){
                this->conversionType = Conversion::AbsoluteVal;
            }
            else if(noiseConversion == "linear"){
                this->conversionType = Conversion::Linear;
            }
        };

        virtual Vec3f GetSample(float u, float v)
        {
            // TODO: Implement on the fly perlin noise sampling here.
            return Vec3i(180,30,180);
        }


    private:
        float scale;

        enum Conversion
        {
            Linear,
            AbsoluteVal,
        } conversionType;
    };
}

#endif