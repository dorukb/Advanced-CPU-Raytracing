#ifndef __DORKTRACER_TEXTURE__
#define __DORKTRACER_TEXTURE__

#include "helperMath.h"
#include <string>

namespace DorkTracer
{ 
    enum class DecalMode
    {
        Replace_kd,
        Blend_kd,
        Replace_ks,
        Replace_bg,
        Replace_normal,
        Bump_normal,
        Replace_all, // replace all components (i.e. diffuse, specular, and ambient) of the surface shading color with this texture’s value.

    };

    class Texture
    {
   
    public:
        enum OperationMode
        {
            Blend,
            Replace
        } operationMode;

        enum class Textures
        {
            Diffuse,
            Specular,
            Bump,
            Normal,
            ReplaceAll,
        } type;

        int id;

        Texture(int texId, std::string mode)
        {
            this->id = texId;
            SetOperationMode(mode);
            SetTextureType(mode);
        };

        virtual Vec3f GetDirectSample(int i, int j) = 0;
        virtual Vec3f GetRGBSample(float u, float v) = 0;
        virtual float GetSampleFromWorldPos(float x, float y, float z){
            return 0;
        }

        virtual float GetHeight() = 0;
        virtual float GetWidth() = 0;
        virtual float GetNormalizer() = 0;
        virtual float GetSampleMultiplier()= 0;
        virtual bool IsGenerated(){return false;}

    private:
        void SetTextureType(std::string mode)
        {
            if(mode == "replace_kd" || mode == "blend_kd")
            {
                type = Textures::Diffuse;
            }
            else if(mode == "replace_ks")
            {
                type = Textures::Specular;
            }   
            else if(mode == "replace_normal")
            {
                type = Textures::Normal;
            }
            else if(mode == "bump_normal")
            {
                type = Textures::Bump;
            }
            else if(mode == "replace_all")
            {
                type = Textures::ReplaceAll;
            }
        };

        void SetOperationMode(std::string& mode)
        {
            if(mode == "blend_kd"){
                operationMode = OperationMode::Blend;
            }
            else{
                operationMode = OperationMode::Replace;
            }
        }
    };
}

#endif