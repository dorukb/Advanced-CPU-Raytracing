#ifndef __DORKTRACER_MATERIAL__
#define __DORKTRACER_MATERIAL__

#include "brdf.h"

namespace DorkTracer{

    class Material
    {
    
    public:
    
        // TODO: remove the need for type checking, add Derived classes.
        enum MaterialType{
            Mirror,
            Dielectric,
            Conductor,
            Emissive,
            Default
        };

        int id;
        BRDF* brdf = nullptr;
        
        Vec3f ambient;
        Vec3f diffuse;
        Vec3f specular;
        Vec3f mirror;
        float phong_exponent;
        MaterialType type;
        float refractiveIndex;
        Vec3f absorptionCoefficient;
        float conductorAbsorptionIndex;
        float roughness;
        Vec3f radiance; // only for emissive mats.
        
        void SetBRDF(BRDF* brdf)
        {
            this->brdf = brdf;
        }
        bool HasBRDF()
        {
            return brdf != nullptr;
        }
        
    private:
        
    };

}

#endif