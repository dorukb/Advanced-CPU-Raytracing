#ifndef __DORKTRACER_MATERIAL__
#define __DORKTRACER_MATERIAL__

#pragma once

#include "helperMath.h"

namespace DorkTracer{

    class Material
    {
    
    public:
    
        // TODO: remove the need for type checking, add Derived classes.
        enum MaterialType{
            Mirror,
            Dielectric,
            Conductor,
            Default
        };

        Vec3f ambient;
        Vec3f diffuse;
        Vec3f specular;
        Vec3f mirror;
        float phong_exponent;
        MaterialType type;
        float refractiveIndex;
        Vec3f absorptionCoefficient;
        float conductorAbsorptionIndex;

    private:
        
    };

}

#endif