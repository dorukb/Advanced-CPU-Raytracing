#ifndef __DORKTRACER_BRDF__
#define __DORKTRACER_BRDF__

#include "helperMath.h"

namespace DorkTracer
{
    class Material;
    
    class BRDF
    {

    public:
        int id;
        float exponent;
        bool isEnergyConserving;

        BRDF(float exponent, bool isEnergyConserving, int id){
            this->exponent = exponent;
            this->isEnergyConserving = isEnergyConserving;
            this->id = id;
        }
        virtual Vec3f apply(Material& mat, Vec3f& kd, Vec3f& ks, Vec3f& w_i, Vec3f& w_o, Vec3f& normal) = 0;

    };
}
#endif