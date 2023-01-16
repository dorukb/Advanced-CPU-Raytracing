#ifndef __DORKTRACER_BRDF__
#define __DORKTRACER_BRDF__

#include "helperMath.h"
#include <random>

namespace DorkTracer{
    class BRDF
    {

    public:
        float exponent;
        bool isEnergyConserving;

        BRDF(float exponent, bool isEnergyConserving){
            this->exponent = exponent;
            this->isEnergyConserving = isEnergyConserving;
        }
        virtual Vec3f apply(Material& mat, Vec3f& w_i, Vec3f& w_o, Vec3f& normal);

    };
}
#endif