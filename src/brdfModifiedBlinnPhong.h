#ifndef __DORKTRACER_MODBLINNPHONG__
#define __DORKTRACER_MODBLINNPHONG__

#include "brdf.h"
#include "material.hpp"

namespace DorkTracer
{
    class BrdfModifiedBlinnPhong : public BRDF
    {

    public:
        BrdfModifiedBlinnPhong(float exponent, bool isEnergyConserving) : BRDF(exponent, isEnergyConserving)
        { };

        Vec3f apply(Material& mat, Vec3f& w_i, Vec3f& w_o, Vec3f& normal) override
        {
            // compute angle theta_i: angle btw normal and w_i
            float angleTheta_i = angleBetweenUnitVectors(w_i, normal);
            if(angleTheta_i <= 0.0001) return Vec3f(0,0,0);

            Vec3f half = (w_i + w_o) / len(w_i + w_o);
            double angleBtwHalfandNormal =  angleBetweenUnitVectors(half, normal);

            if(isEnergyConserving) // or normalized
            {
                Vec3f kdTerm = mat.diffuse * (1.0f / M_PI);
                double ksConservationTerm = (exponent+8) / (8*M_PI);
                double cosTerm =  std::pow(std::cos(angleBtwHalfandNormal), exponent);
                Vec3f ksTerm = mat.specular * (ksConservationTerm * cosTerm);
                return kdTerm + ksTerm;
            }
            else return mat.diffuse + mat.specular * std::pow(std::cos(angleBtwHalfandNormal), exponent);
        }
            // Vec3f perfectReflectionDir = makeUnit((normal * 2.0f * dot(normal, w_i)) - w_i);
    };
}
#endif