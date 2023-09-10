#include "brdfModifiedBlinnPhong.h"
#include "material.hpp"
#include "helperMath.h"

namespace DorkTracer
{
    BrdfModifiedBlinnPhong::BrdfModifiedBlinnPhong(float exponent, bool isEnergyConserving, int id) 
        : BRDF(exponent, isEnergyConserving, id)
    { };

    Vec3f BrdfModifiedBlinnPhong::apply(Material& mat, Vec3f& kd, Vec3f& ks, Vec3f& w_i, Vec3f& w_o, Vec3f& normal)
    {
        // compute angle theta_i: angle btw normal and w_i
        float angleTheta_i = angleBetweenUnitVectors(w_i, normal);
        if(angleTheta_i >= 90.0f) return Vec3f(0,0,0);

        Vec3f half = (w_i + w_o) / len(w_i + w_o);
        double angleBtwHalfandNormal =  angleBetweenUnitVectors(half, normal);

        if(isEnergyConserving) // or normalized
        {
            Vec3f kdTerm = kd * (1.0f / M_PI);
            double ksConservationTerm = (exponent+8) / (8*M_PI);
            double cosTerm =  std::pow(cosDeg(angleBtwHalfandNormal), exponent);
            Vec3f ksTerm = ks * (ksConservationTerm * cosTerm);
            return kdTerm + ksTerm;
        }
        else return kd + ks * std::pow(cosDeg(angleBtwHalfandNormal), exponent);
    }
            // Vec3f perfectReflectionDir = makeUnit((normal * 2.0f * dot(normal, w_i)) - w_i);
}