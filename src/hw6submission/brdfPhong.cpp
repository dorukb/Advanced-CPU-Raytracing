#include "brdfPhong.h"
#include "material.hpp"
#include "helperMath.h"

namespace DorkTracer
{
    BrdfPhong::BrdfPhong(float exponent, bool isEnergyConserving, int id) 
        : BRDF(exponent, isEnergyConserving, id)
    { };

    Vec3f BrdfPhong::apply(Material& mat, Vec3f& kd, Vec3f& ks, Vec3f& w_i, Vec3f& w_o, Vec3f& normal)
    {
         // compute angle theta_i: angle btw normal and w_i
        float angleTheta_i = angleBetweenUnitVectors(w_i, normal);
        if(angleTheta_i >= 90.0f || angleTheta_i < 0) return Vec3f(0,0,0);

        Vec3f perfectReflectionDir = makeUnit((normal * 2.0f * dot(normal, w_i)) - w_i);
        double angleR =  angleBetweenUnitVectors(perfectReflectionDir, w_o);
        return kd + ks * (std::pow(cosDeg(angleR), exponent) / cosDeg(angleTheta_i));
    }
}