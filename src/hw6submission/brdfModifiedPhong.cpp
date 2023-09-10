#include "brdfModifiedPhong.h"
#include "material.hpp"
#include "helperMath.h"


#include <iostream>

namespace DorkTracer
{
    BrdfModifiedPhong::BrdfModifiedPhong(float exponent, bool isEnergyConserving, int id) 
        : BRDF(exponent, isEnergyConserving, id)
    { };

    Vec3f BrdfModifiedPhong::apply(Material& mat, Vec3f& kd, Vec3f& ks, Vec3f& w_i, Vec3f& w_o, Vec3f& normal)
    {
         // compute angle theta_i: angle btw normal and w_i
        float angleTheta_i = angleBetweenUnitVectors(w_i, normal);
        if(angleTheta_i >= 90.0f || angleTheta_i < 0) return Vec3f(0,0,0);

        Vec3f perfectReflectionDir = makeUnit((normal * 2.0f * dot(normal, w_i)) - w_i);
        double angleR =  angleBetweenUnitVectors(perfectReflectionDir, w_o);
         if(isnan(angleR))
        {
            std::cout<<"nan angleR!!" << std::endl;
            return Vec3f{0,0,0};
        }  
        if(isEnergyConserving)
        {
            Vec3f kdTerm = kd * (1.0f / M_PI);
            double ksConservationTerm = (exponent+2) / (2*M_PI);
            double cosTerm =  std::pow(cosDeg(angleR), exponent);
            Vec3f ksTerm = ks * (ksConservationTerm * cosTerm);
            Vec3f res =  kdTerm + ksTerm;
            if(isnan(res.x) || isnan(res.y) || isnan(res.z))
            {
                std::cout<<"nan res!!" << std::endl;
                return Vec3f{0,0,0};
            }  
            return res;
        }
        else return kd + ks * std::pow(cosDeg(angleR), exponent);
    }
}