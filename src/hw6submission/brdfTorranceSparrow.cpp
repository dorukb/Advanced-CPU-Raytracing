#include "brdfTorranceSparrow.h"
#include "helperMath.h"
#include "material.hpp"
#include <algorithm>

namespace DorkTracer
{
    // this BRDF is normalized/energy conserving by default.
    BrdfTorranceSparrow::BrdfTorranceSparrow(float exponent, bool kdFresnel, int id) 
        : BRDF(exponent, true, id)
    { 
        this->kdFresnel = kdFresnel;
    };

    Vec3f BrdfTorranceSparrow::apply(Material& mat, Vec3f& kd, Vec3f& ks, Vec3f& w_i, Vec3f& w_o, Vec3f& normal)
    {
        // compute angle theta_i: angle btw normal and w_i
        float angleTheta_i = angleBetweenUnitVectors(w_i, normal);
        if(angleTheta_i >= 90.0f) return Vec3f(0,0,0);

        // energy conserving by default
        Vec3f half = (w_i + w_o) / len(w_i + w_o);

        double d = probabilityOfAngle(dot(half, normal), exponent);
        double f = fresnelReflectance(dot(half, w_o), mat.refractiveIndex);
        double g = geometryTerm(normal, half, w_o, w_i);

        double kdCoeff = (1.0f/M_PI);
        if(kdFresnel){
            kdCoeff *= (1-f);
        }
        Vec3f kdTerm = kd * kdCoeff;

        double costheta = dot(normal, w_i);
        double cosphi = dot(normal, w_o);
        Vec3f ksTerm = ks * ((d * f * g) / (4*costheta*cosphi));

        return kdTerm + ksTerm;
    }


    // The "D" term, simulates orientation distribution of the micro-facet surface.
    double BrdfTorranceSparrow::probabilityOfAngle(double cosAlpha, double exponent)
    {
        return (exponent+2) * std::pow(cosAlpha, exponent) / (2*M_PI);
    }
    double BrdfTorranceSparrow::geometryTerm(Vec3f& n, Vec3f& h, Vec3f& wo, Vec3f& wi)
    {
        double ndoth = dot(n,h);
        double ndotwo = dot(n,wo);
        double ndotwi = dot(n,wi);
        double wodoth = dot(wo, h);
        return std::min(1.0, std::min(2.0f * ndoth * ndotwo / wodoth, 2.0 * ndoth * ndotwi / wodoth));
    }
    double BrdfTorranceSparrow::fresnelReflectance(double cosbeta, double refractiveIdx)
    {
        double r0 = std::pow(refractiveIdx-1, 2) / std::pow(refractiveIdx+1, 2);
        return r0 + (1.0-r0) * std::pow((1.0 - cosbeta) , 5.0);
    }
}