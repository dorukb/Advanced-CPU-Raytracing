#ifndef __DORKTRACER_TORRANCESPARROW__
#define __DORKTRACER_TORRANCESPARROW__

#include "brdf.h"

namespace DorkTracer
{
    class BrdfTorranceSparrow : public BRDF
    {
    public:
        bool kdFresnel;
        BrdfTorranceSparrow(float exponent, bool kdFresnel, int id);
        Vec3f apply(Material& mat, Vec3f& kd, Vec3f& ks, Vec3f& w_i, Vec3f& w_o, Vec3f& normal) override;
    private:
        double probabilityOfAngle(double cosAlpha, double exponent);
        double geometryTerm(Vec3f& n, Vec3f& wh, Vec3f& wo, Vec3f& wi);
        double fresnelReflectance(double cosbeta, double refractiveIdx);
    };
}

#endif