#ifndef __DORKTRACER_PHONG__
#define __DORKTRACER_PHONG__

#include "brdf.h"

namespace DorkTracer
{
    class BrdfPhong : public BRDF
    {
    public:
        BrdfPhong(float exponent, bool isEnergyConserving, int id);
        Vec3f apply(Material& mat, Vec3f& kd, Vec3f& ks, Vec3f& w_i, Vec3f& w_o, Vec3f& normal) override;
    };
}
#endif