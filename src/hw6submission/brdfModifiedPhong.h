#ifndef __DORKTRACER_MODPHONG__
#define __DORKTRACER_MODPHONG__

#include "brdf.h"

namespace DorkTracer
{
    class BrdfModifiedPhong : public BRDF
    {
    public:
        BrdfModifiedPhong(float exponent, bool isEnergyConserving, int id);
        Vec3f apply(Material& mat, Vec3f& kd, Vec3f& ks, Vec3f& w_i, Vec3f& w_o, Vec3f& normal) override;
    };
}
#endif