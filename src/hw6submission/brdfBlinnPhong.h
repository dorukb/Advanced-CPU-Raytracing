#ifndef __DORKTRACER_ORIGBLINNPHONG__
#define __DORKTRACER_ORIGBLINNPHONG__

#include "brdf.h"

namespace DorkTracer
{
    class BrdfBlinnPhong : public BRDF
    {
    public:
        BrdfBlinnPhong(float exponent, bool isEnergyConserving, int id);
        Vec3f apply(Material& mat, Vec3f& kd, Vec3f& ks, Vec3f& w_i, Vec3f& w_o, Vec3f& normal) override;
    };
}
#endif