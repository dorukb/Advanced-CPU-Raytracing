#ifndef __DORKTRACER_SPOTLIGHT__
#define __DORKTRACER_SPOTLIGHT__

#include "light.h"
#include "helperMath.h"
#define DEG2RAD (M_PI / 180.0f)

namespace DorkTracer
{
    class SpotLight : public Light
    {

    public:
        Vec3f pos;
        Vec3f dir;
        Vec3f intensity;
        float coverageAngle;
        float falloffAngle;

        SpotLight(int id, Vec3f pos, Vec3f dir, Vec3f intensity, float coverageAngle, float falloffAngle)
        : Light(id)
        {
            this->pos = pos;
            this->dir = makeUnit(dir);
            this->intensity = intensity;
            this->coverageAngle = coverageAngle;
            this->falloffAngle = falloffAngle;

            this->cosHalfCoverage = std::cos((coverageAngle * DEG2RAD / 2.0f));
            this->cosHalfFalloff = std::cos((falloffAngle * DEG2RAD  / 2.0f));
        }

        Vec3f GetIrradiance(Vec3f point)
        {
            float distToPoint = len(point - pos);
            Vec3f toPoint = (point - pos) / distToPoint;
            double alpha = angleBetweenUnitVectors(dir, toPoint);

            if(alpha <= 0 || alpha > (coverageAngle/2.0f))
            {
                return Vec3f(0,0,0);
            }

            // regular attenuation.
            float distSqr = distToPoint * distToPoint;
            Vec3f irradiance = intensity / distSqr;

            if(alpha > (falloffAngle/2.0f))
            {
                // fall-off zone, apply fall off depending on alpha.
                double cosAlpha = std::cos(alpha * DEG2RAD);
                double s = std::pow((cosAlpha-cosHalfCoverage) / (cosHalfFalloff - cosHalfCoverage), 4.0f);
                irradiance = irradiance * (float) s;
            }

            return irradiance;
        }

    private:
        double cosHalfFalloff;
        double cosHalfCoverage;

    };
}

#endif