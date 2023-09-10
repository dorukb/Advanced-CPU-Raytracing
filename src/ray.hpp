#ifndef __DORKTRACER_RAY__
#define __DORKTRACER_RAY__

#include "helperMath.h"

namespace DorkTracer
{
    class Shape;

    struct HitInfo
    {
        bool hasHit;
        int matId;

        float minT;
        Vec3f normal;
        Vec3f hitPoint;
        Vec2f hitUV;
        Shape* hitShape;
    };

    struct Ray
    {
        Vec3f origin;
        Vec3f dir;
        HitInfo hitInfo;
        float refractiveIndexOfCurrentMedium;

        float motionBlurTime;
        float lightSampleX, lightSampleY;
        Vec3f throughput{1.0f, 1.0f, 1.0f};
    };
}

#endif