#ifndef __DORKTRACER_RAY__
#define __DORKTRACER_RAY__

#include "helperMath.h"

namespace DorkTracer
{
    
    struct HitInfo
    {
        bool hasHit;
        int matId;

        float minT;
        Vec3f normal;
        Vec3f hitPoint;
    };

    struct Ray
    {
        Vec3f origin;
        Vec3f dir;
        HitInfo hitInfo;
        float refractiveIndexOfCurrentMedium;

        float motionBlurTime;
        float lightSampleX, lightSampleY;
    };
}

#endif