#ifndef __DORKTRACER_RAY__
#define __DORKTRACER_RAY__
#pragma once

#include "helperMath.h"

namespace DorkTracer{
    
    struct HitInfo{
        bool hasHit;
        float minT;
        Vec3f normal;
        int matId;
    };

    struct Ray{
        Vec3f origin;
        Vec3f dir;
        HitInfo hitInfo;
        float refractiveIndexOfCurrentMedium;
    };
}

#endif