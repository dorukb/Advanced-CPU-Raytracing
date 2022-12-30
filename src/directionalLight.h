#ifndef __DORKTRACER_DIRECTIONALLIGHT__
#define __DORKTRACER_DIRECTIONALLIGHT__

#include "helperMath.h"

namespace DorkTracer
{
    class DirectionalLight
    {

    public:
        int id;
        Vec3f dir;
        Vec3f radiance;

        DirectionalLight(int id, Vec3f dir, Vec3f radiance)
        {
            this->id = id;
            this->dir = makeUnit(dir);
            this->radiance = radiance;
        }
        
    };

}

#endif