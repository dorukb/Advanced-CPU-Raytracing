#ifndef __DORKTRACER_POINTLIGHT__
#define __DORKTRACER_POINTLIGHT__

#include "light.h"
#include "helperMath.h"

namespace DorkTracer{

    class PointLight : public Light
    {

    public:
        Vec3f position;
        Vec3f intensity;
        PointLight(int id, Vec3f pos, Vec3f intensity)
        : Light(id)
        {
            this->position = pos;
            this->intensity = intensity;
        };
    };
}

#endif