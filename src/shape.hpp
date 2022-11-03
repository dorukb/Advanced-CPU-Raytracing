#ifndef __DORKTRACER_SHAPE__
#define __DORKTRACER_SHAPE__

#pragma once

#include "ray.hpp"
#include "material.hpp"
#include <vector>

namespace DorkTracer{

    class Shape
    {
        
    public:
        Shape(){}

        virtual bool intersect(Ray& ray){};

        Material& getMaterial(){
            return *mat;
        }

    private:

        Material* mat;
    };

    struct Face
    {
        int v0_id;
        int v1_id;
        int v2_id;
        Vec3f n;
    };
}

#endif