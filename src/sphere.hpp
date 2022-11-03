#ifndef __DORKTRACER_SPHERE__
#define __DORKTRACER_SPHERE__

#pragma once

#include "shape.hpp"

namespace DorkTracer
{
    class Sphere : public Shape
    {
    
    public:
        int material_id;
        int center_vertex_id;
        float radius;
        
        std::vector<Vec3f> vertex_data;

        Sphere(std::vector<Vec3f>& vertex_data);
        virtual bool intersect(Ray& ray);

    private:
        
    };

}

#endif