#ifndef __DORKTRACER_MESH__
#define __DORKTRACER_MESH__

#pragma once

#include "shape.hpp"
#include <vector>

namespace DorkTracer{

    class Mesh : public Shape
    {
    
    public:
        // Scene& scene;

        int material_id;
        bool useOwnVertices;
        std::vector<Face> faces;
        std::vector<Vec3f> vertices;

        Mesh();
        virtual bool intersect(Ray& ray);

    private:
        
    };

}

#endif