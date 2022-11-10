#ifndef __DORKTRACER_BVH__
#define __DORKTRACER_BVH__

#include <stdint.h>
#include "shape.hpp"

namespace DorkTracer{

    class Mesh;

    class BVH
    {

    public:
    
        BVH* left;
        BVH* right;
        BoundingBox bbox;
        uint32_t firstFace, faceCount;

        BVH(){
            left = right = nullptr;
        }
        
        bool IntersectBVH(Ray&, Mesh& ownerMesh);
    };
}

#endif