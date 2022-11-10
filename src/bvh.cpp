#include "bvh.hpp"
#include "mesh.hpp"
#include <iostream>

bool DorkTracer::BVH::IntersectBVH (DorkTracer::Ray& ray, DorkTracer::Mesh& ownerMesh)
{    
    if(bbox.doesIntersectWith(ray) == false){
        // this ray entirely misses this bounding box.
        return false;
    }

    if(left == nullptr && right == nullptr && faceCount > 0){
        // Leaf node.
        for(int i = firstFace; i < firstFace+faceCount; i++){
            ownerMesh.IntersectFace(ray, i);
        }
    }   
    else{
        left->IntersectBVH(ray, ownerMesh);
        right->IntersectBVH(ray, ownerMesh);
    }
    return ray.hitInfo.hasHit;
}