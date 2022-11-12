#include "bvh.hpp"
#include "mesh.hpp"
#include <iostream>

bool DorkTracer::BVH::IntersectBVH (DorkTracer::Ray& ray, DorkTracer::Mesh& ownerMesh)
{    
    if(bbox.doesIntersectWith(ray) == false){
        // this ray entirely misses this bounding box.
        return false;
    }

    bool hasHit = false;
    if(left == nullptr && right == nullptr && faceCount > 0){
        // Leaf node.
        for(int i = firstFace; i < firstFace+faceCount; i++){
            bool faceHit = ownerMesh.IntersectFace(ray, i);
            if(faceHit){
                hasHit = true;
            }
        }
    }   
    else{
        bool hitLeft = left->IntersectBVH(ray, ownerMesh);
        bool hitRight = right->IntersectBVH(ray, ownerMesh);
        if(hitLeft || hitRight){
            hasHit = true;
        }
    }
    // return ray.hitInfo.hasHit;
    return hasHit;
}