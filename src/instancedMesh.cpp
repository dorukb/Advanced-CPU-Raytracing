#include "instancedMesh.hpp"
#include <iostream>

using namespace DorkTracer;

DorkTracer::InstancedMesh::InstancedMesh(Mesh* baseMesh)
{
    this->baseMesh = baseMesh;
}

void DorkTracer::InstancedMesh::SetMaterial(int matId){
    this->material_id = matId;
}


bool DorkTracer::InstancedMesh::Intersect(Ray& ray)
{    
    Vec3f rayOriginCache = ray.origin;
    Vec3f rayDirCache = ray.dir;
    bool hasHit = false;
    
    // check if it intersects our transformed top-most level BoundingBox.
    // if(this->bbox.doesIntersectWith(ray))
    // {
        
        // Transform the ray into our local space.
        Vec4f rayOrigin(ray.origin, 1.0f);
        Vec4f rayDir(ray.dir, 0.0f);
        
        ray.origin =  Matrix::ApplyTransform(this->inverseTransform, rayOrigin);
        ray.dir = Matrix::ApplyTransform(this->inverseTransform, rayDir);
        
        // if so, delegate to underlying mesh after transforming ray to this instance's space.
        hasHit = this->baseMesh->bvh[0].IntersectBVH(ray, *(this->baseMesh));
        if(hasHit){
            // remember to change this to our id. otherwise baseMesh material will be used for shading.
            ray.hitInfo.matId = this->material_id;
            ray.hitInfo.normal = makeUnit(Matrix::ApplyTransform(this->inverseTransposeTransform, Vec4f(ray.hitInfo.normal, 0.0f)));
        }
    // }
    // undo origin&dir changes to local space.
    ray.origin = rayOriginCache;
    ray.dir = rayDirCache;

    return hasHit;    
}
