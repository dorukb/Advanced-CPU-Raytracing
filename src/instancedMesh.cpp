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
    bool hasHit = false;
    
    Vec3f rayOriginCache = ray.origin;
    Vec3f rayDirCache = ray.dir;

    if(this->hasMotionBlur)
    {
        ray.origin = ray.origin + motionBlurVector * ray.motionBlurTime;
    }

    // check if it intersects our transformed top-most level BoundingBox.
    if(this->bbox.doesIntersectWith(ray))
    {
        ray.origin = rayOriginCache;

        ray.origin =  Matrix::ApplyTransformToPoint(this->inverseTransform, ray.origin);
        ray.dir = Matrix::ApplyTransformToVector(this->inverseTransform, ray.dir);
        
        if(this->hasMotionBlur)
        {
            ray.origin = ray.origin + motionBlurVector * ray.motionBlurTime;
        }

        // if(dot(motionBlurVector, motionBlurVector) > 0.1)
        // {
        //     if(ray.motionBlurTime < -1.0f){
        //         float motionBlurTime = GetRandom01MotionBlurTime();
        //         ray.motionBlurTime = motionBlurTime;
        //     }
        //     ray.origin = ray.origin - motionBlurVector * ray.motionBlurTime;
        // }

        // if so, delegate to underlying mesh after transforming ray to this instance's space.
        hasHit = this->baseMesh->bvh[0].IntersectBVH(ray, *(this->baseMesh));
        if(hasHit){
            // remember to change this to our id. otherwise baseMesh material will be used for shading.
            ray.hitInfo.hitPoint = rayOriginCache + rayDirCache * ray.hitInfo.minT;
            ray.hitInfo.matId = this->material_id;
            ray.hitInfo.hitShape = this;
            ray.hitInfo.normal = makeUnit(Matrix::ApplyTransformToVector(this->inverseTransposeTransform, ray.hitInfo.normal));
        }
        
        // undo origin&dir changes to local space.
        ray.origin = rayOriginCache;
        ray.dir = rayDirCache;
    }

    return hasHit;    
}
