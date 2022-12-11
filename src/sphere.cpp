#include "sphere.hpp"

using namespace DorkTracer;

Sphere::Sphere(std::vector<Vec3f>& vertex_data)
{
    this->vertex_data = vertex_data;
}

bool DorkTracer::Sphere::Intersect(Ray& r){
    Vec3f center = this->vertex_data[center_vertex_id - 1];

    Vec3f worldCenter =  Matrix::ApplyTransform(this->transform, Vec4f(center,1.0f));
    // Transform the ray into our local space.
    Vec3f rayOriginCache = r.origin;
    Vec3f rayDirCache = r.dir;
    
    Vec4f rayOrg(r.origin, 1.0f);
    Vec4f rayDir(r.dir, 0.0f);
    
    Vec3f transformedOrigin = Matrix::ApplyTransform(this->inverseTransform, rayOrg);
    Vec3f transformedDir = Matrix::ApplyTransform(this->inverseTransform, rayDir);

    r.origin = transformedOrigin;
    // r.dir = makeUnit(transformedDir);
    r.dir = transformedDir;

    Vec3f oc = r.origin - center;
    float t,t1, t2;

    float c = dot(oc,oc) - (radius * radius);
    float b =  2 * dot(r.dir, oc);
    float a = dot(r.dir, r.dir);
    float delta = b*b - (4 * a * c);

    if(delta < 0.0){
        r.origin = rayOriginCache;
        r.dir = rayDirCache;
        return false;
    }
    else
    {
        delta = sqrtf(delta);
        a = 2.0 * a;
        float t1 = (-b + delta) / a;
        float t2 = (-b - delta) / a;
        t = t1 < t2 ? t1 : t2;

        // Find the minimum t which is > 0.0f
        if(t1 < t2){
            if(t1 > 0.0f){
                t = t1;
            }
            else t = t2;
        }
        else if(t2 < t1){
            if(t2 > 0.0f){
                t = t2;
            }
            else t = t1;
        }
        // todo: consider ignoring negative t before selecting the min.

        Vec3f localhitPoint = r.origin + r.dir * t;
        // // undo ray transformation
        r.origin = rayOriginCache;
        r.dir = rayDirCache;

        if(t < r.hitInfo.minT && t > 0.0f){
            r.hitInfo.minT = t;
            r.hitInfo.matId = material_id;
            r.hitInfo.hitShape = this;
            r.hitInfo.hasHit = true;
            // calculate sphere normal at hit point
            Vec3f worldHitPoint = r.origin + r.dir * t;
            r.hitInfo.hitPoint = worldHitPoint;

            r.hitInfo.normal = makeUnit(localhitPoint - center);
            r.hitInfo.normal = makeUnit(Matrix::ApplyTransform(this->inverseTransposeTransform, Vec4f(r.hitInfo.normal, 0.0f)));

            // r.hitInfo.normal = makeUnit(worldHitPoint - worldCenter);
            // r.hitInfo.normal = worldHitPoint - worldCenter;

            return true;
        }
        else return false;
    }
}