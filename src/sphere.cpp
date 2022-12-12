#include "sphere.hpp"
using namespace DorkTracer;

Sphere::Sphere(std::vector<Vec3f>& vertex_data)
{
    this->vertex_data = vertex_data;
}

bool DorkTracer::Sphere::Intersect(Ray& r)
{
    Vec3f center = this->vertex_data[center_vertex_id - 1];

    Vec3f rayOriginCache = r.origin;
    Vec3f rayDirCache = r.dir;
    
    // Transform the ray into our local space.
    r.origin = Matrix::ApplyTransformToPoint(this->inverseTransform, r.origin);
    r.dir = Matrix::ApplyTransformToVector(this->inverseTransform, r.dir);

    Vec3f oc = r.origin - center;
    float t, t1, t2;

    float c = dot(oc,oc) - (radius * radius);
    float b =  2 * dot(r.dir, oc);
    float a = dot(r.dir, r.dir);
    float delta = b*b - (4 * a * c);

    if(delta < 0.0f){
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

        Vec3f localhitPoint = r.origin + r.dir * t;

        // undo ray transformation
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

            // compute hit point UV.
            Vec3f p = worldHitPoint - center;
            float phi = std::atan2(p.z, p.x);
            float theta = std::acos(p.y / radius);

            float u = (-phi + M_PI) / (2.0f * M_PI);
            float v = theta / M_PI;

            r.hitInfo.hitUV.x = u;
            r.hitInfo.hitUV.y = v;

            r.hitInfo.normal = makeUnit(localhitPoint - center);
            r.hitInfo.normal = makeUnit(Matrix::ApplyTransform(this->inverseTransposeTransform, Vec4f(r.hitInfo.normal, 0.0f)));

            return true;
        }
        else return false;
    }
}