#include "sphere.hpp"

using namespace DorkTracer;

Sphere::Sphere(std::vector<Vec3f>& vertex_data)
{
    this->vertex_data = vertex_data;
}

bool DorkTracer::Sphere::Intersect(Ray& r){
    Vec3f center = this->vertex_data[center_vertex_id - 1];
    Vec3f oc = r.origin - center;

    float t,t1, t2;

    float c = dot(oc,oc) - (radius * radius);
    float b =  2 * dot(r.dir, oc);
    float a = dot(r.dir, r.dir);
    float delta = b*b - (4 * a * c);

    if(delta < 0.0) return false;
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

        if(t < r.hitInfo.minT && t > 0.0f){
            r.hitInfo.minT = t;
            r.hitInfo.matId = material_id;
            r.hitInfo.hasHit = true;
            // calculate sphere normal at hit point
            Vec3f hitPoint = r.origin + r.dir * t;
            r.hitInfo.normal = makeUnit(hitPoint - center);
            return true;
        }
        else return false;
    }
}