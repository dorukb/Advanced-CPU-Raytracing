#ifndef __DORKTRACER_SPHERE__
#define __DORKTRACER_SPHERE__


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
        virtual bool Intersect(Ray& ray);

    private:
        void GetTangentAndBitangentAroundPoint(Vec3f& p, float radius, float phi, float theta, Vec3f& tan, Vec3f& bitan);
        
    };

}

#endif