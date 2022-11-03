#ifndef __DORKTRACER_MESH__
#define __DORKTRACER_MESH__

#include "shape.hpp"
#include <vector>

namespace DorkTracer{

    class Mesh : public Shape
    {
    
    public:
        int material_id;
        std::vector<Face> faces;
        std::vector<Vec3f> vertices;
        Material* material;

        Mesh(std::vector<Vec3f>& vertices);
        virtual void Intersect(Ray& ray);
        bool IntersectFace(Ray& ray, Face& face);
        bool IsBackface(Face& face, Vec3f& rayDir);
        bool DoesIntersectTriangle(Ray& ray, Vec3f& v0, Vec3f& v1, Vec3f& v2, float& t);
        void SetVertices(std::vector<Vec3f>& vertices);
    private:
        
    };

}

#endif