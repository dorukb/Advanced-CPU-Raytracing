#ifndef __DORKTRACER_MESH__
#define __DORKTRACER_MESH__

#include "shape.hpp"
#include <vector>
#include "bvh.hpp"

namespace DorkTracer{

    class Mesh : public Shape
    {
    
    public:
        int material_id;
        std::vector<Face> faces;
        std::vector<Vec3f> vertices;
        Material* material;

        BoundingBox bbox;
        std::vector<BVH> bvh;

        Mesh(std::vector<Vec3f>& vertices);
        virtual bool Intersect(Ray& ray);
        bool IntersectFace(Ray& ray, Face& face);
        bool IntersectFace(Ray& ray, uint32_t faceIdx);

        bool IsBackface(Face& face, Vec3f& rayDir);
        bool DoesIntersectTriangle(Ray& ray, Vec3f& v0, Vec3f& v1, Vec3f& v2, float& t);
        void SetVertices(std::vector<Vec3f>& vertices);

        void ConstructBVH();
        void RecomputeBoundingBox(uint32_t nodeIdx);
        void RecursiveBVHBuild(uint32_t nodeIdx);

    private:
        uint nextFreeNodeIdx = 0;
        
    };

}

#endif