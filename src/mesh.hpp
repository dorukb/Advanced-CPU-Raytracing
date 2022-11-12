#ifndef __DORKTRACER_MESH__
#define __DORKTRACER_MESH__

#include "shape.hpp"
#include <vector>
#include "bvh.hpp"

namespace DorkTracer{

    class Mesh : public Shape
    {
    
    public:
        std::vector<Face> faces;
        std::vector<Vec3f> vertices;
        
        BoundingBox bbox;
        std::vector<BVH> bvh;

        Mesh(std::vector<Vec3f>& vertices);

        virtual bool Intersect(Ray& ray);
        int GetMaterial();
        void SetMaterial(int matId);

        bool IntersectFace(Ray& ray, Face& face);
        bool IntersectFace(Ray& ray, uint32_t faceIdx);
        void ConstructBVH();

    private:
        int material_id;
        uint nextFreeNodeIdx = 0;
        
        bool IsBackface(Face& face, Vec3f& rayDir);
        void RecomputeBoundingBox(uint32_t nodeIdx);
        void RecursiveBVHBuild(uint32_t nodeIdx);
        bool DoesIntersectTriangle(Ray& ray, Vec3f& v0, Vec3f& v1, Vec3f& v2, float& t);
    };

}

#endif