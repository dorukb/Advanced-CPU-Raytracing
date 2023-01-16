#ifndef __DORKTRACER_MESH__
#define __DORKTRACER_MESH__

#include "shape.hpp"
#include <vector>
#include "bvh.hpp"
#include <random>

namespace DorkTracer{

    class Mesh : public Shape
    {
    
    public:
        std::vector<Face> faces;
        
        BoundingBox bbox;
        std::vector<BVH> bvh;

        Mesh(std::vector<Vec3f>& vertices, std::vector<Vec2f>& uv);

        virtual bool Intersect(Ray& ray);

        bool IntersectFace(Ray& ray, Face& face);
        bool IntersectFace(Ray& ray, uint32_t faceIdx);
        void ConstructBVH();
        void SetAccessOffsets(int vertexOffset, int textureOffset);
        Vec2f& GetUv(int idx);
        Vec3f& GetVertex(int idx);
        
        void AddVertex(Vec3f vert){
            vertices.push_back(vert);
        }

    private:
        std::vector<Vec3f> vertices;
        std::vector<Vec2f> uv;

        uint nextFreeNodeIdx = 0;
        int vertexOffset, textureOffset;
        
        bool IsBackface(Face& face, Vec3f& rayDir);
        void RecomputeBoundingBox(uint32_t nodeIdx);
        void RecursiveBVHBuild(uint32_t nodeIdx);
        void GetTangentAndBitangentForTriangle(Vec3f& vert0, Vec3f& vert1, Vec3f& vert2, Vec2f& v0_uv, Vec2f& v1_uv, Vec2f& v2_uv, Vec3f& tan, Vec3f& bitan);
        // bool DoesIntersectTriangle(Ray& ray, Vec3f& v0, Vec3f& v1, Vec3f& v2, float& t, float currMinT);
        float inline GetFloorForTiledUV(float x);

    };

}

#endif