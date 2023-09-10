#ifndef __DORKTRACER_MESHLIGHT__
#define __DORKTRACER_MESHLIGHT__

#include "mesh.hpp"
#include <random>
#include "matrix.hpp"

namespace DorkTracer
{
    class MeshLight : public Mesh
    {
    public:
        MeshLight(std::vector<Vec3f>& vertices, std::vector<Vec2f>& uv, Vec3f radiance) 
        : Mesh(vertices, uv)
        {
            this->radiance = radiance;
            this->sampleGen = std::mt19937();
            this->uniform01 = std::uniform_real_distribution<>(0.0, 1.0);
        }
        void setupFaceSelectionDistro(int faceCount)
        {
            this->faceDistro = std::uniform_int_distribution<>(0, faceCount);
        }

        Vec3f radiance;
        // double area;
        void getSample(Vec3f& pos, Vec3f& normal, double& weight)
        {
            Face& face =  faces[faceDistro(sampleGen)];

            double selectionWeight = face.area / this->surfaceArea;

            double rand1 = uniform01(sampleGen);
            double rand2 = uniform01(sampleGen);
            Vec3f& a = GetVertex(face.v0_id);
            Vec3f& b = GetVertex(face.v1_id);
            Vec3f& c = GetVertex(face.v2_id);

            Vec3f q,p;
            q = b * (1-rand2) + c * rand2;
            pos = a * (1- std::sqrt(rand1)) + q * std::sqrt(rand1);

            // Assign the "out" variables.
            pos = Matrix::ApplyTransformToPoint(this->transform,pos);
            normal = face.n;
            weight = selectionWeight;
        }
    private:
    
        std::mt19937 sampleGen;
        std::uniform_int_distribution<> faceDistro;
        std::uniform_real_distribution<> uniform01;
    };
}
    
#endif