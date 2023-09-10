#ifndef __DORKTRACER_INSTANCED_MESH__
#define __DORKTRACER_INSTANCED_MESH__

#include "mesh.hpp"

namespace DorkTracer{

    class InstancedMesh : public Shape
    {
    
    public:
        
        Mesh* baseMesh;
        BoundingBox bbox;

        InstancedMesh(Mesh* baseMesh);

        void SetMaterial(int matId);
        bool Intersect(Ray& ray);

    private:
        int material_id;
        uint nextFreeNodeIdx = 0;
    };

}

#endif