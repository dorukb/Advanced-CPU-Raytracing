#include "mesh.hpp"

using namespace DorkTracer;

Mesh::Mesh(std::vector<Vec3f>& vertices){
    this->vertices = vertices;
}
void DorkTracer::Mesh::Intersect(Ray& ray)
{    
    bool bfcEnabled;

    // BFC should be disabled for refractive materials, they have both "faces" by default.
    if(material->type == Material::Dielectric){
        bfcEnabled = false;
    }else{
        bfcEnabled = true;
    }

    for(int i = 0; i < faces.size(); i++)
    {
        if(bfcEnabled && IsBackface(faces[i], ray.dir)){
            continue;
        }

        IntersectFace(ray, faces[i]);
    }

}

bool DorkTracer::Mesh::IntersectFace(Ray& ray, Face& face)
{
    Vec3f& v0 = vertices[face.v0_id-1];
    Vec3f& v1 = vertices[face.v1_id-1];
    Vec3f& v2 = vertices[face.v2_id-1];
    float newT  = 99999;
    bool hasIntersected = DoesIntersectTriangle(ray, v0, v1, v2, newT);

    // If its the closest intersection so far, update ray HitInfo.
    if(hasIntersected && newT < ray.hitInfo.minT){
        ray.hitInfo.minT = newT;
        ray.hitInfo.hasHit = true;
        ray.hitInfo.normal = face.n;
        ray.hitInfo.matId = this->material_id;
    }
    return hasIntersected;
}
bool DorkTracer::Mesh::IsBackface(Face& face, Vec3f& rayDir)
{
    return dot(face.n, rayDir) > 0.0f;
}

bool DorkTracer::Mesh::DoesIntersectTriangle(Ray& ray, Vec3f& v0, Vec3f& v1, Vec3f& v2, float& t)
{
    float matrixA[3][3] = {v0.x - v1.x, v0.x - v2.x, ray.dir.x,
                           v0.y - v1.y, v0.y - v2.y, ray.dir.y,
                           v0.z - v1.z, v0.z - v2.z, ray.dir.z};
    float detA = determinant(matrixA);
    if (detA == 0)
        return false;

    // Cramers Rule
    float matrixBeta[3][3] = {v0.x - ray.origin.x, v0.x - v2.x, ray.dir.x,
                              v0.y - ray.origin.y, v0.y - v2.y, ray.dir.y,
                              v0.z - ray.origin.z, v0.z - v2.z, ray.dir.z};
    float beta = determinant(matrixBeta) / detA;
    if (beta < 0)
        return false;

    float matrixGama[3][3] = {v0.x - v1.x, v0.x - ray.origin.x, ray.dir.x,
                              v0.y - v1.y, v0.y - ray.origin.y, ray.dir.y,
                              v0.z - v1.z, v0.z - ray.origin.z, ray.dir.z};
    float gama = determinant(matrixGama) / detA;
    if (gama < 0 || gama + beta > 1)
        return false;

    float matrixT[3][3] = {v0.x - v1.x, v0.x - v2.x, v0.x - ray.origin.x,
                           v0.y - v1.y, v0.y - v2.y, v0.y - ray.origin.y,
                           v0.z - v1.z, v0.z - v2.z, v0.z - ray.origin.z};

    // t is out parameter.
    t = determinant(matrixT) / detA;
    return t > 0.0f;
}