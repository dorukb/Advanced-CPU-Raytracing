#include "mesh.hpp"
#include <iostream>

using namespace DorkTracer;

DorkTracer::Mesh::Mesh(std::vector<Vec3f>& vertices)
{
    this->vertices = vertices;
}
int DorkTracer::Mesh::GetMaterial(){
    return this->material_id;
}
void DorkTracer::Mesh::SetMaterial(int matId){
    this->material_id = matId;
}

void DorkTracer::Mesh::ConstructBVH()
{
    // step 1. compute bbox for all objects contained in this mesh. <-- done at parsing.

    // init bvh tree
    int n = this->faces.size();
    this->bvh.resize(n*2-1);

    // create root node.
    BVH& root  = this->bvh[0];
    root.left = root.right = nullptr;
    root.bbox = this->bbox;
    root.firstFace = 0;
    root.faceCount = this->faces.size();

    this->nextFreeNodeIdx = 1;
    
    // -> can choose Middle for simplicity,
    // -> OR sort the objects and choose the median object’s centroid,
    // -> OR Surface Area Heuristic: choose the point such that the total surface areas in two parts are about the same.
    
    // Step 4. assign objects to each part
    // -> use the centroid for decision criteria
    // -> Can be done by swapping the elements in the container (e.g. vector)
    
    RecursiveBVHBuild(0);
}   

void DorkTracer::Mesh::RecursiveBVHBuild(uint32_t nodeIdx)
{
    BVH& node = this->bvh[nodeIdx];

    //base condition
    if(node.faceCount < 2) return;

    // step 2. select a splitting axis.  alternate axes: do: x,y,z,x,y,z,x . . .
    // step 3. select a splitting Point on that axis.
    float lenX = node.bbox.maxCorner.x - node.bbox.minCorner.x;
    float lenY = node.bbox.maxCorner.y - node.bbox.minCorner.y;
    float lenZ = node.bbox.maxCorner.z - node.bbox.minCorner.z;
    float splitPosition;
    uint32_t axisNum = 0;
    if(lenX > lenY){
        if(lenX > lenZ)
        {
            // x greatest
            splitPosition = node.bbox.minCorner.x + lenX * 0.5f;
            axisNum = 0;
        }
        else{
            // z greatest
            splitPosition = node.bbox.minCorner.z + lenZ * 0.5f;
            axisNum = 2;
        }
    }
    else{
        if(lenY > lenZ){
            // y greatest
            splitPosition = node.bbox.minCorner.y + lenY * 0.5f;
            axisNum = 1;
        }
        else{
            // z greatest
            splitPosition = node.bbox.minCorner.z + lenZ * 0.5f;
            axisNum = 2;
        }
    }

    // partition into two pieces, similar to quicksort.
    int i = node.firstFace;
    int j = i + node.faceCount-1;
    while(i <= j){
        if(this->faces[i].center[axisNum] < splitPosition){
            i++;
        }
        else{
            std::swap(this->faces[i], this->faces[j]);
            j--;
        }
    }
    int leftCount = i - node.firstFace;

    bool isOneHalfEmpty = leftCount == 0 || leftCount == node.faceCount;
    if(isOneHalfEmpty) return;
    
    // Set child nodes.
    int leftIndex = this->nextFreeNodeIdx;
    this->nextFreeNodeIdx++;

    int rightIdx = this->nextFreeNodeIdx;
    this->nextFreeNodeIdx++;

    this->bvh[leftIndex].firstFace = node.firstFace;
    this->bvh[leftIndex].faceCount = leftCount;

    node.left = &(this->bvh[leftIndex]);
    
    this->bvh[rightIdx].firstFace = i;
    this->bvh[rightIdx].faceCount = node.faceCount - leftCount;
    node.right = &(this->bvh[rightIdx]);

    // mark current node as "interior", it has no faces/geometry.
    node.faceCount = 0;


    // Update the bounding boxes of the left & right pieces.
    RecomputeBoundingBox(leftIndex);
    RecomputeBoundingBox(rightIdx);

    // continue.
    RecursiveBVHBuild(leftIndex);
    RecursiveBVHBuild(rightIdx);
}

void DorkTracer::Mesh::RecomputeBoundingBox(uint32_t nodeIdx)
{
    BVH& node = this->bvh[nodeIdx];
    node.bbox.minCorner.x = node.bbox.minCorner.y = node.bbox.minCorner.z = INFINITY;
    node.bbox.maxCorner.x = node.bbox.maxCorner.y = node.bbox.maxCorner.z = -INFINITY;

    uint32_t first = node.firstFace;
    for (uint32_t i = 0; i < node.faceCount; i++)
    {
        Face& face = this->faces[first + i];

        node.bbox.minCorner.x = std::min(node.bbox.minCorner.x, face.bbox.minCorner.x);
        node.bbox.minCorner.y = std::min(node.bbox.minCorner.y, face.bbox.minCorner.y);
        node.bbox.minCorner.z = std::min(node.bbox.minCorner.z, face.bbox.minCorner.z);

        node.bbox.maxCorner.x = std::max(node.bbox.maxCorner.x, face.bbox.maxCorner.x);
        node.bbox.maxCorner.y = std::max(node.bbox.maxCorner.y, face.bbox.maxCorner.y);
        node.bbox.maxCorner.z = std::max(node.bbox.maxCorner.z, face.bbox.maxCorner.z);
    }
}

bool DorkTracer::Mesh::Intersect(Ray& ray)
{    
    // Transform the ray into our local space.

    // TODO: improve this. this is inefficient.
    Vec3f rayOriginCache = ray.origin;
    Vec3f rayDirCache = ray.dir;

    Vec4f rayOrigin(ray.origin, 1.0f);
    Vec4f rayDir(ray.dir, 0.0f);
    
    ray.origin =  Matrix::ApplyTransform(this->inverseTransform, rayOrigin);

    Vec3f transformedDir = Matrix::ApplyTransform(this->inverseTransform, rayDir);
    ray.dir = makeUnit(transformedDir);
    // note that BVH intersection test actually hinders BFC perf gain.
    // BFC actually slows down the computation.
    // observed 15.5s NO bfc, 16.7s YES bfc.
    if(this->bbox.doesIntersectWith(ray)){
        bool hasHit = this->bvh[0].IntersectBVH(ray, *(this));
        // undo origin&dir changes to local space.
        ray.origin = rayOriginCache;
        ray.dir = rayDirCache;
        if(hasHit){
            ray.hitInfo.normal = makeUnit(Matrix::ApplyTransform(this->inverseTransposeTransform, Vec4f(ray.hitInfo.normal, 0.0f)));
        }

        // return ray.hitInfo.hasHit;
        return hasHit;
    }
    else{
        ray.origin = rayOriginCache;
        ray.dir = rayDirCache;
        return false;
    } 
}

bool DorkTracer::Mesh::IntersectFace(Ray& ray, uint32_t faceIdx){
    return IntersectFace(ray, this->faces[faceIdx]);
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
        // ray.hitInfo.normal = makeUnit(Matrix::ApplyTransform(this->inverseTransposeTransform, Vec4f(face.n, 0.0f)));

        ray.hitInfo.matId = this->GetMaterial();
        return true;
    }
    else return false;
    // return hasIntersected;
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