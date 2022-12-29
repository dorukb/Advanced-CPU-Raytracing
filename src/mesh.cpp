#include "mesh.hpp"
#include <iostream>
#include <random>

using namespace DorkTracer;

DorkTracer::Mesh::Mesh(std::vector<Vec3f>& vertices, std::vector<Vec2f>& uv)
{
    this->vertices = vertices;
    this->uv = uv;
    this->vertexOffset = 0;
    this->textureOffset = 0;
}
int DorkTracer::Mesh::GetMaterial(){
    return this->material_id;
}
void DorkTracer::Mesh::SetMaterial(int matId){
    this->material_id = matId;
}

Vec2f& DorkTracer::Mesh::GetUv(int idx){
    return this->uv[idx - 1 + textureOffset];
}
Vec3f& DorkTracer::Mesh::GetVertex(int idx){
    return this->vertices[idx - 1 + vertexOffset];
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
    Vec3f rayOriginCache = ray.origin;
    Vec3f rayDirCache = ray.dir;
    
    ray.origin = Matrix::ApplyTransformToPoint(this->inverseTransform, ray.origin);
    ray.dir = Matrix::ApplyTransformToVector(this->inverseTransform, ray.dir);

    if(this->hasMotionBlur)
    {
        ray.origin = ray.origin + motionBlurVector * ray.motionBlurTime;
    }

    if(this->bbox.doesIntersectWith(ray)){
        bool hasHit = this->bvh[0].IntersectBVH(ray, *(this));
        // undo origin&dir changes to local space.
        ray.origin = rayOriginCache;
        ray.dir = rayDirCache;
        if(hasHit){
            ray.hitInfo.hitPoint = ray.origin + ray.dir *ray.hitInfo.minT;
            ray.hitInfo.normal = makeUnit(Matrix::ApplyTransformToVector(this->inverseTransposeTransform, ray.hitInfo.normal));
        }
        return hasHit;
    }
    else{
        ray.origin = rayOriginCache;
        ray.dir = rayDirCache;
        return false;
    } 
}

void DorkTracer::Mesh::SetAccessOffsets(int vertexOffset, int textureOffset)
{
    this->vertexOffset = vertexOffset;
    this->textureOffset = textureOffset;
}
float inline MakeGreyscale(Vec3f& color){
    return (color.x + color.y + color.z) / 3.0f;// maybe divide by 3 to normalize?
}
bool DorkTracer::Mesh::IntersectFace(Ray& ray, uint32_t faceIdx){
    return IntersectFace(ray, this->faces[faceIdx]);
}
bool DorkTracer::Mesh::IntersectFace(Ray& ray, Face& face)
{
    Vec3f& v0 = GetVertex(face.v0_id);
    Vec3f& v1 = GetVertex(face.v1_id);
    Vec3f& v2 = GetVertex(face.v2_id);
    float newT  = INFINITY;

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
    float t = determinant(matrixT) / detA;
    bool isClosestIntersection = t > 0.0f && t < ray.hitInfo.minT;
    if(isClosestIntersection)
    {
        ray.hitInfo.minT = t;
        ray.hitInfo.hasHit = true;
        ray.hitInfo.normal = face.n;
        ray.hitInfo.hitPoint = ray.origin + ray.dir * ray.hitInfo.minT;

        // Calculate and store UV value for this hit point.
        if(uv.size() > 0)
        {
            Vec2f& v0_uv = GetUv(face.v0_id);
            Vec2f& v1_uv = GetUv(face.v1_id);
            Vec2f& v2_uv = GetUv(face.v2_id);

            // u(B,Y) = u_a + B(u_b - u_a) + Y(u_c - u_a)
            // v(B,Y) = v_a + B(v_b - v_a) + Y(v_c - u_a)
            float u = v0_uv.x + beta * (v1_uv.x - v0_uv.x) + gama * (v2_uv.x - v0_uv.x);
            float v = v0_uv.y + beta * (v1_uv.y - v0_uv.y) + gama * (v2_uv.y - v0_uv.y);

            // Support tiling.
            u = GetFloorForTiledUV(u);
            v = GetFloorForTiledUV(v);

            ray.hitInfo.hitUV.x = u;
            ray.hitInfo.hitUV.y = v;

            // Can have either normal or bump map, not both.
            if(HasNormalMap())
            {
                Vec3f sampledNormal = this->normalMap->GetRGBSample(u,v);
                sampledNormal = sampledNormal / (127.5f) - Vec3f(1,1,1);
                sampledNormal = makeUnit(sampledNormal);

                // convert sampled normal vector from canonical to local space.
                Vec3f tan,bitan;
                GetTangentAndBitangentForTriangle(v0, v1, v2, v0_uv, v1_uv, v2_uv, tan, bitan);
                ray.hitInfo.normal = GetTransformedNormal(tan, bitan, face.n, sampledNormal);
                ray.hitInfo.normal = makeUnit(Matrix::ApplyTransformToVector(inverseTransposeTransform, ray.hitInfo.normal));
            }
            else if(HasBumpMap())
            {
                Vec3f tan,bitan;

                // Vec3f v0t = Matrix::ApplyTransformToPoint(transform, v0);
                // Vec3f v1t = Matrix::ApplyTransformToPoint(transform, v1);
                // Vec3f v2t = Matrix::ApplyTransformToPoint(transform, v2);
                GetTangentAndBitangentForTriangle(v0, v1, v2, v0_uv, v1_uv, v2_uv, tan, bitan);

                // GetTangentAndBitangentForTriangle(v0t, v1t, v2t, v0_uv, v1_uv, v2_uv, tan, bitan);
               // tan = P_u (dP/du), bitan = P_v (dP/dv)
                // Vec3f N = makeUnit(cross(bitan, tan));
                Vec3f N = face.n;

                if(this->bumpMap->IsGenerated())
                {
                    Vec3f gradient;
                    float eps = 0.001;
                    Vec3f& p = ray.hitInfo.hitPoint;
                    float bf = this->bumpMap->GetSampleMultiplier();

                    float hxyz = this->bumpMap->GetSampleFromWorldPos(p.x, p.y, p.z) * bf;
                    gradient.x = (this->bumpMap->GetSampleFromWorldPos(p.x + eps, p.y, p.z) * bf - hxyz) / eps;
                    gradient.y = (this->bumpMap->GetSampleFromWorldPos(p.x, p.y + eps, p.z) * bf - hxyz) / eps;
                    gradient.z = (this->bumpMap->GetSampleFromWorldPos(p.x, p.y, p.z + eps) * bf - hxyz) / eps;

                    Vec3f gParallel = N * dot(gradient, N);
                    Vec3f surfaceGradient = gradient - gParallel;

                    Vec3f newNormal = N - surfaceGradient;
                    ray.hitInfo.normal = makeUnit(newNormal);
                    ray.hitInfo.normal = makeUnit(Matrix::ApplyTransformToVector(inverseTransposeTransform, ray.hitInfo.normal));

                }
                else
                {
                    // Bumped surface equation:
                    // q(u,v) = p(u,v) + h(u,v) * n(u,v) then normal: n' = (dq/dv) x (dq/du)
                    // now calculate the partial derivatives.
                    // q_u = P_u + h_u * n(u,v) + ignoredTerm. P_u = tan, n(u,v) = N.
                    float width = this->bumpMap->GetWidth();
                    float height = this->bumpMap->GetHeight();

                    int i = (int)(u * (width-1));
                    int j = (int)(v * (height-1));

                    // float normalizer = 1.0f;
                    // Calculate h_u using forward differences method.
                    int nextI = i+1;
                    int nextJ = j+1;
                    if(i == width-1) nextI = i;
                    if(j == height-1) nextJ = j;

                    Vec3f sampleColor = this->bumpMap->GetDirectSample(i,j);
                    float h_uv = MakeGreyscale(sampleColor);

                    sampleColor = this->bumpMap->GetDirectSample(nextI, j);
                    float hDeltaU = MakeGreyscale(sampleColor);
                
                    sampleColor = this->bumpMap->GetDirectSample(i, nextJ);
                    float hDeltaV = MakeGreyscale(sampleColor);

                    float bumpFactor = this->bumpMap->GetSampleMultiplier();
                    Vec3f q_u = tan + N * ((hDeltaU-h_uv) * bumpFactor);
                    Vec3f q_v = bitan + N * ((hDeltaV-h_uv)  * bumpFactor);

                    Vec3f newNormal = cross(q_v, q_u);
                    ray.hitInfo.normal = makeUnit(newNormal);
                    
                     if(newNormal.x * N.x <= 0 &&
                        newNormal.y * N.y <= 0 &&
                        newNormal.z * N.z <= 0  ){
                        ray.hitInfo.normal = ray.hitInfo.normal * -1;
                    }
                    else if(std::abs(newNormal.y - N.y) > 0.9f ||
                        std::abs(newNormal.x - N.x) > 0.9f ||
                        std::abs(newNormal.z - N.z) > 0.9f){
                        ray.hitInfo.normal = ray.hitInfo.normal * -1;
                    }
                    ray.hitInfo.normal = makeUnit(Matrix::ApplyTransformToVector(inverseTransposeTransform, ray.hitInfo.normal));

                }            
              
            }
            
        }
        else{
            ray.hitInfo.normal = makeUnit(Matrix::ApplyTransformToVector(inverseTransposeTransform, ray.hitInfo.normal));
        }


        ray.hitInfo.matId = this->GetMaterial();
        ray.hitInfo.hitShape = this;
        return true;
    }
    else return false;
}
bool DorkTracer::Mesh::IsBackface(Face& face, Vec3f& rayDir)
{
    return dot(face.n, rayDir) > 0.0f;
}

// bool DorkTracer::Mesh::DoesIntersectTriangle(Ray& ray, Vec3f& v0, Vec3f& v1, Vec3f& v2, float& t, float currMinT)
// {
 
// }
float inline DorkTracer::Mesh::GetFloorForTiledUV(float x)
{
    if(x > 1.0001f){
        x = x - std::floor(x);
        if(x < 0.0001) x = 1.0f;
    }
    return x;
}
void DorkTracer::Mesh::GetTangentAndBitangentForTriangle(Vec3f& vert0, Vec3f& vert1, Vec3f& vert2, Vec2f& v0_uv, Vec2f& v1_uv, Vec2f& v2_uv, Vec3f& tan, Vec3f& bitan)
{
    Vec3f e1 = makeUnit(vert1-vert0);
    Vec3f e2 = makeUnit(vert2-vert1);


    float v0u = GetFloorForTiledUV(v0_uv.x);
    float v0v = GetFloorForTiledUV(v0_uv.y);

    float v1u = GetFloorForTiledUV(v1_uv.x);
    float v1v = GetFloorForTiledUV(v1_uv.y);

    float v2u = GetFloorForTiledUV(v2_uv.x);
    float v2v = GetFloorForTiledUV(v2_uv.y);

    float u1 = v1u - v0u;
    float v1 = v1v - v0v;
    float u2 = v2u - v1u;
    float v2 = v2v - v1v;

    float det = 1.0f / (u1*v2 - v1*u2);

    tan.x = det * (v2*e1.x - v1*e2.x);
    tan.y = det * (v2*e1.y - v1*e2.y);
    tan.z = det * (v2*e1.z - v1*e2.z);

    bitan.x = -det*u2*e1.x + det*u1*e2.x;
    bitan.y = -det*u2*e1.y + det*u1*e2.y;
    bitan.z = -det*u2*e1.z + det*u1*e2.z;

    tan = makeUnit(tan);
    bitan = makeUnit(bitan);
}