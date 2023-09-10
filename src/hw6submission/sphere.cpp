#include "sphere.hpp"
using namespace DorkTracer;

Sphere::Sphere(std::vector<Vec3f>& vertex_data)
{
    this->vertex_data = vertex_data;
}

float inline MakeGreyscale(Vec3f& color01){
    return color01.x + color01.y + color01.z; // maybe divide by 3 to normalize?
}

bool DorkTracer::Sphere::Intersect(Ray& r)
{

    Vec3f center = this->vertex_data[center_vertex_id - 1];

    Vec3f rayOriginCache = r.origin;
    Vec3f rayDirCache = r.dir;
    
    
    // Transform the ray into our local space.
    r.origin = Matrix::ApplyTransformToPoint(this->inverseTransform, r.origin);
    r.dir = Matrix::ApplyTransformToVector(this->inverseTransform, r.dir);

    if(this->hasMotionBlur)
    {
        r.origin = r.origin + motionBlurVector * r.motionBlurTime;
    }
    
    Vec3f oc = r.origin - center;
    float t, t1, t2;

    float c = dot(oc,oc) - (radius * radius);
    float b =  2 * dot(r.dir, oc);
    float a = dot(r.dir, r.dir);
    float delta = b*b - (4 * a * c);

    if(delta < 0.0f){
        r.origin = rayOriginCache;
        r.dir = rayDirCache;
        return false;
    }
    else
    {
        delta = sqrtf(delta);
        a = 2.0 * a;
        float t1 = (-b + delta) / a;
        float t2 = (-b - delta) / a;
        t = t1 < t2 ? t1 : t2;

        // Find the minimum t which is > 0.0f
        if(t1 < t2){
            if(t1 > 0.0f){
                t = t1;
            }
            else t = t2;
        }
        else if(t2 < t1){
            if(t2 > 0.0f){
                t = t2;
            }
            else t = t1;
        }

        Vec3f localhitPoint = r.origin + r.dir * t;

        // undo ray transformation
        r.origin = rayOriginCache;
        r.dir = rayDirCache;

        if(t < r.hitInfo.minT && t > 0.0f)
        {
            r.hitInfo.minT = t;
            r.hitInfo.matId = material_id;
            r.hitInfo.hitShape = this;
            r.hitInfo.hasHit = true;

            // calculate sphere normal at hit point
            Vec3f worldHitPoint = r.origin + r.dir * t;
            r.hitInfo.hitPoint = worldHitPoint;

            // compute hit point UV.
            Vec3f p = localhitPoint - center;
            float phi = std::atan2(p.z, p.x);
            float theta = std::acos(p.y / radius);

            float u = (-phi + M_PI) / (2.0f * M_PI);
            float v = theta / M_PI;

            r.hitInfo.hitUV.x = u;
            r.hitInfo.hitUV.y = v;


            if(HasNormalMap())
            {
                Vec3f sampledNormal = this->normalMap->GetRGBSample(u,v);
                sampledNormal = sampledNormal / (127.5f) - Vec3f(1,1,1);
                sampledNormal = makeUnit(sampledNormal);

                // convert normal map vector from canonical to local space.

                Vec3f tan,bitan;
                GetTangentAndBitangentAroundPoint(p, radius, phi, theta, tan, bitan);
                Vec3f normal = makeUnit(cross(bitan, tan));

                r.hitInfo.normal = GetTransformedNormal(tan, bitan, normal, sampledNormal);
            }
            else if(HasBumpMap())
            {
                Vec3f tan,bitan;
                GetTangentAndBitangentAroundPoint(p, radius, phi, theta, tan, bitan);
                // tan = P_u (dP/du), bitan = P_v (dP/dv)
                Vec3f N = makeUnit(cross(bitan, tan));

                if(this->bumpMap->IsGenerated())
                {
                    Vec3f gradient;
                    float eps = 0.001;
                    float hxyz = this->bumpMap->GetSampleFromWorldPos(p.x, p.y, p.z);
                    gradient.x = (this->bumpMap->GetSampleFromWorldPos(p.x + eps, p.y, p.z) - hxyz) / eps;
                    gradient.y = (this->bumpMap->GetSampleFromWorldPos(p.x, p.y + eps, p.z) - hxyz) / eps;
                    gradient.z = (this->bumpMap->GetSampleFromWorldPos(p.x, p.y, p.z + eps) - hxyz) / eps;

                    Vec3f gParallel = N * dot(gradient, N);
                    Vec3f surfaceGradient = gradient - gParallel;

                    Vec3f newNormal = N - surfaceGradient;
                    r.hitInfo.normal = makeUnit(newNormal);
                }
                else{
                    // Bumped surface equation:
                    // q(u,v) = p(u,v) + h(u,v) * n(u,v) then normal: n' = (dq/dv) x (dq/du)
                    // now calculate the partial derivatives.
                    // q_u = P_u + h_u * n(u,v) + ignoredTerm. P_u = tan, n(u,v) = N.
                    float width = this->bumpMap->GetWidth();
                    float height = this->bumpMap->GetHeight();

                    int i = (int)(u * width);
                    int j = (int)(v * height);

                    float normalizer = this->bumpMap->GetNormalizer();
                    float bumpFactor = this->bumpMap->GetSampleMultiplier();

                    // Calculate h_u using forward differences method.
                    Vec3f sampleColor = this->bumpMap->GetDirectSample(i+1, j) / normalizer;
                    float h1 = MakeGreyscale(sampleColor) * bumpFactor;
                    
                    sampleColor = this->bumpMap->GetDirectSample(i,j) / normalizer;
                    float h_uv = MakeGreyscale(sampleColor) * bumpFactor;

                    sampleColor = this->bumpMap->GetDirectSample(i, j+1) / normalizer;
                    float h2 = MakeGreyscale(sampleColor) * bumpFactor;

                    Vec3f q_u = tan + N * (h1-h_uv);
                    Vec3f q_v = bitan + N * (h2-h_uv);

                    Vec3f newNormal = cross(q_v, q_u);
                    r.hitInfo.normal = makeUnit(newNormal);
                }
                
            }
            else
            {
                r.hitInfo.normal = makeUnit(localhitPoint - center);
            }

            r.hitInfo.normal = makeUnit(Matrix::ApplyTransformToVector(inverseTransposeTransform, r.hitInfo.normal));
            return true;
        }
        else return false;
    }
}
void DorkTracer::Sphere::GetTangentAndBitangentAroundPoint(Vec3f& p, float radius, float phi, float theta, Vec3f& tan, Vec3f& bitan)
{
    tan.x = 2 * M_PI * p.z;
    tan.y = 0;
    tan.z = -2 * M_PI * p.x;

    bitan.x = M_PI * p.y * std::cos(phi);
    bitan.y = -radius * M_PI * std::sin(theta);
    bitan.z = M_PI * p.y * std::sin(phi);

    tan = makeUnit(tan);
    bitan = makeUnit(bitan);
}