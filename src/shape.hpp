#ifndef __DORKTRACER_SHAPE__
#define __DORKTRACER_SHAPE__

#include <vector>
#include <random>

#include "ray.hpp"
#include "material.hpp"
#include "matrix.hpp"
#include "texture.h"

namespace DorkTracer{

    class Shape
    {
    protected:
        int material_id;

    public:
        bool isInstance = false;
        int id;
        Matrix transform;
        Matrix inverseTransform;
        Matrix inverseTransposeTransform;

        bool hasMotionBlur = false;
        Vec3f motionBlurVector;
        std::mt19937 motionBlurRandomGenerator;
        std::uniform_real_distribution<> motionBlurRandomDistro01;

        Texture* diffuseTex;
        Texture* specularTex;
        Texture* normalMap;
        Texture* bumpMap;
        Texture* replaceAll;
        
        int GetMaterial();
        void SetMaterial(int matId);

        Shape() : transform(4,4), inverseTransform(4,4), inverseTransposeTransform(4,4)
        {
            this->motionBlurVector = Vec3f();
            this->motionBlurRandomGenerator = std::mt19937();
            this->motionBlurRandomDistro01 = std::uniform_real_distribution<>(0.0f, 1.0f);
            this->diffuseTex = specularTex = replaceAll = normalMap = bumpMap = nullptr;
        }

        virtual bool Intersect(Ray&){return false;};
        
        float GetRandom01MotionBlurTime()
        {
            return this->motionBlurRandomDistro01(this->motionBlurRandomGenerator);
        }
        bool HasNormalMap(){
            return normalMap != nullptr;
        }
        bool HasBumpMap(){
            return bumpMap != nullptr;
        }
        bool HasDiffuseTexture(){
            return diffuseTex != nullptr;
        }      
        bool HasSpecularTexture(){
            return specularTex != nullptr;
        }
        bool HasReplaceAllTexture(){
            return replaceAll != nullptr;
        }
    };

    struct BoundingBox
    {
        Vec3f minCorner;
        Vec3f maxCorner;

        // bool doesIntersectWith(Ray& ray){return false;}

        bool doesIntersectWith(Ray& ray)
        {
            float tx1 = (minCorner.x - ray.origin.x) / ray.dir.x; 
            float tx2 = (maxCorner.x - ray.origin.x) / ray.dir.x;

            float tmin = tx1;
            float tmax = tx2;
            if(tx1 > tx2){
                tmin = tx2;
                tmax = tx1;
            }

            float ty1 = (minCorner.y - ray.origin.y) / ray.dir.y;
            float ty2 = (maxCorner.y - ray.origin.y) / ray.dir.y;
            tmin = std::fmax(tmin, std::fmin(ty1, ty2));
            tmax = std::fmin(tmax, std::fmax(ty1, ty2));

            float tz1 = (minCorner.z - ray.origin.z) / ray.dir.z;
            float tz2 = (maxCorner.z - ray.origin.z) / ray.dir.z;
            tmin = std::fmax(tmin, std::fmin(tz1, tz2));
            tmax = std::fmin(tmax, std::fmax(tz1, tz2));
            return tmax > 0 && tmax >= tmin && tmin < ray.hitInfo.minT;
        }
    };

    struct Face
    {
        int v0_id;
        int v1_id;
        int v2_id;
        Vec3f n;
        Vec3f center;
        BoundingBox bbox;
        double area;
    };


}

#endif