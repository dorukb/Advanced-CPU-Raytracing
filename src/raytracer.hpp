#ifndef __HW1_RAYTRACER__
#define __HW1_RAYTRACER__

#include "camera.hpp"
#include "mesh.hpp"
#include "scene.h"
#include <random>

namespace DorkTracer
{
    class Raytracer
    {

    public:
        Scene scene;
        Raytracer(Scene& scene);
        Vec3f RenderPixel(int i, int j, Camera& cam);

    private:
        
        std::mt19937 roughnessRandomGenerator;
        std::uniform_real_distribution<> roughnessRandomDistro;

        std::mt19937 dofLensSampleGenerator;
        std::uniform_real_distribution<> dofLensSampleDistro;

        Vec3f PerPixel(int i, int j, Camera& cam);
        Ray GenerateRay(int i, int j, Camera& cam);
        Ray GenerateSecondaryRay(Ray& original, Vec3f& reflectDir, Vec3f& newOrigin);

        void IntersectObjects(Ray& ray);

        Vec3f GetAmbient(Vec3f& reflectance, Vec3f& ambientLightColor);
        Vec3f GetDiffuse(Shape* shape, Vec3f& reflectance, Vec3f& w_i, Ray& ray, Vec3f& receivedIrradiance);
        Vec3f GetSpecular(Vec3f& reflectance, float phongExp, Vec3f& w_in, Vec3f& w_out, Vec3f& normal, Vec3f& receivedIrradiance);
        bool IsInShadow(Ray& ray, Vec3f& lightPos);

        Vec3f PerformShading(Ray& ray, Vec3f& eyePos, int recursionDepth);
        Vec3f ComputeMirrorReflection(Ray& ray, Material& mat, Vec3f& w_o, int recursionDepth);
        Vec3f ComputeDielectricFresnelReflectionAndRefraction(Ray& ray, Material& mat, Vec3f& w_o, float n1, float n2, int recDepth);
        Vec3f ComputeConductorFresnelReflection(Ray& ray, Material& mat, Vec3f& w_o, int recDepth);
        Vec3f Reflect(Vec3f& normal, Vec3f& w_o, float roughness);
        Vec3f BeersLaw(float x, Vec3f c, Vec3f L_0);
        
        float GetRandom();
        float GetLensSample();
    };
}

#endif