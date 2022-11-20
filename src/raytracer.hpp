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
        
        std::mt19937 gRandomGenerator;
        std::uniform_real_distribution<> roughnessRandomDistro01;

        Vec3f PerPixel(int i, int j, Camera& cam);
        Ray GenerateRay(int i, int j, Camera& cam);

        void IntersectObjects(Ray& ray);

        Vec3f GetAmbient(Vec3f& reflectance, Vec3f& ambientLightColor);
        Vec3f GetDiffuse(Vec3f& reflectance, Vec3f& w_i, Vec3f& normal, Vec3f& receivedIrradiance);
        Vec3f GetSpecular(Vec3f& reflectance, float phongExp, Vec3f& w_in, Vec3f& w_out, Vec3f& normal, Vec3f& receivedIrradiance);
        bool IsInShadow(Vec3f& intersectionPoint, Vec3f& lightPos);

        Vec3f PerformShading(Ray& ray, Vec3f& eyePos, int recursionDepth);
        Vec3f ComputeMirrorReflection(Material& mat, Vec3f& w_o, Vec3f& normal, Vec3f& intersectionPoint, int recursionDepth);
        Vec3f ComputeDielectricFresnelReflectionAndRefraction(Material& mat, Vec3f x, Vec3f& w_o, Vec3f normal, float n1, float n2, int recDepth);
        Vec3f Reflect(Vec3f& normal, Vec3f& w_o, float roughness);
        Vec3f BeersLaw(float x, Vec3f c, Vec3f L_0);
        Vec3f ComputeConductorFresnelReflection(Vec3f& w_o, Vec3f& n, Vec3f intPoint, Material& conductorMat, int recDepth);
        float GetRandom01();
    };
}

#endif