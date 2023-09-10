#ifndef __HW1_RAYTRACER__
#define __HW1_RAYTRACER__

#include "camera.hpp"
#include "mesh.hpp"
#include "scene.h"
#include <random>
#include "rendererParams.h"

namespace DorkTracer
{
    class Raytracer
    {

    public:
        Scene scene;
        Camera* activeCamera;
        Raytracer(Scene& scene);
        Vec3f RenderPixel(int i, int j, Camera& cam);
        void EnablePathTracing(RendererParams params);

    private:
        bool pathTracingEnabled = false;
        RendererParams rendererParams;

        std::mt19937 randGen;
        std::uniform_real_distribution<> roughnessRandomDistro;
        std::uniform_real_distribution<> normalizedDistro;

        std::mt19937 dofLensSampleGenerator;
        std::uniform_real_distribution<> dofLensSampleDistro;

        Vec3f PerPixel(int i, int j, Camera& cam);
        Ray GenerateRay(int i, int j, Camera& cam);
        Ray GenerateSecondaryRay(Ray& original, Vec3f& reflectDir, Vec3f& newOrigin);

        void IntersectObjects(Ray& ray);

        Vec3f GetAmbient(Vec3f& reflectance, Vec3f& ambientLightColor);
        Vec3f GetDiffuse(Ray& ray, Shape* shape, Material& mat, Vec3f& w_i, Vec3f& receivedIrradiance);
        Vec3f GetSpecular(Ray& ray, Shape* shape, Material& mat, Vec3f& w_in, Vec3f& w_out, Vec3f& receivedIrradiance);
        Vec3f SampleDirectLighting(Ray& ray, Material& mat, Vec3f& w_o, int lightIDToSkip);

        bool IsInShadow(Ray& ray, Vec3f& lightPos);
        bool IsInShadowDirectional(Ray& originalRay, Vec3f& dir);
        bool CastShadowRay(Ray& shadowRay, float lightSourceT);

        Vec3f PerformShading(Ray& ray, Vec3f& eyePos, int recursionDepth);
        Vec3f ComputeMirrorReflection(Ray& ray, Material& mat, Vec3f& w_o, int recursionDepth);
        Vec3f ComputeDielectricFresnelReflectionAndRefraction(Ray& ray, Material& mat, Vec3f& w_o, float n1, float n2, int recDepth);
        Vec3f ComputeConductorFresnelReflection(Ray& ray, Material& mat, Vec3f& w_o, int recDepth);
        Vec3f Reflect(Vec3f& normal, Vec3f& w_o, float roughness);
        Vec3f BeersLaw(float x, Vec3f c, Vec3f L_0);

        Vec3f GetDiffuseReflectanceCoeff(Ray& ray, Shape* shape, Material& mat);
        Vec3f GetSpecularReflectanceCoeff(Ray& ray, Shape* shape, Material& mat);
        Vec3f ComputeGlobalIllumination(Ray& originalRay,Material& orgMat, Vec3f& w_o, int recDepth, int& hitLightId);
        Vec3f Shade(Ray& ray, Material& mat, Vec3f& w_i, Vec3f& w_o, Vec3f& receivedIrradiance);

        float GetNormalizedRandom();
        float GetLensSample();
    };
}

#endif