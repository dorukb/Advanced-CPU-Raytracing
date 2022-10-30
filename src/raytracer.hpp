#ifndef __HW1_RAYTRACER__
#define __HW1_RAYTRACER__

#pragma once

#include "camera.hpp"
#include "parser.h"
using namespace parser;

namespace DorkTracer
{
    struct HitInfo{
        bool hasHit;
        float minT;
        Vec3f normal;
        int matId;
    };

    struct Ray{
        Vec3f origin;
        Vec3f dir;
        HitInfo hitInfo;
        float refractiveIndexOfCurrentMedium;
    };

    class Raytracer
    {

    public:
        parser::Scene scene;
        Raytracer(parser::Scene& scene);
        Vec3i RenderPixel(int i, int j, int camIndex);

    private:
        Vec3i PerPixel(int i, int j, Camera& cam);
        Ray GenerateRay(int i, int j, Camera& cam);

        void IntersectFace(Ray& ray, Face& face, std::vector<Vec3f>& vertices, int matId);
        bool DoesIntersectTriangle(Ray& ray, Vec3f& v0, Vec3f& v1, Vec3f& v2, float& t);
        void IntersectSphere(Ray& r, Sphere& s);
        void IntersectObjects(Ray& ray);

        Vec3f GetAmbient(Vec3f& reflectance, Vec3f& ambientLightColor);
        Vec3f GetDiffuse(Vec3f& reflectance, Vec3f& w_i, Vec3f& normal, Vec3f& receivedIrradiance);
        Vec3f GetSpecular(Vec3f& reflectance, float phongExp, Vec3f& w_in, Vec3f& w_out, Vec3f& normal, Vec3f& receivedIrradiance);
        bool IsInShadow(Vec3f& intersectionPoint, Vec3f& lightPos);
        bool IsBackface(Face& face, Vec3f& rayDir);
        Vec3f PerformShading(Ray& ray, Vec3f& eyePos, int recursionDepth);
        Vec3f ComputeMirrorReflection(Vec3f reflectance, Vec3f& w_o, Vec3f& normal, Vec3f& intersectionPoint, int recursionDepth);
        Vec3f ComputeDielectricFresnelReflectionAndRefraction(Material& mat, Vec3f x, Vec3f& w_o, Vec3f normal, float n1, float n2, int recDepth);
        Vec3f Reflect(Vec3f& normal, Vec3f& w_o);
        Vec3f BeersLaw(float x, Vec3f c, Vec3f L_0);
        Vec3f ComputeConductorFresnelReflection(Vec3f& w_o, Vec3f& n, Vec3f intPoint, Material& conductorMat, int recDepth);

    };
}

#endif