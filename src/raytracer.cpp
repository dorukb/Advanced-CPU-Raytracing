#include "raytracer.hpp"
#include <iostream>
#include <random>

using namespace DorkTracer;

Raytracer::Raytracer(Scene& scene){
    this->gRandomGenerator = std::mt19937();
    this->roughnessRandomDistro01 =  std::uniform_real_distribution<>(-0.5f, 0.5);
    this->scene = scene;
}

float Raytracer::GetRandom01()
{
    return this->roughnessRandomDistro01(this->gRandomGenerator);
}
Vec3f Raytracer::RenderPixel(int i, int j, Camera& cam)
{
    return PerPixel(i,j, cam);
}

Vec3f Raytracer::PerPixel(int coordX, int coordY, Camera& cam)
{
    Ray ray = GenerateRay(coordX, coordY, cam);
    
    Vec3f bgColor = {scene.background_color.x, scene.background_color.y, scene.background_color.z};
    // Intersect with all Objects in the scene to find the closest intersection point, if any.
    IntersectObjects(ray);
    
    if(ray.hitInfo.hasHit)
    {
       return PerformShading(ray, cam.position, scene.max_recursion_depth);
    }
    else return bgColor;
}

Vec3f Raytracer::PerformShading(Ray& ray, Vec3f& eyePos, int recursionDepth)
{
    Vec3f intersectionPoint = ray.origin + ray.dir * ray.hitInfo.minT;
    // Actually perform shading, depending on Material type.

    Vec3f color{0,0,0};
    Material& mat = scene.materials[ray.hitInfo.matId-1];
    Vec3f w_o = makeUnit(eyePos - intersectionPoint);
    float refractiveIndexOfVacuum = 1.00001;
    bool travellingInsideAnObject = ray.refractiveIndexOfCurrentMedium > refractiveIndexOfVacuum;

    if(!travellingInsideAnObject)
    {
        color = color + GetAmbient(mat.ambient, scene.ambient_light);
        // Repeat for each light source.
        for(int l = 0; l < scene.point_lights.size(); l++)
        {
            // Skip all calculation if in shadow, no color "contribution".
            PointLight& light = scene.point_lights[l];
            if(IsInShadow(intersectionPoint, light.position)){
                continue;
            }

            Vec3f w_in = makeUnit(light.position - intersectionPoint);
            
            float distToLight = len(light.position - intersectionPoint);
            Vec3f receivedIrradiance = light.intensity / (distToLight * distToLight);

            color = color + GetDiffuse(mat.diffuse, w_in, ray.hitInfo.normal, receivedIrradiance);
            color = color + GetSpecular(mat.specular, mat.phong_exponent, w_in, w_o, ray.hitInfo.normal, receivedIrradiance);
        }    
    }
 

    if (mat.type == Material::Default)
    {   // we are done with shading, default/standard only have Diffuse Specular and Ambient components.
        return color; 
    }
    else if(mat.type == Material::Mirror)
    {
        // Compute radiance along the ideal reflection ray
        color = color + ComputeMirrorReflection(mat, w_o, ray.hitInfo.normal, intersectionPoint, recursionDepth);
    }
    else if(mat.type == Material::Dielectric)
    {
        // Both reflection & transmission
        float n1 = ray.refractiveIndexOfCurrentMedium;
        color = color + ComputeDielectricFresnelReflectionAndRefraction(mat, intersectionPoint, w_o, ray.hitInfo.normal, n1, mat.refractiveIndex, recursionDepth);
    }
    else if(mat.type == Material::Conductor)
    {
        // Reflection & absorption, no tranmission.
        color = color + ComputeConductorFresnelReflection(w_o, ray.hitInfo.normal, intersectionPoint, mat, recursionDepth);
    }
    // Vec3i c = clamp(color);
    // color.x = c.x;
    // color.y = c.y;
    // color.z = c.z;
    
    return color;
}

Vec3f Raytracer::ComputeConductorFresnelReflection(Vec3f& w_o, Vec3f& n, Vec3f intPoint, Material& conductorMat, int recDepth)
{
    if(recDepth == 0){
        return Vec3f{0,0,0}; 
    }

    Vec3f d = -w_o;
    float cosTheta = -dot(d, n);
    float n2 = conductorMat.refractiveIndex;
    float k2 = conductorMat.conductorAbsorptionIndex;

    // compute common intermediate terms
    float n2k2 = n2*n2 + k2*k2;
    float n2cosTheta2 = 2 * n2 * cosTheta;
    float cosThetaSqr = cosTheta * cosTheta;

    float rs = (n2k2 - n2cosTheta2 + cosThetaSqr) / (n2k2 + n2cosTheta2 + cosThetaSqr);
    float rp = (n2k2 * cosThetaSqr - n2cosTheta2 + 1) / (n2k2 * cosThetaSqr + n2cosTheta2 + 1);

    float reflectRatio = 0.5 * (rs + rp);
    // transmit ratio = 0.

    // Cast reflected ray
    if(reflectRatio > 0.01)
    {
        Vec3f reflectedRaysColor;
        {
            Vec3f w_reflected = Reflect(n, w_o, conductorMat.roughness);

            Ray ray;
            ray.origin = intPoint + n * scene.shadow_ray_epsilon;
            ray.dir = makeUnit(w_reflected);
            ray.hitInfo.hasHit = false;
            ray.hitInfo.minT = INFINITY;
            ray.refractiveIndexOfCurrentMedium = 1.0f;

            IntersectObjects(ray);
            if(ray.hitInfo.hasHit)
            {
                // no attenuation on reflect case, since we are not INSIDE an object?
                reflectedRaysColor = conductorMat.mirror * PerformShading(ray, ray.origin, recDepth-1);
            }
            else reflectedRaysColor = Vec3f{0,0,0};
        }

        return reflectedRaysColor * reflectRatio;
    }
    else{
        return Vec3f{0,0,0};
    }
 
}
// Ratio of reflected over refracted light amount, for dielectric - dielectric interfacing.
// cosTheta = dot(w_i, normal)
// n1: refractive index of the medium incoming ray is currently in
// n2: refractive index of the medium ray is entering, or equivalently that it will be transmitted in.
// float GetFresnelReflectionRatioDielectric(Vec3f& w_i, Vec3f& normal, float n1, float n2)

Vec3f Raytracer::ComputeDielectricFresnelReflectionAndRefraction(Material& mat, Vec3f x, Vec3f& w_o, Vec3f normal, float n1, float n2, int recDepth)
{
    if(recDepth == 0){
        return Vec3f{0,0,0}; 
    }

    Vec3f d = -w_o;
    float cosTheta = -dot(d, normal);
    bool isEntering = cosTheta > 0.f;
    float raysN = n1;
    float objN = n2;

    if(!isEntering){
        // std::cout << "was not entering??"<<std::endl;
        // swap indices for correctness.
        float tmp = n1;
        n1 = n2;
        n2 = 1.0f;
        objN = 1.0f;
        cosTheta = abs(cosTheta);

        normal = -normal;
        // cosTheta = -dot(w_i, normal);
    }
    // compute cosPhi, refracted/transmitted angle. formula on pg.14 of slides.
    float r = n1/n2;
    float sinThetaSqr = 1 - (cosTheta*cosTheta);
    float criticalTerm = r*r*sinThetaSqr;
    if(criticalTerm > 1)
    {
        // then inside of square root is negative, total internal reflection case.
        // perform reflection only, Fr = 1, Ft = 0.
        Vec3f w_r =  Reflect(normal, w_o, mat.roughness);

        // Do reflection?
        Ray ray;
        ray.origin = x + normal * scene.shadow_ray_epsilon;
        ray.dir = w_r;
        ray.hitInfo.hasHit = false;
        ray.hitInfo.minT = INFINITY;
        ray.refractiveIndexOfCurrentMedium = raysN; // it didnt change medium, reflected back into the one its coming from.

        // todo: attenuate if travelling inside an object.
        IntersectObjects(ray);
        if(ray.hitInfo.hasHit)
        {
            if(ray.refractiveIndexOfCurrentMedium > 1.00001f)
            {
                // non-vacuum, attenuate
                return BeersLaw(ray.hitInfo.minT, mat.absorptionCoefficient, PerformShading(ray, ray.origin, recDepth-1));
            }
            else return PerformShading(ray, ray.origin, recDepth-1);
        }
        else return Vec3f{0,0,0};
    }
    else
    {
        // both reflection and refraction, perform both, find out the ratio too.
        float cosPhi = sqrt(1-criticalTerm);
        float n2cosTheta= n2*cosTheta;
        float n1cosPhi = n1*cosPhi;

        float rparallel = (n2cosTheta - n1cosPhi) / (n2cosTheta + n1cosPhi);
        float rperp = (n1*cosTheta - n2*cosPhi) / (n1*cosTheta + n2*cosPhi);
        float rReflect = (rparallel * rparallel + rperp * rperp) / 2;
        float rRefract = 1 - rReflect;

       Vec3f reflectedRaysColor;
        {
            Vec3f w_reflected = Reflect(normal, w_o, mat.roughness);

            // Again do reflection here
            Ray ray;
            ray.origin = x + normal * scene.shadow_ray_epsilon;
            ray.dir = w_reflected;
            ray.hitInfo.hasHit = false;
            ray.hitInfo.minT = INFINITY;
            ray.refractiveIndexOfCurrentMedium = raysN; // it didnt change medium, reflected back into the one its coming from.

            IntersectObjects(ray);
            reflectedRaysColor = Vec3f{0,0,0};

            if(ray.hitInfo.hasHit)
            {
                reflectedRaysColor = PerformShading(ray, ray.origin, recDepth-1);
                // no attenuation on reflect case, since we are not INSIDE an object?
                // what about total internal reflection though?
                if(ray.refractiveIndexOfCurrentMedium > 1.00001f)
                {
                    // non-vacuum, attenuate
                    reflectedRaysColor = BeersLaw(ray.hitInfo.minT, mat.absorptionCoefficient, reflectedRaysColor);
                }
            }
        }
       
        // BUT also do refraction, 1 ray split into 2 rays in total.
        // wt = (d + ncosTheta)(n1/n2) - ncosPhi, page 16 of slides.
        Vec3f refractedRaysColor;
        {
            Vec3f w_refracted = (d + normal * cosTheta) * r - normal * cosPhi;
            // also perturb the refraction ray?

            Vec3f u,v;
            GetOrthonormalBasis(w_refracted, u ,v);

            float psi1 = GetRandom01();
            float psi2 = GetRandom01();

            w_refracted = makeUnit(w_refracted + (u*psi1 + v*psi2)* mat.roughness);
            

            Ray ray;
            ray.dir = makeUnit(w_refracted);
            ray.origin = x + (-normal) * scene.shadow_ray_epsilon;
            ray.hitInfo.hasHit = false;
            ray.hitInfo.minT = INFINITY;
            ray.refractiveIndexOfCurrentMedium = objN; // it changed medium, transmitted to other one.

            if(isEntering){

            ray.refractiveIndexOfCurrentMedium = objN; // it changed medium, transmitted to other one.
            }
            else{

            ray.refractiveIndexOfCurrentMedium = 1.0f; // it changed medium, transmitted to other one.
            }
            // Maybe only intersect with current object if its entering.
            // cuz it has to intersect without leaving if  this is a non-infinite object, or a plane.
            
            IntersectObjects(ray);
            refractedRaysColor = Vec3f{0,0,0};

            if(ray.hitInfo.hasHit)
            {
                refractedRaysColor = PerformShading(ray, ray.origin, recDepth-1);
                if(ray.refractiveIndexOfCurrentMedium > 1.001f)
                {
                    // non-vacuum, attenuate
                    refractedRaysColor = BeersLaw(ray.hitInfo.minT, mat.absorptionCoefficient, refractedRaysColor);
                }
                else{

                }
            }
        }
        return reflectedRaysColor * rReflect + refractedRaysColor * rRefract;
    }


}
Vec3f Raytracer::BeersLaw(float x, Vec3f c, Vec3f L_0){
    Vec3f res;
    res.x = L_0.x * std::exp(-c.x * x);
    res.y = L_0.y * std::exp(-c.y * x);
    res.z = L_0.z * std::exp(-c.z * x);
    return res;
}
Vec3f Raytracer::Reflect(Vec3f& normal, Vec3f& w_o, float roughness)
{
    // assign a generator per thread, this is costly to do here, and stupid.

    Vec3f r = makeUnit((normal * 2.0f * dot(normal, w_o)) - w_o);
    if(roughness > 0.001)
    {
        Vec3f u,v;
        GetOrthonormalBasis(r, u ,v);

        float psi1 = GetRandom01();
        float psi2 = GetRandom01();

        Vec3f rPrime = makeUnit(r + (u*psi1 + v*psi2)*roughness);
        return rPrime;
    }
    else return r;
}

Vec3f Raytracer::ComputeMirrorReflection(Material& mat, Vec3f& w_o, Vec3f& normal, Vec3f& intersectionPoint, int recursionDepth)
{
    if (recursionDepth == 0)
    {
        return Vec3f{0,0,0};
    }
    else 
    {
        Vec3f w_r =  Reflect(normal, w_o, mat.roughness);

        Ray ray;
        ray.origin = intersectionPoint + normal * scene.shadow_ray_epsilon;
        ray.dir = w_r;
        ray.hitInfo.hasHit = false;
        ray.hitInfo.minT = INFINITY;
        ray.refractiveIndexOfCurrentMedium = 1.0f;

        IntersectObjects(ray);
        if(ray.hitInfo.hasHit)
        {
            return mat.mirror * PerformShading(ray, ray.origin, recursionDepth-1);
        }
        else return Vec3f{0,0,0};
    }
}

Vec3f Raytracer::GetAmbient(Vec3f& reflectance, Vec3f& ambientLightColor){
    return ambientLightColor * reflectance;
}

Vec3f Raytracer::GetDiffuse(Vec3f& reflectance, Vec3f& w_i, Vec3f& normal, Vec3f& receivedIrradiance){
    float costheta = std::max(0.0f, dot(w_i, normal));
    return reflectance * receivedIrradiance * costheta;
}

Vec3f Raytracer::GetSpecular(Vec3f& reflectance, float phongExp, Vec3f& w_in, Vec3f& w_out, Vec3f& normal, Vec3f& receivedIrradiance)
{
    Vec3f half = (w_in + w_out) / len(w_in + w_out);
    float cosAlpha = std::max(0.0f, dot(normal, half));
    return reflectance * receivedIrradiance * std::pow(cosAlpha, phongExp);
}

bool Raytracer::IsInShadow(Vec3f& intersectionPoint, Vec3f& lightPos)
{
    // intersect with all objects in the scene, 
    // Return true if intersects with any objects before the light source.
    float lightSourceT = len(lightPos - intersectionPoint);
    Ray shadowRay;
    shadowRay.dir = makeUnit(lightPos - intersectionPoint);

    shadowRay.origin = intersectionPoint + shadowRay.dir * scene.shadow_ray_epsilon;
    shadowRay.hitInfo.hasHit = false;
    shadowRay.hitInfo.minT = lightSourceT + 0.1f;

    for(int i = 0; i < scene.meshes.size(); i++){
        Shape* mesh = scene.meshes[i];
        mesh->Intersect(shadowRay);

        if(shadowRay.hitInfo.hasHit && shadowRay.hitInfo.minT < lightSourceT)
        {
            // valid hit, in shadow.
            return true;
        }
    }
    // Intersect with all Spheres
    for(int i = 0; i < scene.spheres.size(); i++){
        Sphere* s = scene.spheres[i];
        s->Intersect(shadowRay);
        if(shadowRay.hitInfo.hasHit && shadowRay.hitInfo.minT < lightSourceT)
        {
            // valid hit, in shadow.
            return true;
        }
    }
    for(int i = 0; i < scene.triangles.size(); i++){
        Shape* s = scene.triangles[i];
        s->Intersect(shadowRay);
        if(shadowRay.hitInfo.hasHit && shadowRay.hitInfo.minT < lightSourceT)
        {
            // valid hit, in shadow.
            return true;
        }
    }

    return false;
}
void Raytracer::IntersectObjects(Ray& ray)
{
    for(int i = 0; i < scene.meshes.size(); i++){
        Shape* mesh = scene.meshes[i];
        mesh->Intersect(ray);
    }

    // Intersect with all Spheres
    for(int i = 0; i < scene.spheres.size(); i++){
        Shape* sphere = scene.spheres[i];
        sphere->Intersect(ray);
    }

     // Intersect with all triangles
    for(int i = 0; i < scene.triangles.size(); i++){
        Shape* s = scene.triangles[i];
        s->Intersect(ray);
    }
}

Ray  Raytracer::GenerateRay(int i, int j, Camera& cam)
{
    Ray ray;
    Vec3f imagePlanePos = cam.GetImagePlanePosition(i,j);

    ray.origin = cam.position;
    ray.dir = imagePlanePos - ray.origin;
    
    ray.hitInfo.hasHit = false;
    ray.hitInfo.minT = INFINITY;
    ray.refractiveIndexOfCurrentMedium = 1.0f;
    return ray;
}