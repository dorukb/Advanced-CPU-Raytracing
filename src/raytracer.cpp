#include "raytracer.hpp"
#include <iostream>
#include <random>

using namespace DorkTracer;

Raytracer::Raytracer(Scene& scene){
    this->scene = scene;

    this->roughnessRandomGenerator = std::mt19937(rand());
    this->roughnessRandomDistro =  std::uniform_real_distribution<>(-0.5f, 0.5f);
    
    this->dofLensSampleGenerator = std::mt19937(rand());
    this->dofLensSampleDistro = std::uniform_real_distribution<>(-1.0f, 1.0f);
}

float Raytracer::GetRandom()
{
    return this->roughnessRandomDistro(this->roughnessRandomGenerator);
}
float Raytracer::GetLensSample(){
    return this->dofLensSampleDistro(this->dofLensSampleGenerator);
}

Vec3f Raytracer::RenderPixel(int i, int j, Camera& cam)
{
    return PerPixel(i,j, cam);
}

Vec3f Raytracer::PerPixel(int coordX, int coordY, Camera& cam)
{
    Ray ray = GenerateRay(coordX, coordY, cam);
    
    // Intersect with all Objects in the scene to find the closest intersection point, if any.
    IntersectObjects(ray);
    
    if(ray.hitInfo.hasHit)
    {
       return PerformShading(ray, cam.position, scene.max_recursion_depth);
    }
    else if(scene.bgTexture != nullptr)
    {
        // sample the bg texture.
        float u= coordX / (float)cam.imageWidth;
        float v = coordY / (float)cam.imageHeight;
        Vec3f sampledColor = scene.bgTexture->GetRGBSample(u, v);
        return sampledColor;
    }
    else if( scene.sphericalEnvLights.size() > 0)
    {
        // Has environment light, sample that.
        return scene.sphericalEnvLights[0]->GetSample(ray.dir);
    }
    else return Vec3f(scene.background_color.x, scene.background_color.y, scene.background_color.z);
}

Vec3f Raytracer::PerformShading(Ray& ray, Vec3f& eyePos, int recursionDepth)
{
    ray.hitInfo.hitPoint = ray.origin + ray.dir * ray.hitInfo.minT;
    // Actually perform shading, depending on Material type.

    Vec3f color{0,0,0};
    Material& mat = scene.materials[ray.hitInfo.matId-1];
    Shape* shape = ray.hitInfo.hitShape;

    Vec3f w_o = makeUnit(eyePos - ray.hitInfo.hitPoint);
    float refractiveIndexOfVacuum = 1.00001;
    bool travellingInsideAnObject = ray.refractiveIndexOfCurrentMedium > refractiveIndexOfVacuum;

    // If Replace_all mode, disable shading, directly output texture color.
    if(shape->HasReplaceAllTexture()){
        return shape->replaceAll->GetRGBSample(ray.hitInfo.hitUV.x, ray.hitInfo.hitUV.y);
    }


    if(!travellingInsideAnObject)
    {
        color = color + GetAmbient(mat.ambient, scene.ambient_light);
        // Repeat for each light source.
        for(int l = 0; l < scene.point_lights.size(); l++)
        {
            // Skip all calculation if in shadow, no color "contribution".
            PointLight& light = scene.point_lights[l];
            if(IsInShadow(ray, light.position)){
                continue;
            }

            Vec3f w_i = makeUnit(light.position - ray.hitInfo.hitPoint);
            
            float distToLight = len(light.position - ray.hitInfo.hitPoint);
            Vec3f receivedIrradiance = light.intensity / (distToLight * distToLight);

            color = color + GetDiffuse(shape, mat.diffuse, w_i, ray, receivedIrradiance);
            color = color + GetSpecular(shape, mat.specular, ray, mat.phong_exponent, w_i, w_o, receivedIrradiance);
        }    

        for(int i = 0; i < scene.areaLights.size(); i++)
        {
            // Skip all calculation if in shadow, no color "contribution".
            AreaLight* areaLight = scene.areaLights[i];
            Vec3f samplePos = areaLight->GetSample();
            // Vec3f samplePos = areaLight->GetSamplePosition(ray.lightSampleX, ray.lightSampleY);

            if(IsInShadow(ray, samplePos)){
                continue;
            }

            // TODO: fix w_i calculation convention across lights.

            Vec3f w_i = samplePos - ray.hitInfo.hitPoint;
            float distToLight = len(w_i);
            float dSqr = distToLight * distToLight;
            w_i = w_i / distToLight; // normalize.

            float lCostheta = dot(areaLight->normal, -w_i);
            if(lCostheta < 0){
                lCostheta = dot(areaLight->normal, w_i);
            }
            Vec3f receivedIrradiance = areaLight->radiance * (areaLight->area * lCostheta / dSqr);
          
            color = color + GetDiffuse(shape, mat.diffuse, w_i, ray, receivedIrradiance);
            color = color + GetSpecular(shape, mat.specular, ray, mat.phong_exponent, w_i, w_o, receivedIrradiance);
        }        
        for(int i = 0; i < scene.sphericalEnvLights.size(); i++)
        {
            SphericalEnvironmentLight* envLight = scene.sphericalEnvLights[i];

            // if(IsInShadow(ray, samplePos)){
            //     continue;
            // }
            Vec3f sampleDir = envLight->GetDirection(ray.hitInfo.normal);
            Vec3f receivedIrradiance = envLight->GetSample(sampleDir);
          
            Vec3f w_i = ray.hitInfo.normal;
            color = color + GetDiffuse(shape, mat.diffuse, w_i, ray, receivedIrradiance);
            color = color + GetSpecular(shape, mat.specular, ray, mat.phong_exponent, w_i, w_o, receivedIrradiance);
        }    

        for(int i = 0; i < scene.directionalLights.size(); i++)
        {
            DirectionalLight* light = scene.directionalLights[i];
            if(IsInShadowDirectional(ray, light->dir)){
                continue;
            }

            Vec3f w_i = -(light->dir);
            color = color + GetDiffuse(shape, mat.diffuse, w_i, ray, light->radiance);
            color = color + GetSpecular(shape, mat.specular, ray, mat.phong_exponent, w_i, w_o, light->radiance);
        }
        
        for(int i = 0; i < scene.spotLights.size(); i++)
        {
            SpotLight* light = scene.spotLights[i];
            if(IsInShadow(ray, light->pos)){
                continue;
            }

            Vec3f w_i = makeUnit(light->pos - ray.hitInfo.hitPoint);
            Vec3f receivedIrradiance = light->GetIrradiance(ray.hitInfo.hitPoint);
            color = color + GetDiffuse(shape, mat.diffuse, w_i, ray, receivedIrradiance);
            color = color + GetSpecular(shape, mat.specular, ray, mat.phong_exponent, w_i, w_o, receivedIrradiance);
        }    

    }
 
    if(mat.type == Material::Mirror)
    {
        // Compute radiance along the ideal reflection ray
        color = color + ComputeMirrorReflection(ray, mat, w_o, recursionDepth);
    }
    else if(mat.type == Material::Dielectric)
    {
        // Both reflection & transmission
        float n1 = ray.refractiveIndexOfCurrentMedium;
        color = color + ComputeDielectricFresnelReflectionAndRefraction(ray, mat, w_o, n1, mat.refractiveIndex, recursionDepth);
    }
    else if(mat.type == Material::Conductor)
    {
        // Reflection & absorption, no tranmission.
        color = color + ComputeConductorFresnelReflection(ray, mat, w_o, recursionDepth);
    }
    // Vec3i c = clamp(color);
    // color.x = c.x;
    // color.y = c.y;
    // color.z = c.z;
    
    return color;
}

Vec3f Raytracer::ComputeConductorFresnelReflection(Ray& originalRay, Material& mat, Vec3f& w_o, int recDepth)
{
    if(recDepth == 0){
        return Vec3f{0,0,0}; 
    }

    Vec3f d = -w_o;
    float cosTheta = -dot(d, originalRay.hitInfo.normal);
    float n2 = mat.refractiveIndex;
    float k2 = mat.conductorAbsorptionIndex;

    // compute common intermediate terms
    float n2k2 = n2*n2 + k2*k2;
    float n2cosTheta2 = 2 * n2 * cosTheta;
    float cosThetaSqr = cosTheta * cosTheta;

    float rs = (n2k2 - n2cosTheta2 + cosThetaSqr) / (n2k2 + n2cosTheta2 + cosThetaSqr);
    float rp = (n2k2 * cosThetaSqr - n2cosTheta2 + 1) / (n2k2 * cosThetaSqr + n2cosTheta2 + 1);

    float reflectRatio = 0.5 * (rs + rp);
    // transmit ratio = 0.

    // Cast reflected ray
    if(reflectRatio > 0.0001)
    {
        Vec3f reflectedRaysColor;
        {
            Vec3f w_reflected = Reflect(originalRay.hitInfo.normal, w_o, mat.roughness);

            Vec3f origin = originalRay.hitInfo.hitPoint + originalRay.hitInfo.normal * scene.shadow_ray_epsilon;
            Ray reflectedRay = GenerateSecondaryRay(originalRay, w_reflected, origin);
            reflectedRay.refractiveIndexOfCurrentMedium = 1.0f;

            IntersectObjects(reflectedRay);
            if(reflectedRay.hitInfo.hasHit)
            {
                // no attenuation on reflect case, since we are not INSIDE an object?
                reflectedRaysColor = mat.mirror * PerformShading(reflectedRay, reflectedRay.origin, recDepth-1);
            }
            else reflectedRaysColor = Vec3f{0,0,0};
        }

        return reflectedRaysColor * reflectRatio;
    }
    else return Vec3f{0,0,0};
 
}
// Ratio of reflected over refracted light amount, for dielectric - dielectric interfacing.
// cosTheta = dot(w_i, normal)
// n1: refractive index of the medium incoming ray is currently in
// n2: refractive index of the medium ray is entering, or equivalently that it will be transmitted in.
// float GetFresnelReflectionRatioDielectric(Vec3f& w_i, Vec3f& normal, float n1, float n2)

Vec3f Raytracer::ComputeDielectricFresnelReflectionAndRefraction(Ray& originalRay, Material& mat, Vec3f& w_o, float n1, float n2, int recDepth)
{
    if(recDepth == 0){
        return Vec3f{0,0,0}; 
    }

    Vec3f d = -w_o;
    Vec3f modifiedNormal = originalRay.hitInfo.normal;
    float cosTheta = -dot(d, modifiedNormal);
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

        modifiedNormal = -modifiedNormal;
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
        Vec3f w_r =  Reflect(modifiedNormal, w_o, mat.roughness);

        // Note: normal might be flipped temporarily, but its not written back to original ray.
        // flipped version must be used throughout this calculation.
        Vec3f newOrigin = originalRay.hitInfo.hitPoint + modifiedNormal * Scene::shadow_ray_epsilon;
        Ray reflectedRay = GenerateSecondaryRay(originalRay, w_r, newOrigin);

        IntersectObjects(reflectedRay);
        Vec3f reflectedRaysColor{0,0,0}; 
        if(reflectedRay.hitInfo.hasHit)
        {
            reflectedRaysColor = PerformShading(reflectedRay, reflectedRay.origin, recDepth-1);
            if(reflectedRay.refractiveIndexOfCurrentMedium > 1.0001){
                // attenuate if travelling inside an object.
                reflectedRaysColor = BeersLaw(reflectedRay.hitInfo.minT, mat.absorptionCoefficient, reflectedRaysColor);
            }
        }
        return reflectedRaysColor;
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


        // Handle Reflected part
        Vec3f w_reflected = Reflect(modifiedNormal, w_o, mat.roughness);
        Vec3f newOrigin = originalRay.hitInfo.hitPoint + modifiedNormal * Scene::shadow_ray_epsilon;
        Ray reflectedRay = GenerateSecondaryRay(originalRay, w_reflected, newOrigin);
        IntersectObjects(reflectedRay);
        if(isEntering)
        {
            reflectedRay.refractiveIndexOfCurrentMedium = objN;
        }
        else
        {
            reflectedRay.refractiveIndexOfCurrentMedium = 1.0f;
        }

        Vec3f reflectedRaysColor = Vec3f{0,0,0};
        if(reflectedRay.hitInfo.hasHit)
        {
            // TODO: mayue attenuate always?
            reflectedRaysColor = PerformShading(reflectedRay, reflectedRay.origin, recDepth-1);
            if(reflectedRay.refractiveIndexOfCurrentMedium > 1.00001f)
            {
                // non-vacuum, attenuate
                reflectedRaysColor = BeersLaw(reflectedRay.hitInfo.minT, mat.absorptionCoefficient, reflectedRaysColor);
            }
        }    
        else if(scene.sphericalEnvLights.size() > 0)
        {
            // Has environment light, sample that.
            SphericalEnvironmentLight* envLight = scene.sphericalEnvLights[0];
            reflectedRaysColor = envLight->GetSample(reflectedRay.dir);
        }
       
        // BUT also do refraction, 1 ray split into 2 rays in total.
        // Handle Refraction part.
        Vec3f refractedRaysColor;
        {
            // wt = (d + ncosTheta)(n1/n2) - ncosPhi, page 16 of slides.
            Vec3f w_refracted = (d + modifiedNormal * cosTheta) * r - modifiedNormal * cosPhi;
           
            // Perturb the refraction if has a rough material.
            if(mat.roughness > 0.001)
            {
                Vec3f u,v;
                GetOrthonormalBasis(w_refracted, u ,v);

                float psi1 = GetRandom();
                float psi2 = GetRandom();

                w_refracted = makeUnit(w_refracted + (u*psi1 + v*psi2)* mat.roughness);
            }
            else w_refracted = makeUnit(w_refracted);

            Vec3f newOrigin = originalRay.hitInfo.hitPoint + (-modifiedNormal) * Scene::shadow_ray_epsilon;
            Ray refractedRay = GenerateSecondaryRay(originalRay, w_refracted, newOrigin);

            if(isEntering)
            {
                refractedRay.refractiveIndexOfCurrentMedium = objN;
            }
            else
            {
                refractedRay.refractiveIndexOfCurrentMedium = 1.0f;
            }
            // Maybe only intersect with current object if its entering.
            // cuz it has to intersect without leaving if  this is a non-infinite object, or a plane.
            
            IntersectObjects(refractedRay);
            refractedRaysColor = Vec3f{0,0,0};

            if(refractedRay.hitInfo.hasHit)
            {
                refractedRaysColor = PerformShading(refractedRay, refractedRay.origin, recDepth-1);
                if(refractedRay.refractiveIndexOfCurrentMedium > 1.001f)
                {
                    // non-vacuum, attenuate
                    refractedRaysColor = BeersLaw(refractedRay.hitInfo.minT, mat.absorptionCoefficient, refractedRaysColor);
                }
            }
            else if(scene.sphericalEnvLights.size() > 0)
            {
                // Has environment light, sample that instead of empty black.
                SphericalEnvironmentLight* envLight = scene.sphericalEnvLights[0];
                refractedRaysColor = envLight->GetSample(reflectedRay.dir);
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
    Vec3f r = makeUnit((normal * 2.0f * dot(normal, w_o)) - w_o);
    if(roughness > 0.001)
    {
        Vec3f u,v;
        GetOrthonormalBasis(r, u ,v);

        float psi1 = GetRandom();
        float psi2 = GetRandom();

        Vec3f rPrime = makeUnit(r + (u*psi1 + v*psi2)*roughness);
        return rPrime;
    }
    else return r;
}

Vec3f Raytracer::ComputeMirrorReflection(Ray& originalRay, Material& mat, Vec3f& w_o, int recursionDepth)
{
    if (recursionDepth == 0)
    {
        return Vec3f{0,0,0};
    }
    else 
    {
        Vec3f w_r =  Reflect(originalRay.hitInfo.normal, w_o, mat.roughness);

        Vec3f origin = originalRay.hitInfo.hitPoint + originalRay.hitInfo.normal * Scene::shadow_ray_epsilon;
        Ray reflectedRay = GenerateSecondaryRay(originalRay, w_r, origin);
        reflectedRay.refractiveIndexOfCurrentMedium = 1.0f;

        IntersectObjects(reflectedRay);
        if(reflectedRay.hitInfo.hasHit)
        {
            return mat.mirror * PerformShading(reflectedRay, reflectedRay.origin, recursionDepth-1);
        }
        else if(scene.sphericalEnvLights.size() > 0)
        {
            // Has environment light, sample that.
            SphericalEnvironmentLight* envLight = scene.sphericalEnvLights[0];

            // Vec3f sampleDir = envLight->GetDirection(originalRay.hitInfo.normal);
            Vec3f receivedIrradiance = envLight->GetSample(reflectedRay.dir);
            return mat.mirror * receivedIrradiance;
        }
        else return Vec3f{0,0,0}; 
    }
}

Vec3f Raytracer::GetAmbient(Vec3f& reflectance, Vec3f& ambientLightColor){
    return ambientLightColor * reflectance;
}

Vec3f Raytracer::GetDiffuse(Shape* shape, Vec3f& k_d, Vec3f& w_i, Ray& ray, Vec3f& receivedIrradiance)
{
    Vec3f reflectance = k_d;
    if(shape->HasDiffuseTexture())
    {
        // TODO: calculate UV for hit point using vertices & texCoords.
        // normalize, obtain k_d in [0,1] range.
        Vec3f textureKd;

        if(shape->diffuseTex->IsGenerated()){
            // such as Perlin noise. Sampling process is different.
            auto hp = ray.hitInfo.hitPoint;
            float perlinSample = shape->diffuseTex->GetSampleFromWorldPos(hp.x, hp.y, hp.z);
            textureKd.x = textureKd.y = textureKd.z = perlinSample;
        }
        else
        {
            Vec2f& uv = ray.hitInfo.hitUV;
            textureKd = shape->diffuseTex->GetRGBSample(uv.x, uv.y) / 255.0f;
        }

        if(shape->diffuseTex->operationMode == Texture::OperationMode::Blend)
        {
            //Equally mix the existing material kd value with the tex sample.
            reflectance = (textureKd + k_d) / 2.0f;
        }
        else // Replace/overwrite k_d with tex value
        {
            reflectance = textureKd;
        }
    }
    float costheta = std::max(0.0f, dot(w_i, ray.hitInfo.normal));
    return reflectance * receivedIrradiance * costheta;
}

Vec3f Raytracer::GetSpecular(Shape* shape, Vec3f& k_s, Ray& ray, float phongExp, Vec3f& w_in, Vec3f& w_out, Vec3f& receivedIrradiance)
{
    Vec3f reflectance = k_s;
    if(shape->specularTex != nullptr)
    {
        // TODO: calculate UV for hit point using vertices & texCoords.
        // normalize, obtain k_d in [0,1] range.
        Vec2f& uv = ray.hitInfo.hitUV;
        Vec3f textureKs = shape->specularTex->GetRGBSample(uv.x, uv.y) / 255.0f;
        reflectance = textureKs;
    }

    Vec3f half = (w_in + w_out) / len(w_in + w_out);
    float cosAlpha = std::max(0.0f, dot(ray.hitInfo.normal, half));
    return reflectance * receivedIrradiance * std::pow(cosAlpha, phongExp);
}
bool Raytracer::IsInShadowDirectional(Ray& originalRay, Vec3f& lightDir)
{
    Ray shadowRay;
    shadowRay.dir = -lightDir;
    shadowRay.origin = originalRay.hitInfo.hitPoint + originalRay.hitInfo.normal * Scene::shadow_ray_epsilon;
    shadowRay.hitInfo.hasHit = false;
    
    //directional lights are INFINITELY far away.
    shadowRay.hitInfo.minT = INFINITY;
    shadowRay.motionBlurTime = originalRay.motionBlurTime;
    return CastShadowRay(shadowRay, INFINITY);
}
bool Raytracer::IsInShadow(Ray& originalRay, Vec3f& lightPos)
{
    // intersect with all objects in the scene, 
    // Return true if intersects with any objects before the light source.
    Ray shadowRay;

    shadowRay.dir = lightPos - originalRay.hitInfo.hitPoint;
    float lightSourceT = len(shadowRay.dir);

    shadowRay.dir = shadowRay.dir / lightSourceT; // cheap normalization.
    shadowRay.origin = originalRay.hitInfo.hitPoint + originalRay.hitInfo.normal * Scene::shadow_ray_epsilon;

    shadowRay.hitInfo.hasHit = false;
    shadowRay.hitInfo.minT = lightSourceT + 0.01f;
    shadowRay.motionBlurTime = originalRay.motionBlurTime;

    return CastShadowRay(shadowRay, lightSourceT);
}
bool Raytracer::CastShadowRay(Ray& shadowRay, float lightSourceT)
{
    for(int i = 0; i < scene.meshes.size(); i++)
    {
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

    // If missed all objects, Intersect with the Spherical Environment Map (or Skybox)
    if(!ray.hitInfo.hasHit)
    {

    }

}


Ray Raytracer::GenerateSecondaryRay(Ray& original, Vec3f& newDir, Vec3f& newOrigin)
{
    Ray ray;
    ray.dir = newDir;
    ray.origin = newOrigin;

    ray.hitInfo.hasHit = false;
    ray.hitInfo.minT = INFINITY;

    // it didnt change medium, reflected back into the one its coming from.
    ray.refractiveIndexOfCurrentMedium = original.refractiveIndexOfCurrentMedium;

    ray.motionBlurTime = original.motionBlurTime;
    // ray.lightSampleX = original.lightSampleX;
    // ray.lightSampleY = original.lightSampleY;
    
    // if(scene.areaLights.size() > 0)
    // {    
    //     ray.lightSampleX = GetRandom();
    //     ray.lightSampleY = GetRandom();
    // }
    // else{
    //     ray.lightSampleX = ray.lightSampleY = 0.0f;
    // }
    
    return ray;
}
Ray Raytracer::GenerateRay(int i, int j, Camera& cam)
{
    Ray ray;
    Vec3f imagePlanePos = cam.GetImagePlanePosition(i,j);

    ray.origin = cam.position;

    bool depthOfFieldEnabled = cam.apertureSize > 0.0001;
    if(depthOfFieldEnabled)
    {
        Vec3f apertureSamplePos = ray.origin;

        float first01 = GetLensSample();
        apertureSamplePos = apertureSamplePos + cam.up * (first01 * cam.apertureSize * 0.5f);

        float second01 = GetLensSample();
        apertureSamplePos = apertureSamplePos + cam.right * (second01 * cam.apertureSize * 0.5f);

        Vec3f dir = makeUnit(ray.origin - imagePlanePos);
        float tFd = cam.focusDistance / dot(dir, cam.gaze);
        Vec3f bentScenePoint = ray.origin + dir * tFd;

        ray.dir = makeUnit(bentScenePoint - apertureSamplePos);
        ray.origin = apertureSamplePos;

    }
    else
    {
        ray.dir = makeUnit(imagePlanePos - ray.origin);
    }

    // if(scene.areaLights.size() > 0)
    // {    
    //     ray.lightSampleX = GetRandom();
    //     ray.lightSampleY = GetRandom();
    // }
    // else{
    //     ray.lightSampleX = ray.lightSampleY = 0.0f;
    // }
    
    ray.hitInfo.hasHit = false;
    ray.hitInfo.minT = INFINITY;
    ray.refractiveIndexOfCurrentMedium = 1.0f;

    // ray.motionBlurTime = -999.0f;
    ray.motionBlurTime = GetRandom() + 0.5f;

    
    return ray;
}