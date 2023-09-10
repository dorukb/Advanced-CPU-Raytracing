#include "raytracer.hpp"
#include <iostream>
#include <random>

using namespace DorkTracer;

Raytracer::Raytracer(Scene& scene){
    this->scene = scene;

    this->randGen = std::mt19937(rand());
    this->roughnessRandomDistro =  std::uniform_real_distribution<>(-0.5f, 0.5f);
    this->normalizedDistro =  std::uniform_real_distribution<>(0.0f, 1.0f);

    this->dofLensSampleGenerator = std::mt19937(rand());
    this->dofLensSampleDistro = std::uniform_real_distribution<>(-1.0f, 1.0f);
}

void Raytracer::EnablePathTracing(RendererParams params)
{
    this->pathTracingEnabled = true;
    this->rendererParams = params;
}

float Raytracer::GetNormalizedRandom()
{
    // Returns a random float in range [0,1].
    return this->normalizedDistro(this->randGen);
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
    // if(recursionDepth == 0) return Vec3f(0,0,0);

    ray.hitInfo.hitPoint = ray.origin + ray.dir * ray.hitInfo.minT;
    // Actually perform shading, depending on Material type.

    Vec3f color{0,0,0};
    Material& mat = scene.materials[ray.hitInfo.matId-1];
    Shape* shape = ray.hitInfo.hitShape;

    Vec3f w_o = makeUnit(eyePos - ray.hitInfo.hitPoint);
    float refractiveIndexOfVacuum = 1.00001;
    bool travellingInsideAnObject = ray.refractiveIndexOfCurrentMedium > refractiveIndexOfVacuum;


    if(mat.type == Material::Emissive)
    {
        return mat.radiance * 2.0f * M_PI;
    }

    // If Replace_all mode, disable shading, directly output texture color.
    if(shape->HasReplaceAllTexture()){
        return shape->replaceAll->GetRGBSample(ray.hitInfo.hitUV.x, ray.hitInfo.hitUV.y);
    }

    int hitLightMeshId = -1; // to prevent sampling the same object light twice if NextEventEstimation is also active.
    if(rendererParams.pathTracingEnabled)
    {
        // send global illumination ray here.
        color = color + ComputeGlobalIllumination(ray, mat, w_o, recursionDepth, hitLightMeshId);
    }

    // TODO: sampleDirectLight should be TRUE iff either
    // 1) NO path tracing
    // 2) path tracing with NextEventEstimation enabled.
    bool sampleDirectLight = !rendererParams.pathTracingEnabled 
                            || (rendererParams.pathTracingEnabled && rendererParams.nextEventEstimationEnabled);

    if(!travellingInsideAnObject && sampleDirectLight)
    {
        color = color + GetAmbient(mat.ambient, scene.ambient_light);
        color = color + SampleDirectLighting(ray, mat, w_o, hitLightMeshId);
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


    if(isnan(color.x) || isnan(color.y) || isnan(color.z))
    {
        std::cout<<"nan color!!" << std::endl;
    }

    return color;
}
Vec3f Raytracer::ComputeGlobalIllumination(Ray& ray, Material& orgMat, Vec3f& w_o, int recDepth, int& hitMeshLightId)
{ 
    // if(recDepth == 0){
    //     return Vec3f{0,0,0}; 
    // }
    if(rendererParams.russianRouletteEnabled){
        // TODO: Implement russian roulette.
        // terminate the ray based on probability & throughput
        float probTest = GetNormalizedRandom();
        // float terminationProb = 1 - ray.throughput;
        float maxThroughput = std::max(ray.throughput.x, std::max(ray.throughput.x, ray.throughput.z));
        if(probTest > maxThroughput && recDepth <= 0){
            // the lower the ray throughput, the higher the chance to eliminate it.
            return Vec3f{0,0,0};
        }
        ray.throughput = ray.throughput / maxThroughput;
    }
    else if(recDepth <= 0){
        return Vec3f{0,0,0};
    }
   

    float rand1 = GetNormalizedRandom(); // [0,1]
    float rand2 = GetNormalizedRandom(); // [0,1]

    float phi = 2 * M_PI * rand1;
    float theta = 0.0f;

    if(rendererParams.sampleImportance){
        theta = std::asin(std::sqrt(rand2)); // Importance
    }
    else{
        theta = std::acos(rand2); // uniform
    }

    // construct ONB -> uvw, whetere v is aligned with surface normal at hit point.
    Vec3f u,v;
    GetOrthonormalBasis(ray.hitInfo.normal, u, v);
    Vec3f newDir = u * std::sin(theta) * std::cos(phi) 
                 + ray.hitInfo.normal * std::cos(theta)
                 + v * std::sin(theta) * std::sin(phi);

    newDir = makeUnit(newDir);

    // Create new ray
    Vec3f newOrigin = ray.hitInfo.hitPoint + ray.hitInfo.normal * 0.0001;
    Ray globalRay = GenerateSecondaryRay(ray, newDir, newOrigin);

    // cast the new ray
    IntersectObjects(globalRay);
    Vec3f globalRaysColor{0,0,0}; 
    if(globalRay.hitInfo.hasHit)
    {
        Material& mat = scene.materials[globalRay.hitInfo.matId-1];
        if(mat.type == Material::Emissive){
            hitMeshLightId = globalRay.hitInfo.hitShape->id;
        }
        // Compute the contribution.
        Vec3f transferredLight = PerformShading(globalRay, globalRay.origin, recDepth-1);
        globalRaysColor = Shade(ray, orgMat, globalRay.dir, w_o, transferredLight) * 2.0f * M_PI;
    }
    return globalRaysColor;
}
Vec3f Raytracer::Shade(Ray& ray, Material& mat, Vec3f& w_i, Vec3f& w_o, Vec3f& Li)
{
    Vec3f color;
    Shape* s = ray.hitInfo.hitShape;
    if(mat.HasBRDF())
    {
        float costheta_i = std::max(0.0f, dot(w_i, ray.hitInfo.normal));
        Vec3f kd = GetDiffuseReflectanceCoeff(ray, s, mat);
        Vec3f ks = GetSpecularReflectanceCoeff(ray, s, mat);
        Vec3f res =  mat.brdf->apply(mat, kd, ks, w_i, w_o, ray.hitInfo.normal);
        ray.throughput = ray.throughput * res;
        return res * Li * costheta_i;
    }
    else return GetDiffuse(ray, s, mat, w_i, Li) + GetSpecular(ray, s, mat, w_i, w_o, Li);
}

Vec3f Raytracer::ComputeConductorFresnelReflection(Ray& originalRay, Material& mat, Vec3f& w_o, int recDepth)
{
    if(recDepth <= 0){
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
    if(recDepth <= 0){
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

            float psi1 = GetNormalizedRandom()- 0.5f;
            float psi2 = GetNormalizedRandom()- 0.5f;

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

        // generate two random floats in range [-0.5, 0.5]
        float psi1 = GetNormalizedRandom()- 0.5f;
        float psi2 = GetNormalizedRandom()- 0.5f;

        Vec3f rPrime = makeUnit(r + (u*psi1 + v*psi2)*roughness);
        return rPrime;
    }
    else return r;
}

Vec3f Raytracer::ComputeMirrorReflection(Ray& originalRay, Material& mat, Vec3f& w_o, int recursionDepth)
{
    if (recursionDepth <= 0)
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

Vec3f Raytracer::GetDiffuseReflectanceCoeff(Ray& ray, Shape* shape, Material& mat)
{
    Vec3f reflectance = mat.diffuse;
    if(shape->HasDiffuseTexture())
    {
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
            reflectance = (textureKd + mat.diffuse) / 2.0f;
        }
        else // Replace/overwrite k_d with tex value
        {
            reflectance = textureKd;
        }
    }
    return reflectance;
}
Vec3f Raytracer::GetSpecularReflectanceCoeff(Ray& ray, Shape* shape, Material& mat)
{
    Vec3f reflectance = mat.specular;
    if(shape->HasSpecularTexture())
    {
        // normalize, obtain k_d in [0,1] range.
        Vec3f textureKs;
        if(shape->diffuseTex->IsGenerated()){
            // such as Perlin noise. Sampling process is different.
            auto hp = ray.hitInfo.hitPoint;
            float perlinSample = shape->diffuseTex->GetSampleFromWorldPos(hp.x, hp.y, hp.z);
            textureKs.x = textureKs.y = textureKs.z = perlinSample;
        }
        else
        {
            Vec2f& uv = ray.hitInfo.hitUV;
            textureKs = shape->diffuseTex->GetRGBSample(uv.x, uv.y) / 255.0f;
        }

        if(shape->diffuseTex->operationMode == Texture::OperationMode::Blend)
        {
            //Equally mix the existing material kd value with the tex sample.
            reflectance = (textureKs + mat.diffuse) / 2.0f;
        }
        else // Replace/overwrite k_d with tex value
        {
            reflectance = textureKs;
        }
    }
    return reflectance;
}
Vec3f Raytracer::GetDiffuse(Ray& ray, Shape* shape, Material& mat, Vec3f& w_i, Vec3f& receivedIrradiance)
{
    Vec3f reflectance = GetDiffuseReflectanceCoeff(ray, shape, mat);
    float costheta = std::max(0.0f, dot(w_i, ray.hitInfo.normal));
    return reflectance * receivedIrradiance * costheta;
}

Vec3f Raytracer::GetSpecular(Ray& ray, Shape* shape, Material& mat, Vec3f& w_i, Vec3f& w_o, Vec3f& receivedIrradiance)
{
    Vec3f ks = GetSpecularReflectanceCoeff(ray, shape, mat);

    Vec3f half = (w_i + w_o) / len(w_i + w_o);
    float cosAlpha = std::max(0.0f, dot(ray.hitInfo.normal, half));
    return ks * receivedIrradiance * std::pow(cosAlpha, mat.phong_exponent);
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
        if(scene.materials[mesh->GetMaterial()-1].type == Material::Emissive){
            // skip light objects
            continue;
        }

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
    ray.throughput = original.throughput;

    return ray;
}
Ray Raytracer::GenerateRay(int i, int j, Camera& cam)
{
    Ray ray;
    Vec3f imagePlanePos = cam.GetImagePlanePosition(i,j);

    ray.origin = cam.position;
    ray.throughput = Vec3f{1.0f, 1.0f, 1.0f};

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
    
    ray.hitInfo.hasHit = false;
    ray.hitInfo.minT = INFINITY;
    ray.refractiveIndexOfCurrentMedium = 1.0f;
    ray.motionBlurTime = GetNormalizedRandom();
    
    return ray;
}

Vec3f Raytracer::SampleDirectLighting(Ray& ray, Material& mat, Vec3f& w_o, int lightIDToSkip)
{
    Vec3f color{0,0,0};

    // Repeat for each light source.
    for(int i = 0; i < scene.point_lights.size(); i++)
    {
        // Skip all calculation if in shadow, no color "contribution".
        PointLight& light = scene.point_lights[i];
        if(IsInShadow(ray, light.position)){
            continue;
        }
        Vec3f w_i = makeUnit(light.position - ray.hitInfo.hitPoint);
        float distToLight = len(light.position - ray.hitInfo.hitPoint);
        Vec3f receivedIrradiance = light.intensity / (distToLight * distToLight);

        color = color + Shade(ray, mat, w_i, w_o, receivedIrradiance);
    }    

    for(int i = 0; i < scene.areaLights.size(); i++)
    {
        AreaLight* areaLight = scene.areaLights[i];
        Vec3f lightSamplePos = areaLight->GetSample();
        if(IsInShadow(ray, lightSamplePos)){
            continue;
        }

        Vec3f w_i = lightSamplePos - ray.hitInfo.hitPoint;
        float distToLight = len(w_i);
        float dSqr = distToLight * distToLight;
        w_i = w_i / distToLight; // normalize.

        float lCostheta = dot(areaLight->normal, -w_i);
        if(lCostheta < 0){
            lCostheta = dot(areaLight->normal, w_i);
        }
        Vec3f receivedIrradiance = areaLight->radiance * (areaLight->area * lCostheta / dSqr);

        color = color + Shade(ray, mat, w_i, w_o, receivedIrradiance);
    }        
    for(int i = 0; i < scene.sphericalEnvLights.size(); i++)
    {
        SphericalEnvironmentLight* envLight = scene.sphericalEnvLights[i];
        Vec3f sampleDir = envLight->GetDirection(ray.hitInfo.normal);
        
        // TODO: should we cast Shadow ray?
        // if(IsInShadow(ray, samplePos)){
        //     continue;
        // }
        // TODO: should we cast Shadow ray?

        Vec3f receivedIrradiance = envLight->GetSample(sampleDir);
        Vec3f w_i = ray.hitInfo.normal;
        color = color + Shade(ray, mat, w_i, w_o, receivedIrradiance);
    }    

    for(int i = 0; i < scene.directionalLights.size(); i++)
    {
        DirectionalLight* light = scene.directionalLights[i];
        if(IsInShadowDirectional(ray, light->dir)){
            continue;
        }
        Vec3f w_i = -(light->dir);
        color = color + Shade(ray, mat, w_i, w_o, light->radiance);
    }
    
    for(int i = 0; i < scene.spotLights.size(); i++)
    {
        SpotLight* light = scene.spotLights[i];
        if(IsInShadow(ray, light->pos)){
            continue;
        }
        Vec3f w_i = makeUnit(light->pos - ray.hitInfo.hitPoint);
        Vec3f receivedIrradiance = light->GetIrradiance(ray.hitInfo.hitPoint);
        color = color + Shade(ray, mat, w_i, w_o, receivedIrradiance);
    }

    for(int i = 0; i < scene.meshLights.size(); i++)
    {
        MeshLight* light = scene.meshLights[i];
        if(light->id == lightIDToSkip) continue;

        Vec3f lightSamplePos, lightNormal;
        double weight;
        light->getSample(lightSamplePos, lightNormal, weight);

        if(IsInShadow(ray, lightSamplePos)){
            continue;
        }     
        Vec3f w_i = lightSamplePos - ray.hitInfo.hitPoint;
        float distToLight = len(w_i);
        float dSqr = distToLight * distToLight;
        w_i = w_i / distToLight; // normalize.
        
        if(isnan(w_i.x) || isnan(w_i.y) || isnan(w_i.z)){
            std::cout<<"nan w_i!!" << std::endl;
        }
        float lCostheta = dot(lightNormal, -w_i);
        if(lCostheta < 0){
            lCostheta = dot(lightNormal, w_i);
        }

        Vec3f rad= light->radiance * weight * 2 * M_PI;
        color = color + Shade(ray, mat, w_i, w_o, rad);

        if(isnan(color.x) || isnan(color.y) || isnan(color.z))
        {
            std::cout<<"nan color!!" << std::endl;
            exit(1);
        }
    }

    return color;
}