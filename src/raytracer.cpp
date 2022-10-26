#include "raytracer.hpp"
#include <iostream>

using namespace parser;
using namespace DorkTracer;

Raytracer::Raytracer(Scene& scene){
    this->scene = scene;
}

Vec3i Raytracer::RenderPixel(int i, int j, int camIndex)
{
    return PerPixel(i,j, this->scene.cameras[camIndex]);
}

Vec3i Raytracer::PerPixel(int coordX, int coordY, Camera& cam)
{
    Ray ray = GenerateRay(coordX, coordY, cam);
    
    // Intersect with all Objects in the scene to find the closest intersection point, if any.
    IntersectObjects(ray);
    
    if(ray.hitInfo.hasHit)
    {

       return clamp(PerformShading(ray, cam.position, scene.max_recursion_depth));
    }
    else return scene.background_color;
}

Vec3f Raytracer::PerformShading(Ray& ray, Vec3f& eyePos, int recursionDepth)
{
    Vec3f intersectionPoint = ray.origin + ray.dir * ray.hitInfo.minT;
    // Actually perform shading, depending on Material type.

    Material& mat = scene.materials[ray.hitInfo.matId-1];
    Vec3f color = GetAmbient(mat.ambient, scene.ambient_light);
    Vec3f w_o = makeUnit(eyePos - intersectionPoint);

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

    if (mat.type == Default)
    {   // we are done with shading, default/standard only have Diffuse Specular and Ambient components.
        return color; 
    }
    else if(mat.type == Mirror)
    {
        // Compute radiance along the ideal reflection ray
        color = color + ComputeMirrorReflection(mat.mirror, w_o, ray.hitInfo.normal, intersectionPoint, recursionDepth);
    }
    else if(mat.type == Dielectric)
    {
        // Both reflection & transmission
        float n1 = ray.refractiveIndexOfCurrentMedium;
        color = color + ComputeDielectricFresnelReflectionAndRefraction(mat, intersectionPoint, w_o, ray.hitInfo.normal, n1, mat.refractiveIndex, recursionDepth);
    }
    else if(mat.type == Conductor)
    {
        // Reflection & absorption, no tranmission.
    }

    return color;
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
        Vec3f w_r =  Reflect(normal, w_o);

        // Do reflection?
        Ray ray;
        ray.origin = x + normal * scene.shadow_ray_epsilon;
        ray.dir = w_r;
        ray.hitInfo.hasHit = false;
        ray.hitInfo.minT = 999999;
        ray.refractiveIndexOfCurrentMedium = raysN; // it didnt change medium, reflected back into the one its coming from.

        IntersectObjects(ray);
        if(ray.hitInfo.hasHit)
        {
            return PerformShading(ray, ray.origin, recDepth-1);
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
            Vec3f w_reflected = Reflect(normal, w_o);

            // Again do reflection here
            Ray ray;
            ray.origin = x + normal * scene.shadow_ray_epsilon;
            ray.dir = makeUnit(w_reflected);
            ray.hitInfo.hasHit = false;
            ray.hitInfo.minT = 999999;
            ray.refractiveIndexOfCurrentMedium = raysN; // it didnt change medium, reflected back into the one its coming from.

            IntersectObjects(ray);
            if(ray.hitInfo.hasHit)
            {
                // no attenuation on reflect case, since we are not INSIDE an object?
                // what about total internal reflection though?
                reflectedRaysColor = PerformShading(ray, ray.origin, recDepth-1);
            }
            else reflectedRaysColor = Vec3f{0,0,0};
        }
       
        // BUT also do refraction, 1 ray split into 2 rays in total.
        // wt = (d + ncosTheta)(n1/n2) - ncosPhi, page 16 of slides.
        Vec3f refractedRaysColor;
        {
            Vec3f w_refracted = (d + normal * cosTheta) * r - normal * cosPhi;
            
            Ray ray;
            ray.dir = makeUnit(w_refracted);
            ray.origin = x + (-normal) * scene.shadow_ray_epsilon;
            ray.hitInfo.hasHit = false;
            ray.hitInfo.minT = 999999;
            ray.refractiveIndexOfCurrentMedium = objN; // it changed medium, transmitted to other one.

            // Maybe only intersect with current object if its entering.
            // cuz it has to intersect without leaving if  this is a non-infinite object, or a plane.

            IntersectObjects(ray);
            
            if(ray.hitInfo.hasHit)
            {
                // todo attenuation??
                // Vec3f newHitPos = ray.origin + ray.dir * ray.hitInfo.minT;
                float travelledDist = ray.hitInfo.minT;
                // Vec3f attenuationFactor =   L_0);
                // refractedRaysColor = attenuationFactor * PerformShading(ray, ray.origin, recDepth-1);
                refractedRaysColor = BeersLaw(travelledDist, mat.absorptionCoefficient, PerformShading(ray, ray.origin, recDepth-1));

            }
            else{
                // std::cout << "refracted ray didnt hit anything"<< std::endl;
                refractedRaysColor = Vec3f{0,0,0}; 
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
Vec3f Raytracer::Reflect(Vec3f& normal, Vec3f& w_o)
{
    return makeUnit((normal * 2.0f * dot(normal, w_o)) - w_o);
}

Vec3f Raytracer::ComputeMirrorReflection(Vec3f reflectance, Vec3f& w_o, Vec3f& normal, Vec3f& intersectionPoint, int recursionDepth)
{
    if (recursionDepth == 0)
    {
        return Vec3f{0,0,0};
    }
    else 
    {
        Vec3f w_r =  Reflect(normal, w_o);

        Ray ray;
        ray.origin = intersectionPoint + normal * scene.shadow_ray_epsilon;
        ray.dir = w_r;
        ray.hitInfo.hasHit = false;
        ray.hitInfo.minT = 999999;
        
        IntersectObjects(ray);
        if(ray.hitInfo.hasHit)
        {
            return reflectance * PerformShading(ray, ray.origin, recursionDepth-1);
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
        Mesh& mesh = scene.meshes[i];
        for(int j = 0; j < mesh.faces.size(); j++){

            if(mesh.useOwnVertices){
                IntersectFace(shadowRay, mesh.faces[j], mesh.vertices, mesh.material_id);
            }
            else{
                IntersectFace(shadowRay, mesh.faces[j], scene.vertex_data, mesh.material_id);
            }

            if(shadowRay.hitInfo.hasHit && shadowRay.hitInfo.minT < lightSourceT)
            {
                // valid hit, in shadow.
                return true;
            }
        }
    }
    
    // Intersect with all Triangles
    for(int i = 0; i < scene.triangles.size(); i++){
        Triangle& tri = scene.triangles[i];
        IntersectFace(shadowRay, tri.indices, scene.vertex_data, tri.material_id);
        if(shadowRay.hitInfo.hasHit && shadowRay.hitInfo.minT < lightSourceT)
        {
            // valid hit, in shadow.
            return true;
        }
    }

    // Intersect with all Spheres
    for(int i = 0; i < scene.spheres.size(); i++){
        Sphere& s = scene.spheres[i];
        IntersectSphere(shadowRay, s);
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
    bool bfcEnabled = false;
    // TODO: BFC should be disabled for refractive materials, they have both "faces" by default.

    for(int i = 0; i < scene.meshes.size(); i++){
        Mesh& mesh = scene.meshes[i];
        for(int j = 0; j < mesh.faces.size(); j++){
            if(bfcEnabled && IsBackface(mesh.faces[j], ray.dir)){
                continue;
            }
            
            // TODO: Improve this by introducing proper OOP.
            Face& face = mesh.faces[j];
            if(mesh.useOwnVertices){
                IntersectFace(ray, face, mesh.vertices, mesh.material_id);
            }
            else{
                IntersectFace(ray, face, scene.vertex_data, mesh.material_id);
            }
        }
    }
    
    // Intersect with all Triangles
    for(int i = 0; i < scene.triangles.size(); i++){
        Triangle& tri = scene.triangles[i];
        if(bfcEnabled && IsBackface(tri.indices, ray.dir))
        {
            continue;
        }
        IntersectFace(ray, tri.indices, scene.vertex_data, tri.material_id);
    }

    // Intersect with all Spheres
    for(int i = 0; i < scene.spheres.size(); i++){
        Sphere& s = scene.spheres[i];
        IntersectSphere(ray, s);
    }
}
bool Raytracer::IsBackface(Face& face, Vec3f& rayDir)
{
    return dot(face.n, rayDir) > 0.0f;
}
void Raytracer::IntersectFace(Ray& ray, Face& face, std::vector<Vec3f>& vertices, int matId)
{
    Vec3f& v0 = vertices[face.v0_id-1];
    Vec3f& v1 = vertices[face.v1_id-1];
    Vec3f& v2 = vertices[face.v2_id-1];
    float newT  = 99999;
    bool hasIntersected = DoesIntersectTriangle(ray, v0, v1, v2, newT);

    // If its the closest intersection so far, update ray HitInfo.
    if(hasIntersected && newT < ray.hitInfo.minT){
        ray.hitInfo.minT = newT;
        ray.hitInfo.hasHit = true;
        ray.hitInfo.normal = face.n;
        ray.hitInfo.matId = matId;
    }
}
bool Raytracer::DoesIntersectTriangle(Ray& ray, Vec3f& v0, Vec3f& v1, Vec3f& v2, float& t)
{
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
    t = determinant(matrixT) / detA;
    return t > 0.0f;
}
void Raytracer::IntersectSphere(Ray& r, Sphere& s)
{
    Vec3f center = scene.vertex_data[s.center_vertex_id - 1];
    Vec3f oc = r.origin - center;

    float t,t1, t2;

    float c = dot(oc,oc) - (s.radius * s.radius);
    float b =  2 * dot(r.dir, oc);
    float a = dot(r.dir, r.dir);
    float delta = b*b - (4 * a * c);

    if(delta < 0.0) return;
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
        // todo: consider ignoring negative t before selecting the min.

        if(t < r.hitInfo.minT && t > 0.0f){
            r.hitInfo.minT = t;
            r.hitInfo.matId = s.material_id;
            r.hitInfo.hasHit = true;
            // calculate sphere normal at hit point
            Vec3f hitPoint = r.origin + r.dir * t;
            r.hitInfo.normal = makeUnit(hitPoint - center);
        }
    }
}
Ray  Raytracer::GenerateRay(int i, int j, Camera& cam)
{
    Ray ray;
    float su, sv;
    Vec3f m, q, s;

    // NearPlane: coords of image plane with Left, Right, Bottom, Top floats respectively.
    float left,right,bottom,top;
    float nx = cam.image_width;
    float ny = cam.image_height;
    float dist = cam.near_distance;

    if(cam.isLookAt){
        float aspect = nx/ny;
        
        top = dist * std::tan((cam.fovY *(M_PI / 180.0f) / 2.0f));
        right = top * aspect;

        bottom = -top;
        left = -right;

        cam.gaze = makeUnit(cam.gazePoint - cam.position);
    }
    else{
        left = cam.near_plane.x;
        right = cam.near_plane.y;
        bottom = cam.near_plane.z;
        top = cam.near_plane.w;
    }

    su = (i + 0.5) * (right - left) / nx;
    sv = (j + 0.5) * (top - bottom) / ny;

    Vec3f e = cam.position;

    // Up = v,  Gaze = −w, u = v ×w
    Vec3f v = cam.up;
    Vec3f w = -cam.gaze;

    Vec3f u = cross(v, w);

    m = e + cam.gaze * dist;
    q = m + u * left+ v * top;
    s = q + u * su + v *-sv;

    ray.origin = e;
    ray.dir = makeUnit(s - e);
    
    ray.hitInfo.hasHit = false;
    ray.hitInfo.minT = 99999999;
    ray.refractiveIndexOfCurrentMedium = 1.0f;
    return ray;

}