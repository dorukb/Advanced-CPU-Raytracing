#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <limits>

#include "happly.h"
#include "scene.h"
#include "matrix.hpp"
#include "instancedMesh.hpp"
#include "HDRImage.h"
#include "LDRImage.h"
#include "spotLight.h"
#include "meshLight.h"

#include "brdfModifiedBlinnPhong.h"
#include "brdfBlinnPhong.h"
#include "brdfPhong.h"
#include "brdfModifiedPhong.h"
#include "brdfTorranceSparrow.h"
#include "brdf.h"


using namespace DorkTracer;

void DorkTracer::Scene::loadFromXml(const std::string &filepath)
{
    tinyxml2::XMLDocument file;
    std::stringstream stream;
    tinyxml2::XMLElement *child, *element;

    auto res = file.LoadFile(filepath.c_str());
    if (res)
    {
        throw std::runtime_error("Error: The xml file cannot be loaded.");
    }

    auto root = file.FirstChild();
    if (!root)
    {
        throw std::runtime_error("Error: Root is not found.");
    }

    this->isMotionBlurEnabled = false;

    //Get BackgroundColor
    background_color = Vec3i(0,0,0);
    element = root->FirstChildElement("BackgroundColor");
    if (element)
    {
        stream << element->GetText() << std::endl;
        stream >> background_color.x >> background_color.y >> background_color.z;
    }

    //Get ShadowRayEpsilon
    element = root->FirstChildElement("ShadowRayEpsilon");
    if (element)
    {
        stream << element->GetText() << std::endl;
        stream >> Scene::shadow_ray_epsilon;
    }

    //Get MaxRecursionDepth
    max_recursion_depth = 0;
    element = root->FirstChildElement("MaxRecursionDepth");
    if (element)
    {
        stream << element->GetText() << std::endl;
        stream >> max_recursion_depth;
    }

    //Get Cameras
    parseCameras(root);

    //Get Lights
    parseLights(root);

    // Get BRDFs, must perform before parsing Materials.
    parseBRDFs(root);

    // Get Materials
    parseMaterials(root);

    // Get Textures
    this->bgTexture = nullptr;
    element = root->FirstChildElement("Textures");
    if(element)
    {
        // Load images if any
        element = element->FirstChildElement("Images");
        if(element){
            element = element->FirstChildElement("Image");
            while (element){
                int id;
                stream << element->Attribute("id") << std::endl;
                stream >> id;

                std::string filename;
                stream << element->GetText() << std::endl;
                stream >> filename;

                Image* image = nullptr;
                if(filename.find(".exr") != std::string::npos)
                {
                    // EXR/ HDR image
                    std::cout<<"Creating HDR image from:" << filename << std::endl;
                    image = new HDRImage(std::string("inputs/" + filename), id);
                }
                else{
                    image = new LDRImage(std::string("inputs/" + filename), id);
                }
                this->images.push_back(image);

                
                element = element->NextSiblingElement("Image");
            }
        }

        // Load texture maps if any
        element = root->FirstChildElement("Textures");
        element = element->FirstChildElement("TextureMap");

        while(element)
        {
            int textureId;
            stream << element->Attribute("id") << std::endl;
            stream >> textureId;

            std::string textureType;
            stream << element->Attribute("type") << std::endl;
            stream >> textureType;

            // The only property all textures share is the decal mode or apply mode.
            child = element->FirstChildElement("DecalMode");
            std::string decalMode;
            stream << child->GetText() << std::endl;
            stream >> decalMode;
            
            if(textureType == "image")
            {
                int imageId;
                child = element->FirstChildElement("ImageId");
                stream << child->GetText() << std::endl;
                stream >> imageId;
                
                child = element->FirstChildElement("Interpolation");
                std::string interpMode = "nearest";
                if(child){
                    stream << child->GetText() << std::endl;
                    stream >> interpMode;
                }

                float normalizer = 255.0f;
                child = element->FirstChildElement("Normalizer");
                if(child){
                    stream << child->GetText() << std::endl;
                    stream >> normalizer;
                }

                float sampleMultiplier = 1.0f;
                child = element->FirstChildElement("BumpFactor");
                if(child){
                    stream << child->GetText() << std::endl;
                    stream >> sampleMultiplier;
                }

                Image* img = nullptr;
                for(int i = 0; i < this->images.size(); i++){
                    if(images[i]->id == imageId){
                        img = images[i];
                        break;
                    }
                }
                if(img == nullptr){
                    std::cout << "Couldnt get image with id: " + imageId << std::endl; 
                }

                ImageTexture* tex = new ImageTexture(textureId, normalizer, sampleMultiplier, decalMode, interpMode, img);
                this->textures.push_back(tex);

                if(decalMode == "replace_background")
                {
                    this->bgTexture = tex;
                    std::cout << "set the background texture." << std::endl;
                }
            }
            else if(textureType == "perlin")
            {
                std::string noiseConversion = "linear";
                child = element->FirstChildElement("NoiseConversion");
                if(child){
                    stream << child->GetText() << std::endl;
                    stream >> noiseConversion;
                }

                float noiseScale = 1.0f;
                child = element->FirstChildElement("NoiseScale");
                if(child){
                    stream << child->GetText() << std::endl;
                    stream >> noiseScale;
                }   
                
                float bumpFactor = 1.0f;
                child = element->FirstChildElement("BumpFactor");
                if(child){
                    stream << child->GetText() << std::endl;
                    stream >> bumpFactor;
                }

                PerlinTexture* tex = new PerlinTexture(textureId, decalMode, noiseScale, noiseConversion);
                tex->bumpFactor = bumpFactor;
                this->textures.push_back(tex);

                if(decalMode == "replace_background")
                {
                    this->bgTexture = tex;
                    std::cout << "set the background texture." << std::endl;
                }
            }
            else if(textureType == "checkerboard")
            {
                // TODO: implement procedural checkboard texture.
                std::cout << "procedural checkerboard texture is not implemented yet." << std::endl;
            }

            element = element->NextSiblingElement("TextureMap");
        }
    }


    // Parse Environment lights after Images are parsed, to obtain the image pointer.
    element = root->FirstChildElement("Lights");
    if(element)
    {
        element = element->FirstChildElement("SphericalDirectionalLight");
        while (element)
        {
            int id, imageId;
            
            stream << element->Attribute("id") << std::endl;
            stream >> id;

            child = element->FirstChildElement("ImageId");
            stream << child->GetText() << std::endl;
            stream >> imageId;
            
            Image* img = nullptr;
            for(int i = 0; i < this->images.size(); i++){
                if(images[i]->id == imageId){
                    img = images[i];
                    break;
                }
            }
            if(img == nullptr){
                std::cout<<"couldnt get img with id: " <<imageId <<std::endl;
            }
            sphericalEnvLights.push_back(new SphericalEnvironmentLight(id, img));
            element = element->NextSiblingElement("SphericalDirectionalLight");
        }
        stream.clear();
    }
    

    // Get VertexData
    element = root->FirstChildElement("VertexData");
    if(element)
    {
        stream << element->GetText() << std::endl;
        Vec3f vertex;
        while (!(stream >> vertex.x).eof())
        {
            stream >> vertex.y >> vertex.z;
            vertex_data.push_back(vertex);
        }
    }
   
    stream.clear();

    // Get TexCoordData (u,v)
    element = root->FirstChildElement("TexCoordData");
    if(element)
    {
        stream << element->GetText() << std::endl;
        Vec2f texCoord;
        while (!(stream >> texCoord.x).eof())
        {
            stream >> texCoord.y;
            texCoords.push_back(texCoord);
        }
        stream.clear();
    }

    // Get Transformations
    element = root->FirstChildElement("Transformations");
    if(element)
    {
        auto child = element->FirstChildElement("Translation");
        while(child)
        {
            std::string idStr = child->Attribute("id");
            int id = stoi(idStr);

            // std::string text = child->GetText();
            stream << child->GetText() << std::endl;
            Vec3f translation;
            stream >> translation.x >> translation.y >> translation.z;
            this->translations.push_back(translation);

            child = child->NextSiblingElement("Translation");
        }

        child = element->FirstChildElement("Scaling");
        while(child)
        {
            std::string idStr = child->Attribute("id");
            int id = stoi(idStr);

            // std::string text = child->GetText();
            stream << child->GetText() << std::endl;
            Vec3f scaling;
            stream >> scaling.x >> scaling.y >> scaling.z;
            this->scalings.push_back(scaling);
            
            child = child->NextSiblingElement("Scaling");
        }

        child = element->FirstChildElement("Rotation");
        while(child)
        {
            std::string idStr = child->Attribute("id");
            int id = stoi(idStr);

            // std::string text = child->GetText();
            stream << child->GetText() << std::endl;

            // rot.w = angle, rot.xyz = axis
            Vec4f rot;
            stream >> rot.w >> rot.x >> rot.y >> rot.z;
            this->rotations.push_back(rot);
            
            child = child->NextSiblingElement("Rotation");
        }

    }
    stream.clear();

    // Get Meshes
    parseMeshes(root, "Mesh");
    parseMeshes(root, "LightMesh");


    // Parse MeshInstances
    element = root->FirstChildElement("Objects");
    element = element->FirstChildElement("MeshInstance");

    while(element)
    {
        bool resetTransform = false;
        if(element->Attribute("resetTransform") != NULL){
            std::string transformText = element->Attribute("resetTransform");
            if(transformText == "true"){
                resetTransform = true;
            }
        }

        int ownId, baseMeshId;
        stream << element->Attribute("id") << std::endl;
        stream >> ownId;
        stream << element->Attribute("baseMeshId") << std::endl;
        stream >> baseMeshId;

        // If instancing is not nested, parentMesh = baseMesh.
        // but if instance has another instance as the parent, then baseMesh will be different.
        DorkTracer::Mesh* parentMesh = nullptr;
        DorkTracer::Mesh* baseMesh = nullptr;

        for(int i = 0; i < meshes.size(); i++){
            if(meshes[i]->id == baseMeshId){
                parentMesh = (Mesh*) meshes[i];
            }
        }

        baseMesh = parentMesh;
        while(baseMesh->isInstance){
            baseMesh = (*(InstancedMesh*)baseMesh).baseMesh;
        }

        DorkTracer::InstancedMesh* meshInstance = new InstancedMesh(baseMesh);
        meshInstance->id = ownId;
        meshInstance->isInstance = true;

        // Setup textures
        child = element->FirstChildElement("Textures");
        if(child){
            std::string texturesInp = child->GetText();
            texturesInp += " ";
            SetupTextures(meshInstance, texturesInp);
        }

        // Get material
        child = element->FirstChildElement("Material");
        if(child != nullptr){
            stream << child->GetText() << std::endl;
            int matId = 0;
            stream >> matId;
            meshInstance->SetMaterial(matId);
        }
        else{
            meshInstance->SetMaterial(baseMesh->GetMaterial());
        }  
        // parse motion blur
        // <MotionBlur>0 0 4</MotionBlur>
        child = element->FirstChildElement("MotionBlur");
        if(child != nullptr)
        {
            stream << child->GetText() << std::endl;
            stream >> meshInstance->motionBlurVector.x >> meshInstance->motionBlurVector.y >> meshInstance->motionBlurVector.z;

            // enable motior blur globally even if only one object is using it.
            this->isMotionBlurEnabled = true; 
            meshInstance->hasMotionBlur = true;
        }
        else
        {
            meshInstance->hasMotionBlur = false;
            meshInstance->motionBlurVector = Vec3f(0.0f,0.0f,0.0f);
        }

        // Get transformations
        child = element->FirstChildElement("Transformations");
        meshInstance->transform.MakeIdentity();
        meshInstance->inverseTransform.MakeIdentity();
        if(child != nullptr)
        {
            // Compute Model, invModel and invTransModel matrices from given transformations.
            std::string input = child->GetText();
            computeTransform(meshInstance, input);

            if(resetTransform == false){
                
                // also need to factor in transform matrices of our baseMesh.
                meshInstance->transform = meshInstance->transform * parentMesh->transform;

                // order must be reversed here.
                meshInstance->inverseTransform = parentMesh->inverseTransform * meshInstance->inverseTransform;
                meshInstance->inverseTransposeTransform = Matrix::GetTranspose(meshInstance->inverseTransform);
            }
        }

        // compute and set transformed bbox.
        meshInstance->bbox = transformBoundingBox(baseMesh->bbox, meshInstance->transform);
        
        meshes.push_back(meshInstance);
        element = element->NextSiblingElement("MeshInstance");
    }


    //Get Triangles, they can be represented with Mesh structure as well.
    element = root->FirstChildElement("Objects");
    element = element->FirstChildElement("Triangle");

    while (element)
    {
        DorkTracer::Mesh* mesh = new DorkTracer::Mesh(this->vertex_data, this->texCoords);

        // Compute Model, invModel and invTransModel matrices from given transformations.
        child = element->FirstChildElement("Transformations");
        mesh->transform.MakeIdentity();
        if(child != nullptr)
        {
            std::string input = child->GetText();
            computeTransform(mesh, input);
        }
        else{
            mesh->inverseTransform.MakeIdentity();
            mesh->inverseTransposeTransform.MakeIdentity();
        }
        
        // Setup textures
        child = element->FirstChildElement("Textures");
        if(child){
            std::string texturesInp = child->GetText();
            texturesInp += " ";
            SetupTextures(mesh, texturesInp);
        }


        child = element->FirstChildElement("Material");
        stream << child->GetText() << std::endl;
        int matId = 0;
        stream >> matId;
        mesh->SetMaterial(matId);

        child = element->FirstChildElement("Indices");
        stream << child->GetText() << std::endl;
        // stream >> triangle.indices.v0_id >> triangle.indices.v1_id >> triangle.indices.v2_id;
        
        DorkTracer::Face face;
        
        stream >> face.v0_id >> face.v1_id >> face.v2_id;

        // Normal, Centroid and BoundingBox computation per face.
        computeFaceProperties(face, mesh);
        mesh->faces.push_back(face);

        // Has only one face, so the bounding boxes are identical.
        mesh->bbox = face.bbox;
        mesh->ConstructBVH();

        meshes.push_back(mesh);
        element = element->NextSiblingElement("Triangle");
    }

    //Get Spheres
    element = root->FirstChildElement("Objects");
    element = element->FirstChildElement("Sphere");
    while (element)
    {    
        Sphere* sphere = new Sphere(this->vertex_data);

        // Compute Model, invModel and invTransModel matrices from given transformations.
        child = element->FirstChildElement("Transformations");
        sphere->transform.MakeIdentity();
        if(child != nullptr)
        {
            std::string input = child->GetText();
            computeTransform(sphere, input);
        }
        else{
            sphere->inverseTransform.MakeIdentity();
            sphere->inverseTransposeTransform.MakeIdentity();
        }

        // Setup textures
        child = element->FirstChildElement("Textures");
        if(child){
            std::string texturesInp = child->GetText();
            texturesInp += " ";
            SetupTextures(sphere, texturesInp);
        }

        child = element->FirstChildElement("Material");
        stream << child->GetText() << std::endl;
        stream >> sphere->material_id;

        child = element->FirstChildElement("Center");
        stream << child->GetText() << std::endl;
        stream >> sphere->center_vertex_id;

        child = element->FirstChildElement("Radius");
        stream << child->GetText() << std::endl;
        stream >> sphere->radius;


        // parse motion blur
        child = element->FirstChildElement("MotionBlur");
        if(child != nullptr)
        {
            stream << child->GetText() << std::endl;
            stream >> sphere->motionBlurVector.x >> sphere->motionBlurVector.y >> sphere->motionBlurVector.z;

            sphere->hasMotionBlur = true;
            // enable motior blur globally even if only one object is using it.
            this->isMotionBlurEnabled = true; 
        }
        else{
            sphere->hasMotionBlur = false;
            sphere->motionBlurVector = Vec3f(0.0f,0.0f,0.0f);
        }


        spheres.push_back(sphere);
        element = element->NextSiblingElement("Sphere");
    }

    // std::cout <<"parsing succesfully completed." << std::endl;
}

void DorkTracer::Scene::computeFaceProperties(DorkTracer::Face& face, Mesh* mesh)
{
    computeFaceCenter(face, mesh);
    computeFaceNormal(face, mesh);
    computeFaceBoundingBox(face, mesh);
    computeFaceArea(face, mesh);
}
void DorkTracer::Scene::computeFaceCenter(DorkTracer::Face& face, Mesh* mesh)
{
    Vec3f& a = mesh->GetVertex(face.v0_id);
    Vec3f& b = mesh->GetVertex(face.v1_id);
    Vec3f& c = mesh->GetVertex(face.v2_id);

    face.center = (a+b+c) /3.0f;
}
void DorkTracer::Scene::computeFaceArea(DorkTracer::Face& face, Mesh* mesh)
{
    Vec3f& a = mesh->GetVertex(face.v0_id);
    Vec3f& b = mesh->GetVertex(face.v1_id);
    Vec3f& c = mesh->GetVertex(face.v2_id);

    // Heron formula for area of triangle with length of all edges known.
    double e1 = len(a-b);
    double e2 = len(a-c);
    double e3 = len(b-c);
    double s = (e1+e2+e3)/2.0f;
    double area = std::sqrt(s * (s-e1) * (s-e2) * (s-e3));

    face.area = area;
    mesh->surfaceArea += area;
    // std::cout <<"a: " << area << " ";
}

void DorkTracer::Scene::SetupTextures(DorkTracer::Shape* shape, std::string& texIds)
{
    size_t last = 0; 
    size_t next = 0;
    while ((next = texIds.find(" ", last)) != std::string::npos) 
    {   
        auto idStr = texIds.substr(last, next-last);
        int id = stoi(idStr);

        Texture* texture = nullptr;
        for(int i = 0; i < this->textures.size(); i++){
            if(textures[i]->id == id){
                texture = textures[i];
                break;
            }
        }
        if(texture == nullptr) break;

        switch(texture->type){
            case Texture::Textures::Diffuse:
                shape->diffuseTex = texture;
                break;
            case Texture::Textures::Specular:
                shape->specularTex = texture;
                break;
            case Texture::Textures::ReplaceAll:
                shape->replaceAll = texture;
                break;
            case Texture::Textures::Normal:
                shape->normalMap = texture;
                break;
            case Texture::Textures::Bump:
                shape->bumpMap = texture;
                break;
        };

        last = next + 1; 
    }
}
void DorkTracer::Scene::computeTransform(DorkTracer::Shape* mesh, std::string transformationInput)
{
    // sample input: <Transformations>s2 s3 r2 t2</Transformations>
    // param: "s2 s3 r2 t2"
            // s:scale, r:rotation, t:translation, numbers are the indices to look up.
    std::string str = transformationInput;
    int idx = 0;

    std::vector<Matrix> invTransformations;
    while(idx < str.size()-1)
    {
        if(str[idx] == 'r'){
            int id = int(str[idx+1]-'0');
            Vec4f rotation = this->rotations[id-1];
            float angle = rotation.w * (M_PI / 180.0f);

            // TODO: Support rotation around arbitray axis.
            Matrix rot;
            Matrix invRot;

            // std::cout<<"r"<<id<<" applied."<<std::endl;
            if(rotation.x >= 0.99 && rotation.y <= 0.001 && rotation.z <= 0.0001){
                rot = Matrix::GetRotationAroundX(angle);
                invRot = Matrix::GetRotationAroundX(-angle);
            }
            if(rotation.y >= 0.99 && rotation.x <= 0.001 && rotation.z <= 0.0001){
                rot = Matrix::GetRotationAroundY(angle);
                invRot = Matrix::GetRotationAroundY(-angle);
            }
            if(rotation.z >= 0.99 && rotation.x <= 0.001 && rotation.y <= 0.0001){
                rot = Matrix::GetRotationAroundZ(angle);
                invRot = Matrix::GetRotationAroundZ(-angle);
            }
            invTransformations.push_back(invRot);
            mesh->transform = rot * mesh->transform;
        }
        else if(str[idx] == 't'){

            int id = int(str[idx+1]-'0');
            Vec3f translation = this->translations[id-1];
            
            Matrix t = Matrix::GetTranslation(translation.x, translation.y , translation.z);
            Matrix invT = Matrix::GetTranslation(-translation.x, -translation.y , -translation.z);
            
            invTransformations.push_back(invT);
            mesh->transform = t * mesh->transform;
        }
        else if(str[idx] == 's'){
            int id = int(str[idx+1]-'0');
            Vec3f scale = this->scalings[id-1];

            // std::cout<<"s"<<id<<" applying."<<std::endl;
            Matrix s = Matrix::GetScale(scale.x, scale.y , scale.z);
            Matrix invS = Matrix::GetScale(1.0f/scale.x, 1.0f/scale.y , 1.0f/scale.z);

            invTransformations.push_back(invS);
            mesh->transform = s * mesh->transform;
        }
        idx += 3;
    }

    mesh->inverseTransform.MakeIdentity();
    for(int i = 0; i < invTransformations.size(); i++){
        Matrix trans = invTransformations[i];
        mesh->inverseTransform = mesh->inverseTransform * trans;

    }
    // Matrix modelMatrix = mesh->transform;
    // Matrix invModel = mesh->inverseTransform;
    mesh->inverseTransposeTransform = Matrix::GetTranspose(mesh->inverseTransform);

    // Matrix idt = modelMatrix * invModel;
}
void DorkTracer::Scene::computeFaceNormal(DorkTracer::Face& face, Mesh* mesh)
{
    Vec3f& a = mesh->GetVertex(face.v0_id);
    Vec3f& b = mesh->GetVertex(face.v1_id);
    Vec3f& c = mesh->GetVertex(face.v2_id);
    Vec3f ab = b - a;
    Vec3f ac = c - a;

    face.n = makeUnit(cross(ab, ac));
}
void DorkTracer::Scene::computeFaceBoundingBox(DorkTracer::Face& face, Mesh* mesh)
{
    Vec3f& a = mesh->GetVertex(face.v0_id);
    Vec3f& b = mesh->GetVertex(face.v1_id);
    Vec3f& c = mesh->GetVertex(face.v2_id);

    face.bbox.minCorner.x = std::min(std::min(a.x, b.x), c.x);
    face.bbox.minCorner.y = std::min(std::min(a.y, b.y), c.y);
    face.bbox.minCorner.z = std::min(std::min(a.z, b.z), c.z);

    face.bbox.maxCorner.x = std::max(std::max(a.x, b.x), c.x);    
    face.bbox.maxCorner.y = std::max(std::max(a.y, b.y), c.y);
    face.bbox.maxCorner.z = std::max(std::max(a.z, b.z), c.z);
}

DorkTracer::BoundingBox DorkTracer::Scene::transformBoundingBox(DorkTracer::BoundingBox original, DorkTracer::Matrix& transform)
{

    // Transformed box is no longer axis aligned, so can not correctly decide on min max using just 2 corners as before.
    // Calculate all 8 corners of the original, transform each,
    // create another aligned bounding box out of min max.
    BoundingBox res;
    res.maxCorner.x = res.maxCorner.y = res.maxCorner.z = -INFINITY;
    res.minCorner.x = res.minCorner.y = res.minCorner.z = INFINITY;

    std::vector<Vec3f> corners;
    corners.push_back(original.maxCorner);
    corners.push_back(original.minCorner);

    Vec3f extents = original.maxCorner - original.minCorner;
    Vec3f corner = original.maxCorner;
    corner.x -= extents.x;
    corners.push_back(corner);

    corner = original.maxCorner;
    corner.y -= extents.y;
    corners.push_back(corner);

    corner = original.maxCorner;
    corner.z -= extents.z;
    corners.push_back(corner);
    
    corner = original.maxCorner;
    corner.x -= extents.x;
    corner.y -= extents.y;
    corners.push_back(corner);
    
    corner = original.maxCorner;
    corner.x -= extents.x;
    corner.z -= extents.z;
    corners.push_back(corner);

    corner = original.maxCorner;
    corner.y -= extents.y;
    corner.z -= extents.z;
    corners.push_back(corner);

    for(int i = 0; i < corners.size(); i++)
    {
        corners[i] = Matrix::ApplyTransformToPoint(transform, corners[i]);
        
        res.maxCorner.x = std::max(corners[i].x, res.maxCorner.x);
        res.maxCorner.y = std::max(corners[i].y, res.maxCorner.y); 
        res.maxCorner.z = std::max(corners[i].z, res.maxCorner.z);

        res.minCorner.x = std::min(corners[i].x, res.minCorner.x);
        res.minCorner.y = std::min(corners[i].y, res.minCorner.y);
        res.minCorner.z = std::min(corners[i].z, res.minCorner.z);
    }   
    
    return res;
}
DorkTracer::Face DorkTracer::Scene::createFace(int v0idx, int v1idx, int v2idx, DorkTracer::Mesh* mesh)
{
    DorkTracer::Face face;
    face.v0_id = v0idx;
    face.v1_id = v1idx;
    face.v2_id = v2idx;
    computeFaceProperties(face, mesh);
    std::cout << " total area of mesh: " << mesh->surfaceArea << std::endl;

    return face;
}
void Scene::updateBBox(BoundingBox* bbox, Face& face)
{
    // calculate bbox of the owner mesh here
    bbox->minCorner.x = std::min(face.bbox.minCorner.x,  bbox->minCorner.x);
    bbox->minCorner.y = std::min(face.bbox.minCorner.y,  bbox->minCorner.y);
    bbox->minCorner.z = std::min(face.bbox.minCorner.z,  bbox->minCorner.z);

    bbox->maxCorner.x = std::max(face.bbox.maxCorner.x, bbox->maxCorner.x);
    bbox->maxCorner.y = std::max(face.bbox.maxCorner.y, bbox->maxCorner.y);
    bbox->maxCorner.z = std::max(face.bbox.maxCorner.z, bbox->maxCorner.z);
}
void DorkTracer::Scene::parseTonemapper(tinyxml2::XMLElement* cameraElm, DorkTracer::Camera* cam)
{ 
    // - TMO tonemapping operator type
    // - TMOOptions key value & burn percentage
    // - Saturation
    // - gamma
    auto child = cameraElm->FirstChildElement("Tonemap");
    if(child != nullptr)
    {
        std::stringstream stream;
        auto childElmPtr = child->FirstChildElement("TMO");
        std::string opType = "Photographic";
        if(childElmPtr != nullptr){
            stream << childElmPtr->GetText() << std::endl;
            stream >> opType;
        }

        childElmPtr = child->FirstChildElement("TMOOptions");
        float keyValue = 0.18f;
        float burnPerct = 1.0f;
        if(childElmPtr != nullptr){
            stream << childElmPtr->GetText() << std::endl;
            stream >> keyValue >> burnPerct;
        }

        childElmPtr = child->FirstChildElement("Saturation");
        float saturation = 1.0f;
        if(childElmPtr != nullptr){
            stream << childElmPtr->GetText() << std::endl;
            stream >> saturation;
        }

        childElmPtr = child->FirstChildElement("Gamma");
        float gamma = 2.2f;
        if(childElmPtr != nullptr){
            stream << childElmPtr->GetText() << std::endl;
            stream >> gamma;
        }
        cam->SetTonemapper(new Tonemapper(opType, keyValue,burnPerct,saturation,gamma));
        
    }
}
void DorkTracer::Scene::parseBRDFs(tinyxml2::XMLNode* root)
{ 
    // Get BRDFs
    tinyxml2::XMLElement* element = root->FirstChildElement("BRDFs");
    if(element)
    {
        // Parse ModifiedBlinnPhongBRDF
        element = element->FirstChildElement("ModifiedBlinnPhong");
        std::stringstream stream;
        while(element)
        {
            int id = -1;
            stream << element->Attribute("id") << std::endl;
            stream >> id;

            bool normalized = false;
            if(element->Attribute("normalized", "true")){
                normalized = true;
            }

            float exponent = 0.0f;
            stream << element->FirstChildElement("Exponent")->GetText() << std::endl;
            stream >> exponent;

            // Create and store the BRDF in scene.
            this->brdfs.push_back(new BrdfModifiedBlinnPhong(exponent, normalized, id));

            element = element->NextSiblingElement("ModifiedBlinnPhong");
        }

        element = root->FirstChildElement("BRDFs");
        element = element->FirstChildElement("OriginalBlinnPhong");
        while(element)
        {
            int id = -1;
            stream << element->Attribute("id") << std::endl;
            stream >> id;

            float exponent = 0.0f;
            stream << element->FirstChildElement("Exponent")->GetText() << std::endl;
            stream >> exponent;

            // Create and store the BRDF in scene.
            this->brdfs.push_back(new BrdfBlinnPhong(exponent, false, id));

            element = element->NextSiblingElement("OriginalBlinnPhong");
        }

        element = root->FirstChildElement("BRDFs");
        element = element->FirstChildElement("OriginalPhong");
        while(element)
        {
            int id = -1;
            stream << element->Attribute("id") << std::endl;
            stream >> id;
            
            float exponent = 0.0f;
            stream << element->FirstChildElement("Exponent")->GetText() << std::endl;
            stream >> exponent;

            // Create and store the BRDF in scene.
            this->brdfs.push_back(new BrdfPhong(exponent, false, id));

            element = element->NextSiblingElement("OriginalPhong");
        }

        element = root->FirstChildElement("BRDFs");
        element = element->FirstChildElement("ModifiedPhong");
        while(element)
        {
            int id = -1;
            stream << element->Attribute("id") << std::endl;
            stream >> id;

            bool normalized = false;
            if(element->Attribute("normalized", "true")){
                normalized = true;
            }

            float exponent = 0.0f;
            stream << element->FirstChildElement("Exponent")->GetText() << std::endl;
            stream >> exponent;

            // Create and store the BRDF in scene.
            this->brdfs.push_back(new BrdfModifiedPhong(exponent, normalized, id));

            element = element->NextSiblingElement("ModifiedPhong");
        } 
        
        element = root->FirstChildElement("BRDFs");
        element = element->FirstChildElement("TorranceSparrow");
        while(element)
        {
            int id = -1;
            stream << element->Attribute("id") << std::endl;
            stream >> id;

            bool kdfresnel = false;
            if(element->Attribute("kdfresnel", "true")){
                kdfresnel = true;
            }

            float exponent = 0.0f;
            stream << element->FirstChildElement("Exponent")->GetText() << std::endl;
            stream >> exponent;

            // Create and store the BRDF in scene.
            this->brdfs.push_back(new BrdfTorranceSparrow(exponent, kdfresnel, id));

            element = element->NextSiblingElement("TorranceSparrow");
        }
    }
}

void DorkTracer::Scene::parseLights(tinyxml2::XMLNode* root)
{
    auto lights = root->FirstChildElement("Lights");
    if(lights == NULL) return;

    std::stringstream stream;

    auto child = lights->FirstChildElement("AmbientLight");
    if(child){
        stream << child->GetText() << std::endl;
        stream >> ambient_light.x >> ambient_light.y >> ambient_light.z;
    }
    else{
        ambient_light = {0,0,0};
    }
    
    auto light = lights->FirstChildElement("PointLight");
    while (light)
    {
        int id;
        stream << light->Attribute("id") << std::endl;
        stream >> id;

        child = light->FirstChildElement("Position");
        stream << child->GetText() << std::endl;
        child = light->FirstChildElement("Intensity");
        stream << child->GetText() << std::endl;

        Vec3f pos, intensity;
        stream >> pos.x >> pos.y >> pos.z;
        stream >> intensity.x >> intensity.y >> intensity.z;

        point_lights.push_back(PointLight(id , pos, intensity));
        light = light->NextSiblingElement("PointLight");
    }

    light = lights->FirstChildElement("AreaLight");
    while (light)
    {
        Vec3f pos,normal, radiance;
        float size;
        int id;
        
        stream << light->Attribute("id") << std::endl;
        stream >> id;
        
        child = light->FirstChildElement("Position");
        stream << child->GetText() << std::endl;
        stream >> pos.x >> pos.y >> pos.z;

        child = light->FirstChildElement("Normal");
        stream << child->GetText() << std::endl;
        stream >> normal.x >> normal.y >> normal.z;
        
        child = light->FirstChildElement("Radiance");
        stream << child->GetText() << std::endl;
        stream >> radiance.x >> radiance.y >> radiance.z;

        child = light->FirstChildElement("Size");
        stream << child->GetText() << std::endl;
        stream >> size;

        areaLights.push_back(new AreaLight(id, pos, normal, radiance, size));

        light = light->NextSiblingElement("AreaLight");
    }

    light = lights->FirstChildElement("DirectionalLight");
    while (light)
    {
        Vec3f radiance, dir;
        int id;
        
        stream << light->Attribute("id") << std::endl;
        stream >> id;
        

        child = light->FirstChildElement("Direction");
        stream << child->GetText() << std::endl;
        stream >> dir.x >> dir.y >> dir.z;
        
        child = light->FirstChildElement("Radiance");
        stream << child->GetText() << std::endl;
        stream >> radiance.x >> radiance.y >> radiance.z;

        directionalLights.push_back(new DirectionalLight(id, dir, radiance));

        light = light->NextSiblingElement("DirectionalLight");
    }

    light = lights->FirstChildElement("SpotLight");
    while (light)
    {
        Vec3f pos, dir, intensity;
        int id;
        float coverageAngle, falloffAngle;
        
        stream << light->Attribute("id") << std::endl;
        stream >> id;
        
        child = light->FirstChildElement("Position");
        stream << child->GetText() << std::endl;
        stream >> pos.x >> pos.y >> pos.z;

        child = light->FirstChildElement("Direction");
        stream << child->GetText() << std::endl;
        stream >> dir.x >> dir.y >> dir.z;
        
        child = light->FirstChildElement("Intensity");
        stream << child->GetText() << std::endl;
        stream >> intensity.x >> intensity.y >> intensity.z;

        child = light->FirstChildElement("CoverageAngle");
        stream << child->GetText() << std::endl;
        stream >> coverageAngle;

        child = light->FirstChildElement("FalloffAngle");
        stream << child->GetText() << std::endl;
        stream >> falloffAngle;

        spotLights.push_back(new SpotLight(id, pos, dir, intensity, coverageAngle, falloffAngle));
        light = light->NextSiblingElement("SpotLight");
    }
}

void DorkTracer::Scene::parseMaterials(tinyxml2::XMLNode* root)
{
    // Get Materials
    tinyxml2::XMLElement* element = root->FirstChildElement("Materials");

    element = element->FirstChildElement("Material");
    DorkTracer::Material material;
    std::stringstream stream;

    while (element)
    {
        int matId = -1;
        stream << element->Attribute("id") << std::endl;
        stream >> material.id;

        if(element->Attribute("BRDF")){
            int brdfID = -1;
            stream << element->Attribute("BRDF") << std::endl;
            stream >> brdfID;
            
            // find & assign the correct BRDF. Assumes BRDF parsing is already DONE.
            BRDF* brdf = nullptr;
            for(int i = 0; i < brdfs.size(); i++){
                if(brdfs[i]->id == brdfID){
                    brdf = brdfs[i];
                    break;
                }
            }
            material.brdf = brdf;
        }

        if(element->Attribute("type", "mirror"))
        {
            material.type = DorkTracer::Material::Mirror;
        }
        else if(element->Attribute("type", "dielectric")){
            material.type = DorkTracer::Material::Dielectric;
        }
        else if(element->Attribute("type", "conductor")){
            material.type = DorkTracer::Material::Conductor;
        }
        else{
            material.type = DorkTracer::Material::Default;
        }

        bool degamma = false;
        float gamma = 2.2f;
        if(element->Attribute("degamma", "true"))
        {
            degamma = true;
        }

        auto child = element->FirstChildElement("AmbientReflectance");
        if(child != NULL)
        {
            stream << child->GetText() << std::endl;
            stream >> material.ambient.x >> material.ambient.y >> material.ambient.z;
            if(degamma)
            {
                material.ambient.x = std::pow(material.ambient.x, gamma);
                material.ambient.y = std::pow(material.ambient.y, gamma);
                material.ambient.z = std::pow(material.ambient.z, gamma);
            }
        }

        child = element->FirstChildElement("DiffuseReflectance");
        if( child != NULL){
            stream << child->GetText() << std::endl;
            stream >> material.diffuse.x >> material.diffuse.y >> material.diffuse.z;
    
            if(degamma)
            {
                material.diffuse.x = std::pow(material.diffuse.x, gamma);
                material.diffuse.y = std::pow(material.diffuse.y, gamma);
                material.diffuse.z = std::pow(material.diffuse.z, gamma);
                std::cout <<"diffuse after degamma: " << material.diffuse.x << std::endl;
            }
        }

        child = element->FirstChildElement("SpecularReflectance");
        if( child != NULL){
            stream << child->GetText() << std::endl;
            stream >> material.specular.x >> material.specular.y >> material.specular.z;
      
            if(degamma)
            {
                material.specular.x = std::pow(material.specular.x, gamma);
                material.specular.y = std::pow(material.specular.y, gamma);
                material.specular.z = std::pow(material.specular.z, gamma);
            }
        }

        child = element->FirstChildElement("MirrorReflectance");
        if( child != NULL){
            stream << child->GetText() << std::endl;
            stream >> material.mirror.x >> material.mirror.y >> material.mirror.z;
     
            if(degamma)
            {
                material.mirror.x = std::pow(material.mirror.x, gamma);
                material.mirror.y = std::pow(material.mirror.y, gamma);
                material.mirror.z = std::pow(material.mirror.z, gamma);
            }
        }
        else
        {
            material.mirror.x = material.mirror.y = material.mirror.z = 0.0f;
        }

        child = element->FirstChildElement("RefractionIndex");
        if( child != NULL){
            stream << child->GetText() << std::endl;
            stream >> material.refractiveIndex;
            // if(degamma)
            // {
            //     material.refractiveIndex= std::pow(material.refractiveIndex, gamma);
            // }
        }else
        {
            material.refractiveIndex = 1.0f;
        }

        child = element->FirstChildElement("AbsorptionCoefficient");
        if( child != NULL){
            stream << child->GetText() << std::endl;
            stream >> material.absorptionCoefficient.x >>  material.absorptionCoefficient.y >> material.absorptionCoefficient.z;
            // if(degamma)
            // {
            //     material.absorptionCoefficient.x = std::pow(material.absorptionCoefficient.x, gamma);
            //     material.absorptionCoefficient.y = std::pow(material.absorptionCoefficient.y, gamma);
            //     material.absorptionCoefficient.z = std::pow(material.absorptionCoefficient.z, gamma);
            // }
        }else
        {
            material.absorptionCoefficient = Vec3f{0,0,0};
        }

        child = element->FirstChildElement("AbsorptionIndex");
        if( child != NULL){
            stream << child->GetText() << std::endl;
            stream >> material.conductorAbsorptionIndex;
        }else
        {
            material.conductorAbsorptionIndex = 0.0f;
        }
        
        child = element->FirstChildElement("PhongExponent");
        if( child != NULL)
        {
            stream << child->GetText() << std::endl;
            stream >> material.phong_exponent;
        }
        else{
            material.phong_exponent = 1.0f;
        }

        child = element->FirstChildElement("Roughness");
        if( child != NULL)
        {
            stream << child->GetText() << std::endl;
            stream >> material.roughness;
        }
        else{
            material.roughness = 0.0f;
        }

        materials.push_back(material);
        element = element->NextSiblingElement("Material");
    }
}

void DorkTracer::Scene::parseMeshes(tinyxml2::XMLNode* root, std::string elemName)
{ 
    auto element = root->FirstChildElement("Objects");
    element = element->FirstChildElement(elemName.c_str());
    std::stringstream stream;

    while (element)
    {
        auto child = element->FirstChildElement("Faces");
        bool hasPlyData = child->Attribute("plyFile") != NULL;

        DorkTracer::Mesh* mesh;

        if(hasPlyData)
        {
            std::vector<Vec3f> meshVertices;
            std::vector<Vec2f> meshUV;
            if(elemName == "Mesh")
            {
                mesh = new DorkTracer::Mesh(meshVertices, meshUV);
            }
            else if(elemName == "LightMesh")
            {
                Vec3f radiance{0,0,0};
                child = element->FirstChildElement("Radiance");
                if(child){
                    stream << child->GetText() << std::endl;
                    stream >> radiance.x >> radiance.y >> radiance.z;
                }
                mesh = (Mesh*)new DorkTracer::MeshLight(meshVertices, meshUV, radiance);
            }
        }
        else{
            if(elemName == "Mesh")
            {
                mesh = new DorkTracer::Mesh(this->vertex_data, this->texCoords);
            }
            else if(elemName == "LightMesh")
            { 
                Vec3f radiance{0,0,0};
                child = element->FirstChildElement("Radiance");
                if(child){
                    stream << child->GetText() << std::endl;
                    stream >> radiance.x >> radiance.y >> radiance.z;
                }
                mesh = (Mesh*)new DorkTracer::MeshLight(this->vertex_data, this->texCoords, radiance);
            }
        }
        mesh->isInstance = false;

        stream << element->Attribute("id") << std::endl;
        stream >> mesh->id;

        // Setup textures
        child = element->FirstChildElement("Textures");
        if(child){
            std::string texturesInp;
            texturesInp = child->GetText();
            texturesInp += " ";
            SetupTextures(mesh, texturesInp);
        }
        // Compute Model, invModel and invTransModel matrices from given transformations.
        child = element->FirstChildElement("Transformations");

        mesh->transform.MakeIdentity();
        mesh->inverseTransform.MakeIdentity();
        mesh->inverseTransposeTransform.MakeIdentity();

        if(child != nullptr)
        {
            std::string input = child->GetText();
            computeTransform(mesh, input);
        }

        // Common ops to both cases.
        child = element->FirstChildElement("Material");
        stream << child->GetText() << std::endl;
        int matId = 0;
        stream >> matId;
        mesh->SetMaterial(matId);

        // parse motion blur
        // <MotionBlur>0 0 4</MotionBlur>
        child = element->FirstChildElement("MotionBlur");
        if(child != nullptr)
        {
            stream << child->GetText() << std::endl;
            stream >> mesh->motionBlurVector.x >> mesh->motionBlurVector.y >> mesh->motionBlurVector.z;

            mesh->hasMotionBlur = true;
            // enable motior blur globally even if only one object is using it.
            this->isMotionBlurEnabled = true; 
        }
        else{
            mesh->hasMotionBlur = false;
            mesh->motionBlurVector = Vec3f(0.0f,0.0f,0.0f);
        }

        child = element->FirstChildElement("Faces");

        int vertexOffset = 0;
        if(child->Attribute("vertexOffset") != NULL){
            stream << child->Attribute("vertexOffset") << std::endl;
            stream >> vertexOffset;
        }
        int textureOffset = 0;
        if(child->Attribute("textureOffset") != NULL){
            stream << child->Attribute("textureOffset") << std::endl;
            stream >> textureOffset;
        }
        mesh->SetAccessOffsets(vertexOffset, textureOffset);

        BoundingBox bbox;
        bbox.maxCorner.x = bbox.maxCorner.y = bbox.maxCorner.z = std::numeric_limits<float>::min();
        bbox.minCorner.x = bbox.minCorner.y = bbox.minCorner.z = std::numeric_limits<float>::max();
        
        if(hasPlyData)
        {
            // Face data is given as plyFile
            std::string filename;
            stream << child->Attribute("plyFile");
            stream >> filename;
            // std::cout <<"Ply Filename:"<<filename<<std::endl;
            // Construct the data object by reading from file
            happly::PLYData plyIn(child->Attribute("plyFile"));

            std::vector<std::array<double, 3>> vPos = plyIn.getVertexPositions();
            for(int i = 0; i < vPos.size(); i++){
                Vec3f vert { vPos[i][0], vPos[i][1], vPos[i][2]};
                mesh->AddVertex(vert);
            }
            vPos.clear();

            // for(auto str : plyIn.getElementNames()){
            //     std::cout << str << std::endl;
            // }
            std::vector<std::vector<int>> fInd = plyIn.getFaceIndices<int>();
            // std::cout<<"has: "<<  fInd.size() << " Triangles."<< std::endl;
            for(int i = 0; i < fInd.size(); i++)
            {
                if(fInd[i].size() == 3){
                    // ply format is 0 base indexed, whereas rest of our project uses 1-based indexing, so add that offset.
                    DorkTracer::Face face = createFace(fInd[i][0]+1, fInd[i][1] + 1, fInd[i][2] + 1, mesh);

                    // calculate bbox of the owner mesh here
                    updateBBox(&bbox, face);
                    mesh->faces.push_back(face);
                }
                else if(fInd[i].size() == 4)
                {
                    // ply format is 0 base indexed, whereas rest of our project uses 1-based indexing, so add that offset.
                    DorkTracer::Face face1 = createFace(fInd[i][0]+1, fInd[i][1] + 1, fInd[i][2] + 1, mesh);
                    updateBBox(&bbox, face1);
                    mesh->faces.push_back(face1);

                    DorkTracer::Face face2 = createFace(fInd[i][2]+1, fInd[i][3] + 1, fInd[i][0] + 1, mesh);
                    updateBBox(&bbox, face2);
                    mesh->faces.push_back(face2);

                }
                else{
                    std::cout<<"A face is assumed to have 3 or 4 indices. can not parse. index count:" << fInd[i].size() << std::endl;
                }
            }
        }
        else
        {
            stream << child->GetText() << std::endl;
            DorkTracer::Face face;
            while (!(stream >> face.v0_id).eof())
            {
                stream >> face.v1_id >> face.v2_id;
                computeFaceProperties(face, mesh);
                
                // calculate bbox of the owner mesh here
                bbox.minCorner.x = std::min(face.bbox.minCorner.x,  bbox.minCorner.x);
                bbox.minCorner.y = std::min(face.bbox.minCorner.y,  bbox.minCorner.y);
                bbox.minCorner.z = std::min(face.bbox.minCorner.z,  bbox.minCorner.z);

                bbox.maxCorner.x = std::max(face.bbox.maxCorner.x, bbox.maxCorner.x);
                bbox.maxCorner.y = std::max(face.bbox.maxCorner.y, bbox.maxCorner.y);
                bbox.maxCorner.z = std::max(face.bbox.maxCorner.z, bbox.maxCorner.z);

                mesh->faces.push_back(face);
            }
        }
        stream.clear();

        mesh->bbox = bbox;
        mesh->ConstructBVH();

        if(elemName == "Mesh")
        {
            meshes.push_back(mesh);
        }
        else if(elemName == "LightMesh")
        {
            MeshLight* meshLight = (MeshLight*) mesh;
            std::cout<<"meshlight has: "  << meshLight->faces.size() <<" many faces."<<std::endl;
            meshLight->setupFaceSelectionDistro(meshLight->faces.size());

            meshLights.push_back(meshLight);
            meshes.push_back(mesh);

            // mark its material as emissive/light?
            Material& mat = this->materials[mesh->GetMaterial()-1];
            mat.type = DorkTracer::Material::Emissive;
            mat.radiance = meshLight->radiance;
            
            std::cout<<"parsed a light mesh with radiance:" << mat.radiance.x << std::endl;
        }
        // mesh.faces.clear();
        // std::cout << " total area of mesh: " << mesh->surfaceArea << std::endl;
        element = element->NextSiblingElement(elemName.c_str());
    }
    stream.clear();
}

void DorkTracer::Scene::parseCameras(tinyxml2::XMLNode* root)
{
    auto element = root->FirstChildElement("Cameras");
    element = element->FirstChildElement("Camera");

    std::stringstream stream;
    DorkTracer::Camera camera;
    while (element)
    {
        bool isLookAt = element->Attribute("type", "lookAt") != NULL;
        Vec3f camPos, upDir;
        float nearDist, width, height;
        std::string imageName;

        auto child = element->FirstChildElement("Position");
        stream << child->GetText() << std::endl;
        stream >> camPos.x >> camPos.y >> camPos.z;

        child = element->FirstChildElement("Up");
        stream << child->GetText() << std::endl;
        stream >> upDir.x >> upDir.y >> upDir.z;

        child = element->FirstChildElement("NearDistance");
        stream << child->GetText() << std::endl;
        stream >> nearDist;

        child = element->FirstChildElement("ImageResolution");
        stream << child->GetText() << std::endl;
        stream >> width >> height;

        child = element->FirstChildElement("ImageName");
        stream << child->GetText() << std::endl;
        stream >> imageName;

        if(isLookAt){

            Vec3f gazePoint;
            float fovY;

            child = element->FirstChildElement("GazePoint");
            if(child == NULL){
                child = element->FirstChildElement("Gaze");
            }
            stream << child->GetText() << std::endl;
            stream >> gazePoint.x >> gazePoint.y >> gazePoint.z;

            child = element->FirstChildElement("FovY");
            stream << child->GetText() << std::endl;
            stream >> fovY;
            camera.SetupLookAt(camPos, gazePoint, upDir, nearDist, fovY, width, height, imageName);
        }
        else{
            Vec3f gazeDir;
            Vec4f nearPlane;

            child = element->FirstChildElement("Gaze");
            stream << child->GetText() << std::endl;
            stream >> gazeDir.x >> gazeDir.y >> gazeDir.z;

            child = element->FirstChildElement("NearPlane");
            stream << child->GetText() << std::endl;
            stream >> nearPlane.x >> nearPlane.y >> nearPlane.z >> nearPlane.w;
            camera.SetupDefault(camPos, gazeDir, upDir, nearPlane, nearDist, width, height, imageName);
        }
        
        int numSamples = 1;
        child = element->FirstChildElement("NumSamples");
        if(child != nullptr){
            stream << child->GetText() << std::endl;
            stream >> numSamples;
        }
        camera.samplesPerPixel = numSamples;
        std::cout << "Num samples: " <<camera.samplesPerPixel << std::endl;


        float focusDistance = 0.0f;
        child = element->FirstChildElement("FocusDistance");
        if(child != nullptr){
            stream << child->GetText() << std::endl;
            stream >> focusDistance;
        }
        camera.focusDistance = focusDistance;

        float apertureSize = 0.0f;
        child = element->FirstChildElement("ApertureSize");
        if(child != nullptr){
            stream << child->GetText() << std::endl;
            stream >> apertureSize;
        }
        camera.apertureSize = apertureSize;

        // Parse Render parameters, mainly for Path Tracing.
        child = element->FirstChildElement("Renderer");
        if(child != nullptr)
        {
            std::string rendererType;
            stream << child->GetText() << std::endl;
            stream >> rendererType;
            if(rendererType == "PathTracing")
            {
                child = element->FirstChildElement("RendererParams");
                bool importanceSampling, russianRoulette,nextEventEstimation;
                importanceSampling = russianRoulette = nextEventEstimation = false;

                if(child != nullptr)
                {
                    stream << child->GetText() << std::endl;
                    std::string param;
                    while(stream >> param) {
                        std::cout << param << std::endl;
                        if(param == "NextEventEstimation"){
                            nextEventEstimation= true;
                        }
                        else if(param == "RussianRoulette"){
                            russianRoulette = true;
                        }
                        else if(param == "ImportanceSampling"){
                            importanceSampling = true;
                        }
                    }
                }
                stream.clear();
                
                RendererParams rendererParams;
                rendererParams.SetParams(true, importanceSampling, nextEventEstimation, russianRoulette);
                camera.SetRendererParams(rendererParams);
            }
            else{
                std::cout << "Unknown renderer type: " << rendererType << std::endl;
            }
        }


        parseTonemapper(element, &camera);

        cameras.push_back(camera);
        element = element->NextSiblingElement("Camera");
    }
}