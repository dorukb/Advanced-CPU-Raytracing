#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "happly.h"
#include "scene.h"
#include <limits>
#include "matrix.hpp"
#include "instancedMesh.hpp"

void DorkTracer::Scene::loadFromXml(const std::string &filepath)
{
    tinyxml2::XMLDocument file;
    std::stringstream stream;

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

    //Get BackgroundColor
    auto element = root->FirstChildElement("BackgroundColor");
    if (element)
    {
        stream << element->GetText() << std::endl;
    }
    else
    {
        stream << "0 0 0" << std::endl;
    }
    stream >> background_color.x >> background_color.y >> background_color.z;

    //Get ShadowRayEpsilon
    element = root->FirstChildElement("ShadowRayEpsilon");
    if (element)
    {
        stream << element->GetText() << std::endl;
    }
    else
    {
        stream << "0.001" << std::endl;
    }
    stream >> shadow_ray_epsilon;

    //Get MaxRecursionDepth
    element = root->FirstChildElement("MaxRecursionDepth");
    if (element)
    {
        stream << element->GetText() << std::endl;
    }
    else
    {
        stream << "0" << std::endl;
    }
    stream >> max_recursion_depth;

    //Get Cameras
    element = root->FirstChildElement("Cameras");
    element = element->FirstChildElement("Camera");
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
        cameras.push_back(camera);
        element = element->NextSiblingElement("Camera");
    }

    //Get Lights
    element = root->FirstChildElement("Lights");
    auto child = element->FirstChildElement("AmbientLight");
    stream << child->GetText() << std::endl;
    stream >> ambient_light.x >> ambient_light.y >> ambient_light.z;
    element = element->FirstChildElement("PointLight");
    PointLight point_light;
    while (element)
    {
        child = element->FirstChildElement("Position");
        stream << child->GetText() << std::endl;
        child = element->FirstChildElement("Intensity");
        stream << child->GetText() << std::endl;

        stream >> point_light.position.x >> point_light.position.y >> point_light.position.z;
        stream >> point_light.intensity.x >> point_light.intensity.y >> point_light.intensity.z;

        point_lights.push_back(point_light);
        element = element->NextSiblingElement("PointLight");
    }

    // Get Materials
    element = root->FirstChildElement("Materials");
    element = element->FirstChildElement("Material");
    DorkTracer::Material material;
    while (element)
    {
        if((element->Attribute("type", "mirror") != NULL))
        {
            material.type = DorkTracer::Material::Mirror;
        }
        else if((element->Attribute("type", "dielectric") != NULL)){
            material.type = DorkTracer::Material::Dielectric;
        }
        else if((element->Attribute("type", "conductor") != NULL)){
            material.type = DorkTracer::Material::Conductor;
        }
        else{
            material.type = DorkTracer::Material::Default;
        }

        child = element->FirstChildElement("AmbientReflectance");
        if( child != NULL){
            stream << child->GetText() << std::endl;
            stream >> material.ambient.x >> material.ambient.y >> material.ambient.z;
        }

        child = element->FirstChildElement("DiffuseReflectance");
        if( child != NULL){
            stream << child->GetText() << std::endl;
            stream >> material.diffuse.x >> material.diffuse.y >> material.diffuse.z;
        }

        child = element->FirstChildElement("SpecularReflectance");
        if( child != NULL){
            stream << child->GetText() << std::endl;
            stream >> material.specular.x >> material.specular.y >> material.specular.z;
        }

        child = element->FirstChildElement("MirrorReflectance");
        if( child != NULL){
            stream << child->GetText() << std::endl;
            stream >> material.mirror.x >> material.mirror.y >> material.mirror.z;
        }else
        {
            material.mirror.x = material.mirror.y = material.mirror.z = 0.0f;
        }

        child = element->FirstChildElement("RefractionIndex");
        if( child != NULL){
            stream << child->GetText() << std::endl;
            stream >> material.refractiveIndex;
        }else
        {
            material.refractiveIndex = 1.0f;
        }

        child = element->FirstChildElement("AbsorptionCoefficient");
        if( child != NULL){
            stream << child->GetText() << std::endl;
            stream >> material.absorptionCoefficient.x >>  material.absorptionCoefficient.y >> material.absorptionCoefficient.z;
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


        materials.push_back(material);
        element = element->NextSiblingElement("Material");
    }

    // Get VertexData
    element = root->FirstChildElement("VertexData");
    stream << element->GetText() << std::endl;
    Vec3f vertex;
    while (!(stream >> vertex.x).eof())
    {
        stream >> vertex.y >> vertex.z;
        vertex_data.push_back(vertex);
    }
    stream.clear();

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
    element = root->FirstChildElement("Objects");
    element = element->FirstChildElement("Mesh");
    while (element)
    {
        child = element->FirstChildElement("Faces");
        bool hasPlyData = child->Attribute("plyFile") != NULL;

        DorkTracer::Mesh* mesh;

        if(hasPlyData)
        {
            std::vector<Vec3f> meshVertices;
            mesh = new DorkTracer::Mesh(meshVertices);
        }
        else{
            mesh = new DorkTracer::Mesh(this->vertex_data);
        }

        stream << element->Attribute("id") << std::endl;
        stream >> mesh->id;

        // Compute Model, invModel and invTransModel matrices from given transformations.
        child = element->FirstChildElement("Transformations");
        mesh->transform.MakeIdentity();
        if(child != nullptr)
        {
            computeTransform(mesh, child);
        }
        else{
            mesh->inverseTransform.MakeIdentity();
            mesh->inverseTransposeTransform.MakeIdentity();
        }

        // Common ops to both cases.
        child = element->FirstChildElement("Material");
        stream << child->GetText() << std::endl;
        int matId = 0;
        stream >> matId;
        mesh->SetMaterial(matId);

        child = element->FirstChildElement("Faces");

        BoundingBox bbox;
        bbox.maxCorner.x = bbox.maxCorner.y = bbox.maxCorner.z = std::numeric_limits<float>::min();
        bbox.minCorner.x = bbox.minCorner.y = bbox.minCorner.z = std::numeric_limits<float>::max();
        
        if(hasPlyData)
        {
         // <Faces plyFile="ply/dragon_remeshed_fixed.ply" />
            // Face data is given as plyFile
            
            std::string filename;
            stream << child->Attribute("plyFile");
            stream >> filename;
            std::cout <<"Ply Filename:"<<filename<<std::endl;
            // Construct the data object by reading from file
            happly::PLYData plyIn(child->Attribute("plyFile"));

            std::vector<std::array<double, 3>> vPos = plyIn.getVertexPositions();
            for(int i = 0; i < vPos.size(); i++){
                Vec3f vert { vPos[i][0], vPos[i][1], vPos[i][2]};
                mesh->vertices.push_back(vert);
            }
            vPos.clear();

            std::vector<std::vector<int>> fInd = plyIn.getFaceIndices<int>();
            std::cout<<"has: "<<  fInd.size() << " faces and "<< mesh->vertices.size() << " vertices."<< std::endl;
            for(int i = 0; i < fInd.size(); i++){

                DorkTracer::Face face;
                if(fInd[i].size() != 3){
                    std::cout<<"A face is assumed to have 3 indices. can not parse. index count:" << fInd[i].size() << std::endl;
                    return;
                }
                // ply format is 0 base indexed, whereas rest of our project uses 1-based indexing, so add that offset.
                face.v0_id = fInd[i][0] + 1;
                face.v1_id = fInd[i][1] + 1;
                face.v2_id = fInd[i][2] + 1;
                computeFaceProperties(face, mesh->vertices);
                // calculate bbox of the owner mesh here, ugly but efficient.
                if(face.bbox.minCorner.x < bbox.minCorner.x){
                    bbox.minCorner.x = face.bbox.minCorner.x;
                }
                if(face.bbox.minCorner.y < bbox.minCorner.y){
                    bbox.minCorner.y = face.bbox.minCorner.y;
                }
                if(face.bbox.minCorner.z < bbox.minCorner.z){
                    bbox.minCorner.z = face.bbox.minCorner.z;
                }

                if(face.bbox.maxCorner.x > bbox.maxCorner.x){
                    bbox.maxCorner.x = face.bbox.maxCorner.x;
                }
                if(face.bbox.maxCorner.y > bbox.maxCorner.y){
                    bbox.maxCorner.y = face.bbox.maxCorner.y;
                }
                if(face.bbox.maxCorner.z > bbox.maxCorner.z){
                    bbox.maxCorner.z = face.bbox.maxCorner.z;
                }
                mesh->faces.push_back(face);
            }
        }
        else
        {
            // regular parsing
            // mesh.useOwnVertices = false;

            stream << child->GetText() << std::endl;
            DorkTracer::Face face;
            while (!(stream >> face.v0_id).eof())
            {
                stream >> face.v1_id >> face.v2_id;
                computeFaceProperties(face, mesh->vertices);
                
                // calculate bbox of the owner mesh here, ugly but efficient.
                if(face.bbox.minCorner.x < bbox.minCorner.x){
                    bbox.minCorner.x = face.bbox.minCorner.x;
                }
                if(face.bbox.minCorner.y < bbox.minCorner.y){
                    bbox.minCorner.y = face.bbox.minCorner.y;
                }
                if(face.bbox.minCorner.z < bbox.minCorner.z){
                    bbox.minCorner.z = face.bbox.minCorner.z;
                }

                if(face.bbox.maxCorner.x > bbox.maxCorner.x){
                    bbox.maxCorner.x = face.bbox.maxCorner.x;
                }
                if(face.bbox.maxCorner.y > bbox.maxCorner.y){
                    bbox.maxCorner.y = face.bbox.maxCorner.y;
                }
                if(face.bbox.maxCorner.z > bbox.maxCorner.z){
                    bbox.maxCorner.z = face.bbox.maxCorner.z;
                }
                mesh->faces.push_back(face);
            }
        }

        stream.clear();

        mesh->bbox = bbox;
        mesh->ConstructBVH();

        meshes.push_back(mesh);
        // mesh.faces.clear();
        element = element->NextSiblingElement("Mesh");
    }
    stream.clear();


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

        DorkTracer::Mesh* baseMesh = nullptr;
        for(int i = 0; i < meshes.size(); i++){
            if(meshes[i]->id == baseMeshId){
                baseMesh = (Mesh*) meshes[i];
            }
        }
        DorkTracer::InstancedMesh* meshInstance = new InstancedMesh(baseMesh);
        meshInstance->id = ownId;

        // Get material
        child = element->FirstChildElement("Material");
        stream << child->GetText() << std::endl;
        int matId = 0;
        stream >> matId;
        meshInstance->SetMaterial(matId);

        // Get transformations
        // Compute Model, invModel and invTransModel matrices from given transformations.
        child = element->FirstChildElement("Transformations");

        if(resetTransform == true)
        {
            meshInstance->transform.MakeIdentity();
        }
        else{
            meshInstance->transform = baseMesh->transform;
        }
        
        if(child != nullptr)
        {
            computeTransform(meshInstance, child);
        }
        else{
            meshInstance->inverseTransform.MakeIdentity();
            meshInstance->inverseTransposeTransform.MakeIdentity();
        }

      
        // if(resetTransform == false)
        // {   
        //     // also need to factor in transform matrices of our baseMesh.
        //     // transform matrices of baseMesh must be computed at this point.
        //     meshInstance->transform = meshInstance->transform * baseMesh->transform;
        //     // order must be reversed here.
        //     meshInstance->inverseTransform = baseMesh->inverseTransform * meshInstance->inverseTransform;
        //     meshInstance->inverseTransposeTransform = Matrix::GetTranspose(meshInstance->inverseTransform);
        // }

        // compute and set transformed bbox.
        BoundingBox bbox = baseMesh->bbox;
        bbox.maxCorner = Matrix::ApplyTransform(meshInstance->transform, Vec4f(bbox.maxCorner, 1.0f));
        bbox.minCorner = Matrix::ApplyTransform(meshInstance->transform, Vec4f(bbox.minCorner, 1.0f));
        meshInstance->bbox = bbox;
        
        meshes.push_back(meshInstance);
        element = element->NextSiblingElement("MeshInstance");
    }


    //Get Triangles, they can be represented with Mesh structure as well.
    element = root->FirstChildElement("Objects");
    element = element->FirstChildElement("Triangle");
    
    // std::shared_ptr<Triangle> p = std::make_shared<Triangle>();
    while (element)
    {
        // DorkTracer::Mesh mesh(this->vertex_data);

        DorkTracer::Mesh* mesh = new DorkTracer::Mesh(this->vertex_data);

        // Compute Model, invModel and invTransModel matrices from given transformations.
        child = element->FirstChildElement("Transformations");
        mesh->transform.MakeIdentity();
        if(child != nullptr)
        {
            computeTransform(mesh, child);
        }
        else{
            mesh->inverseTransform.MakeIdentity();
            mesh->inverseTransposeTransform.MakeIdentity();
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
        computeFaceProperties(face, mesh->vertices);
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
            computeTransform(sphere, child);
        }
        else{
            sphere->inverseTransform.MakeIdentity();
            sphere->inverseTransposeTransform.MakeIdentity();
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

        spheres.push_back(sphere);
        element = element->NextSiblingElement("Sphere");
    }

    // std::cout <<"parsing succesfully completed." << std::endl;
}

void DorkTracer::Scene::computeFaceProperties(DorkTracer::Face& face, std::vector<Vec3f>& vertices)
{
    computeFaceCenter(face, vertices);
    computeFaceNormal(face, vertices);
    computeFaceBoundingBox(face, vertices);
}
void DorkTracer::Scene::computeFaceCenter(DorkTracer::Face& face, std::vector<Vec3f>& vertices)
{
    Vec3f a,b,c,ab,ac;
    a = vertices[face.v0_id -1];
    b = vertices[face.v1_id -1];
    c = vertices[face.v2_id -1];
    
    face.center = (a+b+c)/3;
}

void DorkTracer::Scene::computeTransform(DorkTracer::Shape* mesh, tinyxml2::XMLElement* child)
{
    // sample input: <Transformations>s2 s3 r2 t2</Transformations>
            // s:scale, r:rotation, t:translation, numbers are the indices to look up.
    std::string str = child->GetText();
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

            std::cout<<"r"<<id<<" applied."<<std::endl;
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

            std::cout<<"t"<<id<<" applied."<<std::endl;
            Matrix t = Matrix::GetTranslation(translation.x, translation.y , translation.z);
            Matrix invT = Matrix::GetTranslation(-translation.x, -translation.y , -translation.z);
            
            invTransformations.push_back(invT);
            mesh->transform = t * mesh->transform;
        }
        else if(str[idx] == 's'){
            int id = int(str[idx+1]-'0');
            Vec3f scale = this->scalings[id-1];

            std::cout<<"s"<<id<<" applying."<<std::endl;
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
void DorkTracer::Scene::computeFaceNormal(DorkTracer::Face& face, std::vector<Vec3f>& vertices)
{
    Vec3f a,b,c,ab,ac;
    a = vertices[face.v0_id -1];
    b = vertices[face.v1_id -1];
    c = vertices[face.v2_id -1];
    ab = b - a;
    ac = c - a;
    face.n = makeUnit(cross(ab, ac));
}
void DorkTracer::Scene::computeFaceBoundingBox(DorkTracer::Face& face, std::vector<Vec3f>& vertices)
{
    Vec3f& a = vertices[face.v0_id -1];
    Vec3f& b = vertices[face.v1_id -1];
    Vec3f& c = vertices[face.v2_id -1];

    face.bbox.minCorner.x = std::min(std::min(a.x, b.x), c.x);
    face.bbox.minCorner.y = std::min(std::min(a.y, b.y), c.y);
    face.bbox.minCorner.z = std::min(std::min(a.z, b.z), c.z);

    face.bbox.maxCorner.x = std::max(std::max(a.x, b.x), c.x);    
    face.bbox.maxCorner.y = std::max(std::max(a.y, b.y), c.y);
    face.bbox.maxCorner.z = std::max(std::max(a.z, b.z), c.z);
}
