#include "parser.h"
#include "tinyxml2.h"
#include <sstream>
#include <stdexcept>
#include "helperMath.h"

void parser::Scene::loadFromXml(const std::string &filepath)
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
    Camera camera;
    while (element)
    {
        camera.isLookAt = element->Attribute("type", "lookAt") != NULL;

        auto child = element->FirstChildElement("Position");
        stream << child->GetText() << std::endl;
        stream >> camera.position.x >> camera.position.y >> camera.position.z;

        if(camera.isLookAt){
            child = element->FirstChildElement("GazePoint");
            stream << child->GetText() << std::endl;
            stream >> camera.gazePoint.x >> camera.gazePoint.y >> camera.gazePoint.z;

            child = element->FirstChildElement("FovY");
            stream << child->GetText() << std::endl;
            stream >> camera.fovY;
            camera.near_plane = Vec4f{0,0,0,0};
            std::cout << "LookAt camera. Using fovY and GazePOint insted of NearPlane&Gaze dir to configure camera. fovY: "<< camera.fovY << std::endl;
        }
        else{
            child = element->FirstChildElement("Gaze");
            stream << child->GetText() << std::endl;
            stream >> camera.gaze.x >> camera.gaze.y >> camera.gaze.z;

            child = element->FirstChildElement("NearPlane");
            stream << child->GetText() << std::endl;
            stream >> camera.near_plane.x >> camera.near_plane.y >> camera.near_plane.z >> camera.near_plane.w;
            camera.fovY = 0.0f;
        }

        child = element->FirstChildElement("Up");
        stream << child->GetText() << std::endl;
        stream >> camera.up.x >> camera.up.y >> camera.up.z;

        child = element->FirstChildElement("NearDistance");
        stream << child->GetText() << std::endl;
        child = element->FirstChildElement("ImageResolution");
        stream << child->GetText() << std::endl;
        child = element->FirstChildElement("ImageName");
        stream << child->GetText() << std::endl;

        stream >> camera.near_distance;
        stream >> camera.image_width >> camera.image_height;
        stream >> camera.image_name;

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

    //Get Materials
    element = root->FirstChildElement("Materials");
    element = element->FirstChildElement("Material");
    Material material;
    while (element)
    {
        if((element->Attribute("type", "mirror") != NULL))
        {
            material.type = Mirror;
        }
        else if((element->Attribute("type", "dielectric") != NULL)){
            material.type = Dielectric;
        }
        else if((element->Attribute("type", "conductor") != NULL)){
            material.type = Conductor;
        }
        else{
            material.type = Default;
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
            // cout <<"has mirror reflectance??"<<endl;
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
            std::cout <<"has AbsorptionCoefficent"<<std::endl;
            stream << child->GetText() << std::endl;
            stream >> material.absorptionCoefficient.x >>  material.absorptionCoefficient.y >> material.absorptionCoefficient.z;
        }else
        {
            material.absorptionCoefficient = Vec3f{0,0,0};
        }

        child = element->FirstChildElement("AbsorptionIndex");
        if( child != NULL){
            std::cout <<"has Absorption Index, so must be a conductor."<<std::endl;
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

    //Get VertexData
    element = root->FirstChildElement("VertexData");
    stream << element->GetText() << std::endl;
    Vec3f vertex;
    while (!(stream >> vertex.x).eof())
    {
        stream >> vertex.y >> vertex.z;
        vertex_data.push_back(vertex);
    }
    stream.clear();

    //Get Meshes
    element = root->FirstChildElement("Objects");
    element = element->FirstChildElement("Mesh");
    Mesh mesh;
    while (element)
    {
        child = element->FirstChildElement("Material");
        stream << child->GetText() << std::endl;
        stream >> mesh.material_id;

        child = element->FirstChildElement("Faces");

        
        if((child->Attribute("plyFile") != NULL))
        {
         // <Faces plyFile="ply/dragon_remeshed_fixed.ply" />
            // Face data is given as plyFile
            std::cout <<"ply file!!"<<std::endl;
            
            std::string filename;
            stream << child->Attribute("plyFile");
            stream >> filename;
            std::cout <<"Filename:"<<filename<<std::endl;
            // Construct the data object by reading from file
            happly::PLYData plyIn(child->Attribute("plyFile"));

            // Get mesh-style data from the object
            mesh.useOwnVertices = true;

            std::vector<std::array<double, 3>> vPos = plyIn.getVertexPositions();
            for(int i = 0; i < vPos.size(); i++){
                Vec3f vert { vPos[i][0], vPos[i][1], vPos[i][2]};
                mesh.vertices.push_back(vert);
            }
            vPos.clear();

            std::vector<std::vector<size_t>> fInd = plyIn.getFaceIndices<size_t>();
            std::cout<<"has: "<<  fInd.size() << " faces."<< std::endl;
            for(int i = 0; i < fInd.size(); i++){

                Face face;
                if(fInd[i].size() != 3){
                    std::cout<<"A face is assumed to have 3 indices. can not parse. index count:" << fInd[i].size() << std::endl;
                    return;
                }
                // ply format is 0 base indexed, whereas rest of our project uses 1-based indexing, so add that offset.
                face.v0_id = fInd[i][0] + 1;
                face.v1_id = fInd[i][1] + 1;
                face.v2_id = fInd[i][2] + 1;
                computeFaceNormal(face, mesh.vertices);
                mesh.faces.push_back(face);
            }
        }
        else
        {
            // regular parsing
            mesh.useOwnVertices = false;

            stream << child->GetText() << std::endl;
            Face face;
            while (!(stream >> face.v0_id).eof())
            {
                stream >> face.v1_id >> face.v2_id;
                computeFaceNormal(face, this->vertex_data);

                mesh.faces.push_back(face);
            }
        }

        stream.clear();

        meshes.push_back(mesh);
        mesh.faces.clear();
        element = element->NextSiblingElement("Mesh");
    }
    stream.clear();

    //Get Triangles
    element = root->FirstChildElement("Objects");
    element = element->FirstChildElement("Triangle");
    Triangle triangle;
    while (element)
    {
        child = element->FirstChildElement("Material");
        stream << child->GetText() << std::endl;
        stream >> triangle.material_id;

        child = element->FirstChildElement("Indices");
        stream << child->GetText() << std::endl;
        stream >> triangle.indices.v0_id >> triangle.indices.v1_id >> triangle.indices.v2_id;

        computeFaceNormal(triangle.indices, this->vertex_data);

        triangles.push_back(triangle);
        element = element->NextSiblingElement("Triangle");
    }

    //Get Spheres
    element = root->FirstChildElement("Objects");
    element = element->FirstChildElement("Sphere");
    Sphere sphere;
    while (element)
    {
        child = element->FirstChildElement("Material");
        stream << child->GetText() << std::endl;
        stream >> sphere.material_id;

        child = element->FirstChildElement("Center");
        stream << child->GetText() << std::endl;
        stream >> sphere.center_vertex_id;

        child = element->FirstChildElement("Radius");
        stream << child->GetText() << std::endl;
        stream >> sphere.radius;

        spheres.push_back(sphere);
        element = element->NextSiblingElement("Sphere");
    }
}

void Scene::computeFaceNormal(Face& face, std::vector<Vec3f>& vertices)
{
    Vec3f a,b,c,ab,ac;
    a = vertices[face.v0_id -1];
    b = vertices[face.v1_id -1];
    c = vertices[face.v2_id -1];
    ab = b - a;
    ac = c - a;
    face.n = makeUnit(cross(ab, ac));
}
