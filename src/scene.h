#ifndef __DORKTRACER_SCENE__
#define __DORKTRACER_SCENE__

#include <string>
#include <vector>
#include "helperMath.h"
#include "camera.hpp"
#include "material.hpp"
#include "mesh.hpp"
#include "shape.hpp"
#include "sphere.hpp"
#include "tinyxml2.h"
#include "areaLight.h"
#include "directionalLight.h"
#include "spotLight.h"
#include "sphericalEnvironmentLight.h"
#include "meshLight.h"

#include "texture.h"
#include "image.h"
#include "texture.h"
#include "imageTexture.h"
#include "perlinTexture.h"
#include "tonemapper.h"
#

namespace DorkTracer{


    struct PointLight
    {
        Vec3f position;
        Vec3f intensity;
    };

    class Scene
    {
        
    public:
        static float shadow_ray_epsilon;

        Vec3i background_color;
        Texture* bgTexture;

        int max_recursion_depth;
        bool isMotionBlurEnabled;
        Vec3f ambient_light;

        std::vector<Camera> cameras;
        std::vector<PointLight> point_lights;
        std::vector<AreaLight*> areaLights;
        std::vector<DirectionalLight*> directionalLights;
        std::vector<SphericalEnvironmentLight*> sphericalEnvLights;
        std::vector<MeshLight*> meshLights;

        std::vector<SpotLight*> spotLights;
        std::vector<Material> materials;

        std::vector<Vec3f> vertex_data;
        std::vector<Vec2f> texCoords;

        std::vector<Shape*> meshes;
        std::vector<Sphere*> spheres;
        std::vector<Shape*> triangles;
        
        std::vector<Vec3f> translations;
        std::vector<Vec3f> scalings;
        std::vector<Vec4f> rotations;

        std::vector<Texture*> textures;
        std::vector<Image*> images;
        std::vector<BRDF*> brdfs;

        void loadFromXml(const std::string &filepath);
        void computeFaceNormal(Face& face, Mesh* mesh);
        DorkTracer::Face createFace(int v0idx, int v1idx, int v2idx, DorkTracer::Mesh* mesh);
        void computeFaceBoundingBox(DorkTracer::Face& face, Mesh* mesh);
        void updateBBox(BoundingBox* bbox, Face& face);

        void computeFaceCenter(Face& face, Mesh* mesh);
        void computeFaceProperties(DorkTracer::Face& face, Mesh* mesh);
        void computeTransform(DorkTracer::Shape* shape, std::string input);
        void computeFaceArea(DorkTracer::Face& face, Mesh* mesh);

        void SetupTextures(DorkTracer::Shape* shape, std::string& texturesInp);
        DorkTracer::BoundingBox transformBoundingBox(DorkTracer::BoundingBox original, DorkTracer::Matrix& transform);
        void parseTonemapper(tinyxml2::XMLElement* cameraElm, DorkTracer::Camera* cam);
        void parseBRDFs(tinyxml2::XMLNode* root);
        void parseMaterials(tinyxml2::XMLNode* root);
        void parseLights(tinyxml2::XMLNode* root);
        void parseMeshes(tinyxml2::XMLNode* root, std::string elemName);


    };
}

#endif

 