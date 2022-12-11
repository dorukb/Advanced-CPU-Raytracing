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
#include "texture.h"

#include "image.h"

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
        int max_recursion_depth;
        bool isMotionBlurEnabled;
        Vec3f ambient_light;

        std::vector<Camera> cameras;
        std::vector<PointLight> point_lights;
        std::vector<AreaLight*> areaLights;
        std::vector<Material> materials;
        std::vector<Vec3f> vertex_data;
        std::vector<Shape*> meshes;
        std::vector<Sphere*> spheres;
        std::vector<Shape*> triangles;
        
        std::vector<Vec3f> translations;
        std::vector<Vec3f> scalings;
        std::vector<Vec4f> rotations;

        std::vector<Texture> textures;
        std::vector<Image*> images;

        void loadFromXml(const std::string &filepath);
        void computeFaceNormal(Face& face, std::vector<Vec3f>& vertices);
        void computeFaceBoundingBox(DorkTracer::Face& face, std::vector<Vec3f>& vertices);
        void computeFaceCenter(Face& face, std::vector<Vec3f>& vertices);
        void computeFaceProperties(DorkTracer::Face& face, std::vector<Vec3f>& vertices);
        void computeTransform(DorkTracer::Shape* shape, std::string input);
        DorkTracer::BoundingBox transformBoundingBox(DorkTracer::BoundingBox original, DorkTracer::Matrix& transform);
    };
}

#endif

 