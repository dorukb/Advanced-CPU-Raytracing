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


namespace DorkTracer{

    struct PointLight
    {
        Vec3f position;
        Vec3f intensity;
    };

    struct Scene
    {

        Vec3i background_color;
        float shadow_ray_epsilon;
        int max_recursion_depth;
        std::vector<Camera> cameras;
        Vec3f ambient_light;
        std::vector<PointLight> point_lights;
        std::vector<Material> materials;
        std::vector<Vec3f> vertex_data;
        std::vector<Shape*> meshes;
        std::vector<Sphere*> spheres;
        std::vector<Shape*> triangles;
        
        std::vector<Vec3f> translations;
        std::vector<Vec3f> scalings;
        std::vector<Vec4f> rotations;

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

 