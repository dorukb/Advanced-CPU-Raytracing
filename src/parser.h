#ifndef __HW1__PARSER__
#define __HW1__PARSER__

#include <string>
#include <vector>
#include "happly.h"

namespace parser
{
    enum MaterialType{
        Mirror,
        Dielectric,
        Conductor,
        Default
    };
    //Notice that all the structures are as simple as possible
    //so that you are not enforced to adopt any style or design.
    struct Vec3f
    {
        float x, y, z;
        Vec3f operator-() const {
            Vec3f result;
            result.x = x * -1.0f;
            result.y = y * -1.0f;
            result.z = z * -1.0f;
            return result;
        }
    };

    struct Vec3i
    {
        int x, y, z;
    };

    struct Vec4f
    {
        float x, y, z, w;
    };

    struct Camera
    {
        Vec3f position;
        Vec3f gaze;
        Vec3f up;
        Vec4f near_plane;

        bool isLookAt;
        Vec3f gazePoint;
        float fovY;

        float near_distance;
        int image_width, image_height;
        std::string image_name;
    };

    struct PointLight
    {
        Vec3f position;
        Vec3f intensity;
    };

    struct Material
    {
        Vec3f ambient;
        Vec3f diffuse;
        Vec3f specular;
        Vec3f mirror;
        float phong_exponent;
        MaterialType type;
        float refractiveIndex;
        Vec3f absorptionCoefficient;
        float conductorAbsorptionIndex;
    };

    struct Face
    {
        int v0_id;
        int v1_id;
        int v2_id;
        Vec3f n;
    };

    struct Mesh
    {
        bool useOwnVertices;
        int material_id;
        std::vector<Face> faces;
        std::vector<Vec3f> vertices;
    };

    struct Triangle
    {
        int material_id;
        Face indices;
    };

    struct Sphere
    {
        int material_id;
        int center_vertex_id;
        float radius;
    };

    struct Scene
    {
        //Data
        Vec3i background_color;
        float shadow_ray_epsilon;
        int max_recursion_depth;
        std::vector<Camera> cameras;
        Vec3f ambient_light;
        std::vector<PointLight> point_lights;
        std::vector<Material> materials;
        std::vector<Vec3f> vertex_data;
        std::vector<Mesh> meshes;
        std::vector<Triangle> triangles;
        std::vector<Sphere> spheres;

        //Functions
        void loadFromXml(const std::string &filepath);
        void computeFaceNormal(Face& face, std::vector<Vec3f>& vertices);

    };
}

#endif
