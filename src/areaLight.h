#include "helperMath.h"
#include <random>

namespace DorkTracer{
    class AreaLight
    {

    public:
        int id;
        Vec3f position;
        Vec3f normal;
        Vec3f radiance;
        float extent;
        float area;

        AreaLight(int id, Vec3f pos, Vec3f normal, Vec3f radiance, float extent)
        {
            this->id = id;
            this->position = pos;
            this->normal = normal;
            this->radiance = radiance;
            this->extent = extent;
            this->area = extent * extent;

            this->sampleGen = std::mt19937();
            this->sampleDistro = std::uniform_real_distribution<>(-0.5f, 0.5f);

            GetOrthonormalBasis(this->normal, u, v);
        }
        
        Vec3f GetSample()
        {
            float offsetU = sampleDistro(sampleGen);
            float offsetV = sampleDistro(sampleGen);

            return position + (u * (extent * offsetU)) + (v * (extent * offsetV));
        };


        Vec3f GetSamplePosition(float randOffsetX, float randOffsetY)
        {
            return position + (u * (extent * randOffsetX)) + (v * (extent * randOffsetY));
        };

    private:
        std::mt19937 sampleGen;
        std::uniform_real_distribution<> sampleDistro;
        Vec3f u,v;
    };

}
