#include <string>

namespace DorkTracer
{
    class Tonemapper
    {
    public:
        std::string operatorType;
        float keyValue;
        float burnPerct;
        float saturation;
        float gamma;

        Tonemapper(std::string opType, float keyVal, float burnPerct, float saturation, float gamma)
        {
            this->operatorType = opType;
            this->keyValue = keyVal;
            this->burnPerct = burnPerct;
            this->saturation = saturation;
            this->gamma = gamma;
        }
        
        Vec3i GetTonemappedColor(Vec3f color)
        {
            // TODO: Implement tonemapping
            return Vec3i(color.x, color.y, color.z);
        }

    };
}
