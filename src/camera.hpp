#ifndef __HW1__CAMERA__
#define __HW1__CAMERA__

#include "helperMath.h"
#include <string>


namespace DorkTracer
{
    class Camera
    {
        public:
            Vec3f position;
            Vec3f gaze;
            Vec3f up;
            Vec3f right;

            float nearDist;
            int imageWidth;
            int imageHeight;
            int samplesPerPixel;
            float focusDistance;
            float apertureSize;
            std::string imageName;
            
            void SetupDefault(Vec3f pos, Vec3f gazeDir, Vec3f upDir, Vec4f nearPlane, float nearDist, int width, int height, std::string imageName);
            void SetupLookAt(Vec3f pos, Vec3f gazePoint, Vec3f upDir, float nearDist, float fovY, int width, int height, std::string imageName);
            Vec3f GetImagePlanePosition(int width, int height);
        private:
            Vec3f m_w, m_middle, m_q, m_u;
            float m_left, m_right, m_bottom,m_top;

            void CalculateImagePlaneParams();
            Vec3f GetOrthonormal(Vec3f& toBeTransformed, Vec3f& relative);
    };

}

#endif