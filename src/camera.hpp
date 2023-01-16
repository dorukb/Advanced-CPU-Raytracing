#ifndef __HW1__CAMERA__
#define __HW1__CAMERA__

#include <string>

#include "helperMath.h"

namespace DorkTracer
{
    class Tonemapper;

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
            bool hasTonemapper = false;
            std::string imageName;
            
            void SetupDefault(Vec3f pos, Vec3f gazeDir, Vec3f upDir, Vec4f nearPlane, float nearDist, int width, int height, std::string imageName);
            void SetupLookAt(Vec3f pos, Vec3f gazePoint, Vec3f upDir, float nearDist, float fovY, int width, int height, std::string imageName);
            Vec3f GetImagePlanePosition(int width, int height);

            void SetTonemapper(Tonemapper* tonemapper);
            void GetTonemappedImage(int width, int height, float* hdrSrc, unsigned char* ldrDest);
        private:
            Vec3f m_w, m_middle, m_q, m_u;
            float m_left, m_right, m_bottom,m_top;
            Tonemapper* tonemapper;

            void CalculateImagePlaneParams();
            Vec3f GetOrthonormal(Vec3f& toBeTransformed, Vec3f& relative);
    };

}

#endif