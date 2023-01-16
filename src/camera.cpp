#include "camera.hpp"
#include "tonemapper.h"
using namespace DorkTracer;

void Camera::SetupDefault(Vec3f pos, Vec3f gazeDir, Vec3f upDir, Vec4f nearPlane, float nearDist, int width, int height, std::string imageName){
    this->position = pos;
    this->nearDist = nearDist;
    this->imageWidth = width;
    this->imageHeight = height;
    this->imageName = imageName;

    m_left = nearPlane.x;
    m_right = nearPlane.y;
    m_bottom = nearPlane.z;
    m_top = nearPlane.w;

    this->gaze = makeUnit(gazeDir);
    Vec3f tempUp = makeUnit(upDir);

    // make sure cam.up is orthogonal to cam.gaze.
    this->up = makeUnit(GetOrthonormal(tempUp, gaze));

    CalculateImagePlaneParams();
}
void Camera::SetupLookAt(Vec3f pos, Vec3f gazePoint, Vec3f upDir, float nearDist, float fovY, int width, int height, std::string imageName){
    this->position = pos;
    this->nearDist = nearDist;
    this->imageWidth = width;
    this->imageHeight = height;
    this->imageName = imageName;

    float aspect = (float)width/height;
    
    m_top = nearDist * std::tan((fovY *(M_PI / 180.0f) / 2.0f));
    m_right = m_top * aspect;

    m_bottom = -m_top;
    m_left = -m_right;

    this->gaze = makeUnit(gazePoint - pos);
    Vec3f tempUp = makeUnit(upDir);

    // make sure cam.up is orthogonal to cam.gaze.
    Vec3f tempRight = makeUnit(cross(tempUp, gaze));
    this->up = makeUnit(cross(gaze, tempRight));
    
    CalculateImagePlaneParams();
}

Vec3f Camera::GetOrthonormal(Vec3f& toBeTransformed, Vec3f& relative){

    // return orthogonal version of toBeTransformed wrto relative.
    float dotvu = dot(toBeTransformed, relative);
    float relativeSqr = dot(relative, relative);
    Vec3f proj = relative * (dotvu / relativeSqr);
    Vec3f orthogonalUp = toBeTransformed - proj;
    return orthogonalUp;
}

void Camera::CalculateImagePlaneParams()
{
    // calculate and cache internal image plane parameters.
    m_middle = position + gaze * nearDist;
    
    m_w = -gaze;
    right = cross(up, m_w);

    // Up = v,  Gaze = −w, u = v ×w

    m_middle = position + gaze * nearDist;
    m_q = m_middle + right * m_left + up * m_top;
}

Vec3f Camera::GetImagePlanePosition(int width, int height)
{
    float su = (width + 0.5) * (m_right - m_left) / imageWidth;
    float sv = (height + 0.5) * (m_top- m_bottom) / imageHeight;

    return m_q + right * su + up *-sv;
}

void Camera::SetTonemapper(Tonemapper* tonemapper)
{
    this->tonemapper = tonemapper;
    this->hasTonemapper = tonemapper != nullptr;
}


void Camera::GetTonemappedImage(int width, int height, float* hdrSrc, unsigned char* ldrDest)
{
    this->tonemapper->Tonemap(width, height, hdrSrc, ldrDest);
}