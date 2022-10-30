#ifndef __helperMath_h__
#define __helperMath_h__

#include <math.h>

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
Vec3f cross(const Vec3f &first, const Vec3f &second);
float dot(const Vec3f &a, const Vec3f &b);

float len(Vec3f a);
Vec3f makeUnit(Vec3f a);
float determinant(float m [3][3]);

Vec3f operator+(const Vec3f& v1, const Vec3f& v2);
Vec3f operator-(const Vec3f& v1, const Vec3f& v2);
Vec3f operator*(const Vec3f& v1, const Vec3f& v2);

Vec3f operator*(const Vec3f& v1, float scaler);
Vec3f operator/(const Vec3f& v1, float scaler);

Vec3i clamp(Vec3f rgb);
int clamp(int x);

#endif // __helper_h__

