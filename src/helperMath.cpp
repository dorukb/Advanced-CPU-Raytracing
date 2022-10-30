#include "helperMath.h"
#include <algorithm>

// using namespace parser;

Vec3f operator+(const Vec3f& v1, const Vec3f& v2){
    Vec3f result;
    result.x = v1.x + v2.x;
    result.y = v1.y + v2.y;
    result.z = v1.z + v2.z;
    return result;
}

Vec3f operator-(const Vec3f& v1, const Vec3f& v2){
    Vec3f result;
    result.x = v1.x - v2.x;
    result.y = v1.y - v2.y;
    result.z = v1.z - v2.z;
    return result;
}


Vec3f operator*(const Vec3f& v1, float scaler){
    Vec3f result;
    result.x = v1.x * scaler;
    result.y = v1.y * scaler;
    result.z = v1.z * scaler;
    return result;
}
Vec3f operator*(const Vec3f& v1, const Vec3f& v2){
    Vec3f result;
    result.x = v1.x * v2.x;
    result.y = v1.y * v2.y;
    result.z = v1.z * v2.z;
    return result;
}

Vec3f operator/(const Vec3f& v1, float scaler){
    Vec3f result;
    result.x = v1.x / scaler;
    result.y = v1.y / scaler;
    result.z = v1.z / scaler;
    return result;
}

float dot(const Vec3f &a, const Vec3f &b)
{
    float result = a.x*b.x + a.y*b.y + a.z*b.z;
    return result;
}

Vec3f cross(const Vec3f &first, const Vec3f &second){
    Vec3f result;
    result.x = first.y * second.z - first.z * second.y;
    result.y = first.z * second.x - first.x * second.z;
    result.z = first.x * second.y - first.y * second.x;
    return result;
}

float len(Vec3f a)
{
    return sqrtf((a.x*a.x)+(a.y*a.y)+(a.z*a.z));
}
Vec3f makeUnit(Vec3f a)
{
    Vec3f result;
    float l = len(a);
    result.x = a.x/l;
    result.y = a.y/l;
    result.z = a.z/l;
    return result;
}

float determinant(float m [3][3])
{
    float firstTerm =  m[0][0] * (m[1][1] * m[2][2] - m[1][2]* m[2][1]);
    float secondTerm = m[1][0] * (m[0][2] * m[2][1] - m[0][1]* m[2][2]);
    float thirdTerm =  m[2][0] * (m[0][1] * m[1][2] - m[1][1]* m[0][2]);
    return firstTerm + secondTerm + thirdTerm;
}

Vec3i clamp(Vec3f rgb)
{
    Vec3i result;
    result.x = clamp((int) rgb.x);
    result.y = clamp((int) rgb.y);
    result.z = clamp((int) rgb.z);
    return result;
}

int clamp(int x)
{
    return std::min(255, std::max(x, 0));
}

