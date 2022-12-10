#ifndef __helperMath_h__
#define __helperMath_h__

#include <math.h>
#include <stdint.h>

struct Vec3i
{
    int x, y, z;
};
struct Vec3f
{
    float x, y, z;
    Vec3f(){
        x = y = z = 0.0f;
    };
    Vec3f(Vec3i v){
        x = v.x;
        y = v.y;
        z = v.z;
    };
    
    Vec3f(float x, float y, float z){
        this->x = x;
        this->y = y;
        this->z = z;
    };

    Vec3f operator-() const {
        Vec3f result;
        result.x = x * -1.0f;
        result.y = y * -1.0f;
        result.z = z * -1.0f;
        return result;
    }
    float operator[](uint32_t idx) const{
        switch(idx){
            case 0:
                return x;
            case 1:
                return y;
            case 2:
                return z;
            _:
                return x;
        }
    }
};


struct Vec2f{
    float x,y;
};

struct Vec4f
{
    float x, y, z, w;
    Vec4f(){}
    Vec4f(Vec3f& v, float wVal){
        x = v.x;
        y = v.y;
        z = v.z;
        w = wVal;
    }
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

Vec3i operator/(const Vec3i& v1, float scaler);
Vec3i clamp(Vec3f rgb);
int clamp(int x);
void GetOrthonormalBasis(Vec3f r, Vec3f& u, Vec3f& v);

#endif // __helper_h__

