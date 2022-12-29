#ifndef __DORKTRACER_PERLINTEXTURE__
#define __DORKTRACER_PERLINTEXTURE__

#include "texture.h"
#include "helperMath.h"
#include <string>
#include <iostream>

#define GRADIENT_COUNT 12
namespace DorkTracer
{
    class PerlinTexture : public Texture
    {

    static int p[512];
    static Vec3f gradients[GRADIENT_COUNT];
    
    public:
        float bumpFactor;

        PerlinTexture(int id, std::string& opMode, float noiseScale, std::string& noiseConversion)
        : Texture(id, opMode)
        {
            this->scale = noiseScale;
            if(noiseConversion == "absval"){
                this->conversionType = Conversion::AbsoluteVal;
            }
            else if(noiseConversion == "linear"){
                this->conversionType = Conversion::Linear;
            }

            // for(int i = 0; i < GRADIENT_COUNT; i++){
            //     gradients[i] = makeUnit(gradients[i]);
            // }
        };
        float GetNormalizer(){
            return 1.0f;
        }
        float GetSampleMultiplier(){
            return bumpFactor;
        }
        virtual Vec3f GetDirectSample(int i, int j)
        {
            // TODO: Fix this inheritance issue, perlin cant support this.
            return Vec3i(180,30,180);
        }
        virtual Vec3f GetRGBSample(float u, float v)
        {
            // TODO: Fix this inheritance issue, perlin cant support this.
            return Vec3i(180,30,180);
        }
        virtual float GetSampleFromWorldPos(float x, float y, float z)
        {
            x *= scale;
            y *= scale;
            z *= scale;

            int X = std::floor(x);
            int Y = std::floor(y);
            int Z = std::floor(z);

            float dx = x - X;
            float dy = y - Y;
            float dz = z - Z;

            X = X & 255;
            Y = Y & 255;
            Z = Z & 255;
            // Select gradients
            int ind0 = p[X+p[Y+p[Z]]] % GRADIENT_COUNT;
            int ind1 = p[X+p[Y+p[Z+1]]] % GRADIENT_COUNT;
            int ind2 = p[X+p[Y+1+p[Z]]] % GRADIENT_COUNT;
            int ind3 = p[X+p[Y+1+p[Z+1]]] % GRADIENT_COUNT;
            int ind4 = p[X+1+p[Y+p[Z]]] % GRADIENT_COUNT;
            int ind5 = p[X+1+p[Y+p[Z+1]]] % GRADIENT_COUNT;
            int ind6 = p[X+1+p[Y+1+p[Z]]] % GRADIENT_COUNT;
            int ind7 = p[X+1+p[Y+1+p[Z+1]]] % GRADIENT_COUNT;

            // Compute dot product of each vector with the corresponding gradient vector
            double c0= dot(gradients[ind0], dx, dy, dz);
            double c1= dot(gradients[ind4], dx-1, dy, dz);
            double c2= dot(gradients[ind2], dx, dy-1, dz);
            double c3= dot(gradients[ind6], dx-1, dy-1, dz);

            double c4= dot(gradients[ind1], dx, dy, dz-1);
            double c5= dot(gradients[ind5], dx-1, dy, dz-1);
            double c6= dot(gradients[ind3], dx, dy-1, dz-1);
            double c7= dot(gradients[ind7], dx-1, dy-1, dz-1);


            double fdx = f(dx);
            double fdy = f(dy);
            double fdz = f(dz);
            double fdx1 = f(dx-1);
            double fdy1 = f(dy-1);
            double fdz1 = f(dz-1);

            double w0 = fdx * fdy * fdz;
            double w1 = fdx1 * fdy * fdz;
            double w2 = fdx * fdy1 * fdz;
            double w3 = fdx1 * fdy1 * fdz;
            double w4 = fdx * fdy * fdz1;
            double w5 = fdx1 * fdy * fdz1;
            double w6 = fdx * fdy1 * fdz1;
            double w7 = fdx1 * fdy1 * fdz1;
            double total = w0*c0 + w1*c1 + w2*c2 + w3*c3 + w4*c4 + w5*c5 + w6*c6 + w7*c7;

            if(total > 1 || total < -1)
            {
                std::cout<<"wrong total? :" << total << std::endl;
            }

            if(conversionType == Conversion::Linear)
            {
                return (total + 1) / 2.0f;
            } 
            else return std::abs(total);
        }

        // Perlin noise is currently generated on-the-fly and so is defined at every point.
        virtual float GetWidth(){
            return INFINITY;
        }
        virtual float GetHeight(){
            return INFINITY;
        }
        virtual bool IsGenerated(){
            return true;
        }

    private:
        float scale;

        enum Conversion
        {
            Linear,
            AbsoluteVal,
        } conversionType;
        
        float dot(const Vec3f &a, float x, float y, float z){
            return a.x*x + a.y*y + a.z*z;
        }

        double mix(double a, double b, double t) {
        return (1-t)*a + t*b;
        }
        
        double f(float x)
        {
            x = std::abs(x);
            if(x > 1)return 0;
            float xSqr = x*x;
            float xCube = xSqr * x;
            return (-6 * xCube*xSqr) + 15 * xCube*x - 10 * xCube + 1;
        }
    };
}

#endif