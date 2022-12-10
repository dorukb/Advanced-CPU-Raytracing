#include "math.h"

struct Gaussian2D{

    float sigma;
    float sigmaSqr;
    float c1;
    float boundaryValue;

    Gaussian2D(float sigma){
        this->sigma = sigma;
        this->sigmaSqr = sigma * sigma;
        this->c1 = 1.0f / (2.0f * M_PI * sigmaSqr);
    }

    float GetWeight(float x, float y)
    {
        float exponent = -0.5 * ((x*x + y*y) / sigmaSqr);
        return c1 * std::exp(exponent);
    };
};

