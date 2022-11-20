#include "math.h"

struct Gaussian2D{

    float sigma;
    float sigmaSqr;
    float c1;

    Gaussian2D(float sigma){
        this->sigma = sigma;
        this->sigmaSqr = sigma * sigma;
        this->c1 = 1 / (2 * M_PI * sigmaSqr);
    }

    float GetWeight(int x, int y)
    {
        float exp = -0.5 * ((x*x + y*y) / sigmaSqr);
        return c1 * pow(M_E, exp);
    };
};

