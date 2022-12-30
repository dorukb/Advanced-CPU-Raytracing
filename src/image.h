#ifndef __DORKTRACER_IMAGE__
#define __DORKTRACER_IMAGE__

#include "stb_image.h"
#include <string>
#include "helperMath.h"
#include "tinyexr.h"

namespace DorkTracer
{
    class Image
    {

    public:
        int id;
        int width;
        int height;

        virtual Vec3f GetSample(int i, int j) = 0;
        virtual int GetSampleGreyscale(int i, int j) = 0;

    protected:
        int channels;
        virtual void LoadImage(std::string filename) = 0;

    };
}
#endif