#ifndef __DORKTRACER_GENERATEDTEXTURE__
#define __DORKTRACER_GENERATEDTEXTURE__

#include "texture.h"

namespace DorkTracer
{
    class GeneratedTexture : Texture
    {
        
    public:
        enum GenerationType{
            Perlin,
            Checkboard,
        };


        GeneratedTexture(GenerationType genType);

        
   
    };
}

#endif