namespace DorkTracer
{ 
    
    enum class Textures
    {
        Diffuse,
        Bump,
        Normal,
    };

    enum class Interpolation
    {
        Nearest_Neighbor,
        Bilinear,
    };

    enum class DecalMode
    {
        Replace_kd,
        Blend_kd,
        Replace_ks,
        Replace_bg,
        Replace_normal,
        Bump_normal,
        Replace_all, // replace all components (i.e. diffuse, specular, and ambient) of the surface shading color with this texture’s value.

    };

    class Texture
    {
   
    public:
        Textures type;
        int id;

        Texture(Textures type)
        {
            this->type = type;
        };
        
    };

  
}
