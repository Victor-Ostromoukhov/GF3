
#include <cstdio>
#include <string>

#include "image_exr.h"


#ifdef EXR
#include <OpenEXR/ImfRgba.h>
#include <OpenEXR/ImfRgbaFile.h>

Image read_image_exr( const char *filename )
{
    using namespace Imf;
    using namespace Imath;

    RgbaInputFile file(filename);
    Box2i dw= file.dataWindow();

    int width= dw.max.x - dw.min.x + 1;
    int height= dw.max.y - dw.min.y + 1;

    std::vector<Rgba> pixels(width * height);
    file.setFrameBuffer(pixels.data() - dw.min.x - dw.min.y * width, 1,  width);
    file.readPixels(dw.min.y, dw.max.y);

    Image image(width, height);
    int offset= 0;
    //~ for(int y= 0; y < image.height(); y++) 
    for(int y= image.height() -1; y >= 0; y--) 
    for(int x= 0; x < image.width(); x++, offset++) 
    {
        image(x, y)= Color32(pixels[offset].r, pixels[offset].g, pixels[offset].b, pixels[offset].a);
    }
    
    printf("loading exr image '%s'...\n", filename);
    
    return image;
}

int write_image_exr( const Image& image, const char *filename )
{
    using namespace Imf;
    using namespace Imath;

    std::vector<Rgba> pixels(image.size());
    int offset= 0;
    //~ for(int y= 0; y < image.height(); y++)
    for(int y= image.height() -1; y >= 0 ; y--)
    for(int x= 0; x < image.width(); x++, offset++)
    {
        Color32 pixel= image(x, y);
        pixels[offset]= Rgba(pixel.r, pixel.g, pixel.b, pixel.a);
    }

    Box2i displayWindow(V2i(0, 0), V2i(image.width() -1, image.height() -1));
    Box2i dataWindow(V2i(0, 0), V2i(image.width() -1, image.height() -1));

    RgbaOutputFile file(filename, displayWindow, dataWindow, WRITE_RGBA);
    file.setFrameBuffer(pixels.data(), 1, image.width());
    file.writePixels(image.height());
    
    printf("writing exr image '%s'...\n", filename);
    return 0;
}

#else
Image read_image_exr( const char *filename )
{
    printf("[error] no openEXR support... can't read '%s' !\n", filename);
    return Image();
}

int write_image_exr( const Image& image, const char *filename )
{
    printf("[error] no openEXR support... can't write '%s'\n", filename);
    return -1;
}
#endif

//! renvoie vrai si le nom de fichier se termine par .hdr.
bool is_exr_image( const char *filename )
{
    return (std::string(filename).rfind(".exr") != std::string::npos);
}

