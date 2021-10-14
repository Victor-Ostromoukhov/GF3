
#include <cstdio>
#include <cstring>

#ifdef TTF
#include "SDL2/SDL_ttf.h"
#endif

#include "window.h"
#include "image.h"
#include "image_io.h"
#include "image_hdr.h"
#include "image_exr.h"


Image read( const char *filename )
{
    Image image;
    if(is_exr_image(filename))
        image= read_image_exr(filename);
    else if(is_pfm_image(filename))
        image= read_image_pfm(filename);
    else if(is_hdr_image(filename))
        image= read_image_hdr(filename);
    else
        image= read_image(filename);
    
    return image;
}

struct cell
{
    int x, y;
    int w, h;
    Image image;
    const char *title;
};

// copie src a la position x, y, dans image
void blit( Image& image, const int dstx, const int dsty, const Image& src )
{
    for(int y= 0; y < src.height(); y++)
    for(int x= 0; x < src.width(); x++)
        image(dstx+x, dsty+y)= src(x, y);
}


int main( int argc, const char *argv[] )
{
    
    if(argc < 2)
    {
        printf("usage: %s output.[hdr|pfm|exr|png] --line 1 [--title \"...\"] image.[hdr|pfm|hdr|png] [--line 2 images --line 3 images] \n", argv[0]);
        return 0;
    }

    int x= 0;
    int y= -1;
    std::vector<cell> cells;
    
    int last_row= -1;
    const char* last_title= nullptr;
    
    const char *output_filename= nullptr;
    if(argc > 2)
        output_filename= argv[1];
    
    // lire les images 
    for(int option= 2; option < argc; )
    {
        int row= 0;
        if(argv[option][0] == '-')
        {
            // nouvelle ligne
            if(strcmp(argv[option], "--line") == 0 && option +1 < argc)
            {
                if(sscanf(argv[option+1], "%d", &row) != 1)
                    break;
                x= 0;
                y++;
                option+= 2;
                
                printf("line %d\n", row);
            }
            else if(strcmp(argv[option], "--title") == 0 && option +1 < argc)
            {
            #ifdef TTF
                last_row= row;
                last_title= argv[option+1];
                printf("title '%s'\n", argv[option+1]);
            #endif
                option+= 2;
            }
            else
            {
                printf("[error] parsing option '%s'...\n", argv[option]);
                return 1;
            }
        }
        else
        {
            const char *filename= argv[option];
            Image image= read(filename);
            if(image.size() == 0)
            {
                printf("[error] loading image'%s'...\n", argv[option]);
                return 1;
            }
            
            cells.push_back( {x, y, image.width(), image.height(), image, last_title} );
            
            if(last_title) last_title= nullptr;
            x++;
            option++;
        }
    }
    
    int width= 0;
    int height= 0;
    int rows= 0;
    int columns= 0;
    
    for(int i= 0; i < cells.size(); i++)
    {
        rows= std::max(rows, cells[i].y +1);
        columns= std::max(columns, cells[i].x +1);
        
        width= std::max(width, cells[i].x*cells[i].w + cells[i].w);
        height= std::max(height, cells[i].y*cells[i].h + cells[i].h);
    }
    
#ifdef TTF
    // titre
    if(TTF_Init() < 0)
    {
        printf("TTF_Init: %s\n", TTF_GetError());
        return 1;
    }
    
    TTF_Font *font= TTF_OpenFont( smart_path("data/roboto.ttf"), 96 );
    if(!font)
    {
        printf("TTF_OpenFont: %s\n", TTF_GetError());
        return 1;
    }
    
    for(int i= 0; i < cells.size(); i++)
    {
        if(!cells[i].title)
            continue;
        
        SDL_Color color= { 255, 255, 255, 255 };
        SDL_Surface *text= TTF_RenderText_Blended(font, cells[i].title, color);
        if(text)
        {
            // composition du texte dans l'image...
            const SDL_PixelFormat format= *text->format;
            assert(format.BitsPerPixel == 32);
            
            int py= 0;
            for(int y= text->h -1; y >= 0; y--, py++)
            {
                Uint8 *pixel= (Uint8 *) text->pixels + py * text->pitch;
                for(int x= 0; x < text->w; x++)
                {
                    Uint8 r= pixel[format.Rshift / 8];
                    Uint8 g= pixel[format.Gshift / 8];
                    Uint8 b= pixel[format.Bshift / 8];
                    Uint8 a= pixel[format.Ashift / 8];
                    
                    //~ grid(x, y)= Color32(float(a) / 255);
                    cells[i].image(x+8, y+8)= Color32(cells[i].image(x+8, y+8) * Color32(1 - float(a) / 255) + Color32(float(a) / 255), 1);
                    pixel= pixel + format.BytesPerPixel;
                }
            }
            
            SDL_FreeSurface(text);
        }
    }
    
    TTF_Quit();
#endif

    // copie les images dans la grille
    Image grid(width, height);
#pragma omp parallel for
    for(int i= 0; i < cells.size(); i++)
        blit(grid, cells[i].x*cells[i].w, (rows -1 - cells[i].y)*cells[i].h, cells[i].image);
    
    printf("grid: %d rows, %d columns, '%s' %dx%d pixels\n", rows, columns, output_filename, width, height);
    
    
    if(is_exr_image(output_filename))
        write_image_exr(grid, output_filename);
    else if(is_pfm_image(output_filename))
        write_image_pfm(grid, output_filename);
    else if(is_hdr_image(output_filename))
        write_image_hdr(grid, output_filename);
    else
        write_image(grid, output_filename);
    
    return  0;
}
