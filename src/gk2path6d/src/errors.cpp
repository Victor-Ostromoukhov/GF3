
#include <cstdio>
#include <cstring>
#include <cfloat>
#include <algorithm>


#ifdef TTF
#include "SDL2/SDL_ttf.h"
#endif

#include "image.h"
#include "image_io.h"
#include "image_hdr.h"
#include "image_exr.h"


struct cell
{
    Image image;
    
    const char *filename;
    const char *title;

    int x, y;
    int nspp;
};


// 
int nspp( const char *filename )
{
    int n= -1;
    
    // suppose que le nom de fichier est de la forme 'prefix-nspp.[pfm|exr|hdr]'
    const char *ext= strrchr(filename, '.');
    if(ext)
    {
        if(strcmp(ext, ".pfm") == 0 || strcmp(ext, ".hdr") == 0 || strcmp(ext, ".exr") == 0)
        {
            const char *last= strrchr(filename, '-');
            if(sscanf(last, "-%d.", &n) == 1)
                return n;
        }
    }
    
    printf("\n[error] filename '%s'... doesn't match 'prefix-nspp.[pfm|hdr|exr]'\n", filename);
    return 0;
}


// mse d'une vignette 
std::pair<double, double> Mlocal( /* ref */ const Image &A, /* image */ const Image &B,
    int mX, int mY,
    int MX, int MY ,
    double &vmin, double &vmax)
{
    double val= 0.0;
    double var= 0.0;
    vmax= -1.0;
    vmin= 100000.0;
    
    int cpt= 0;
    for(int i=mX; i <= MX; ++i)
    for(int j=mY; j <= MY; ++j)
    {
        Color32 d= A(i, j) - B(i, j);
        double v= (std::abs(d.r) + std::abs(d.g) + std::abs(d.b)) / 3;
        val += v;
        vmax= std::max(vmax, v );
        vmin= std::min(vmin, v );
        var += v*v;
        cpt++;
    }
    
    val /= double(cpt);
    var /= double(cpt);
    var= (var - val*val);
    return std::make_pair(val, var);
}


// comparaison des mse des samplers
struct mse_less
{
    const std::vector<double>& mse;
    mse_less( const std::vector<double>& values ) : mse(values) {}
    
    bool operator() ( const int& a, const int& b ) const
    {
        return mse[a] < mse[b];
    }
};


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


// parula / matplotlib palette, cf http://www.gnuplotting.org/matlab-colorbar-parula-with-gnuplot/
Color32 palette[]= 
{
    Color32(0x00, 0x72, 0xbd, 255) * (1 / float(255)),
    Color32(0xd9, 0x53, 0x19, 255) * (1 / float(255)),
    Color32(0xed, 0xb1, 0x20, 255) * (1 / float(255)),
    Color32(0x7e, 0x2f, 0x8e, 255) * (1 / float(255)),
    Color32(0x77, 0xac, 0x30, 255) * (1 / float(255)),
    Color32(0x4d, 0xbe, 0xee, 255) * (1 / float(255)),
    Color32(0xa2, 0x14, 0x2f, 255) * (1 / float(255))
};


int main( int argc, const char *argv[] )
{
    if(argc < 2)
    {
        printf("usage: %s output_prefix reference.[hdr|pfm|exr] --line 1 [--title \"...\"] image.[hdr|pfm|hdr|png] [--line 2 images --line 3 images] \n", argv[0]);
        return 0;
    }
    
    if(TTF_Init() < 0)
    {
        printf("TTF_Init: %s\n", TTF_GetError());
        return 1;
    }
    
    TTF_Font *font= TTF_OpenFont( "data/roboto.ttf", 64 );
    if(!font)
    {
        printf("TTF_OpenFont: %s\n", TTF_GetError());
        return 1;
    }
    
    //
    const char *prefix= "label";
    if(argc > 2)
        prefix= argv[1];
        
    Image reference;
    if(argc > 3)
    {
        if(is_exr_image(argv[2]))
            reference= read_image_exr(argv[2]);
        else if(is_pfm_image(argv[2]))
            reference= read_image_pfm(argv[2]);
        else if(is_hdr_image(argv[2]))
            reference= read_image_hdr(argv[2]);
    }
    if(reference.size() == 0)
    {
        printf("[error] loading reference image '%s'...\n", argv[2]);
        return 1;
    }
    
    // charge les images
    int rows= 0;
    int columns= 0;
    std::vector< std::vector<cell> > grid;
    {
        std::vector<cell> cells;
        
        int x= 0;
        int y= -1;
        int last_row= -1;
        const char *last_title= nullptr;
        
        for(int option= 3; option < argc; )
        {
            if(argv[option][0] == '-')
            {
                // nouvelle ligne
                if(strcmp(argv[option], "--line") == 0 && option +1 < argc)
                {
                    int row= 0;
                    if(sscanf(argv[option+1], "%d", &row) != 1)
                        break;
                    
                    x= 0;
                    //~ y++;
                    y= row -1;
                    option+= 2;
                    
                    printf("line %d\n", row);
                }
                else if(strcmp(argv[option], "--title") == 0 && option +1 < argc)
                {
                    last_title= argv[option+1];
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
                int n= nspp(filename);
                cells.push_back( {Image(), filename, last_title, x, y, n} );
                
                printf("x %d y %d: '%s' %s\n", x, y, last_title, filename);
                
                x++;
                option++;
            }
        }
        
        // chargement parallele des images
    #pragma omp parallel for
        for(int i= 0; i < int(cells.size()); i++)
            cells[i].image= read(cells[i].filename);
        
        // dimensions de la grille
        for(auto const & cell : cells)
        {
            rows= std::max(rows, cell.y +1);
            columns= std::max(columns, cell.x +1);
        }
        
        // cree la grille
        grid.resize(rows);
        for(int i= 0; i < rows; i++)
            grid[i].resize(columns);
        
        // copie les images a leur place
        for(auto const & cell : cells)
        {
            int x= cell.x;
            int y= cell.y;
            assert(x < columns);
            assert(y < rows);
            grid[y][x]= cell;   // move ?
        }
    }
    
    
    // evalue les erreurs... et classe les samplers
    const int size= 1;
#if 1
    for(int c= 0; c < columns; c++)
    {
        assert(grid[0][c].image.size());
        int width= grid[0][c].image.width();
        int height= grid[0][c].image.height();
        
        Image score(width, height);
        
    #pragma omp parallel for
        for(int y= size; y < height - size; y++)
        for(int x= size; x < width - size; x++)
        {
            std::vector<double> mse(rows, DBL_MAX);
            std::vector<int> index(rows);
            
            // erreurs des differents samplers (1 par ligne de la grille)
            for(int r= 0; r < rows; r++)
            {
                if(grid[r][c].image.size() == 0)
                    continue;
                
                assert(grid[r][c].image.width() == width);
                assert(grid[r][c].image.height() == height);
                
                if(grid[r][c].image(x, y).grey() == 0)
                    continue;
                
                double v0, v1;
                std::pair<double, double> stat= Mlocal(reference, grid[r][c].image, x - size, y - size, x + size, y + size, v0, v1);      // Mean
                mse[r]= stat.first;
            }
            
            // classer les samplers
            for(unsigned int i= 0; i < index.size(); i++)
                index[i]= i;
            
            std::sort(index.begin(), index.end(), mse_less(mse));
            
            if(mse[index[0]] < DBL_MAX)
                score(x, y)= palette[index[0]];
        }
        
        
        
        // legende
        for(int i= 0; i < rows; i++)
        {
            if(!grid[i][0].title)
                continue;
            
            // pave de couleur
            for(int y= 0; y < 64; y++)
            for(int x= 0; x < 64; x++)
                score(x, i*64 + y+8)= Color32(palette[i], 1);
            
            // texte
            SDL_Color color= { 255, 255, 255, 255 };
            SDL_Surface *text= TTF_RenderText_Blended(font, grid[i][0].title, color);
            if(text)
            {
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
                        
                        score(x+64+8, i*64 + y+8)= Color32(score(x+64+8, i*64 + y+8) * Color32(1 - float(a) / 255) + Color32(float(a) / 255), 1);
                        pixel= pixel + format.BytesPerPixel;
                    }
                }
                
                SDL_FreeSurface(text);
            }
        }
        
        char tmp[1024];
        assert(grid[0][c].nspp != 0);
        sprintf(tmp, "%s_size%d-%05d.png", prefix, size, grid[0][c].nspp);
        //~ printf("writing '%s'...\n", tmp);
        write_image(score, tmp);
    }
    
    TTF_Quit();
#endif

    // mse global
    {
        for(int r= 0; r < rows; r++)
        {
            if(!grid[r][0].title)
                continue;
            
            printf("mse %s\n", grid[r][0].title);
            std::string filename;
            {
                char prefix[1024];
                strcpy(prefix, grid[r][0].filename);
                // suppose que le nom de fichier est de la forme 'prefix-nspp.[pfm|exr|hdr]'
                char *ext= strrchr(prefix, '.');
                if(ext)
                {
                    if(strcmp(ext, ".pfm") == 0 || strcmp(ext, ".hdr") == 0 || strcmp(ext, ".exr") == 0)
                    {
                        char *last= strrchr(prefix, '-');
                        int n= 0;
                        if(sscanf(last, "-%d.", &n) == 1)
                            *last= 0;
                    }
                }
                
                filename= change_extension(change_prefix("mse_", prefix), ".txt");
            }          
            
            //~ printf("mse file '%s'\n", filename.c_str());
            FILE *out= fopen(filename.c_str(), "wt");
            if(out)
            {
                printf("writing '%s'...\n", filename.c_str());
                
                fprintf(out, "# %s\n", grid[r][0].title);
                for(int c= 0; c < columns; c++)
                {
                    if(grid[r][c].image.size() == 0)
                        continue;
                        
                    double v0, v1;
                    std::pair<double, double> stat= Mlocal(reference, grid[r][c].image, 0, 0, grid[r][c].image.width()-1, grid[r][c].image.height()-1, v0, v1);      // Mean
                    
                    double mse= stat.first;
                    fprintf(out, "%d %lf\n", grid[r][c].nspp, mse);
                }
                fclose(out);
            }
        }
    }
    
    return  0;
}
