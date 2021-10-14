#include <cstdio>

#include "../includes/SeedTileHelper.hpp"

int main( int argc, char *argv[] )
{
    if(argc < 3)
    {
        printf("usage: %s output seedtile_12 [seedtile_23]\n", argv[0]);
        return 0;
    }
    
    const char *filename= "seeds.dat";
    if(argc > 1)
        filename= argv[1];
    
    SeedTile output;
    if(argc > 2)
        output.loadTile(argv[2]);
    
    assert(output.dimension > 0);
    for(int option= 3; option < argc; option++)
    {
        SeedTile input;
        input.loadTile(argv[option]);
        assert(input.dimension > 0);
        assert(input.size == output.size);
        
        SeedTile tmp(output.size, output.dimension + input.dimension);
        for(int y= 0; y < output.size; y++)
        for(int x= 0; x < output.size; x++)
        {
            for(int d= 0; d < output.dimension; d++)
                tmp(x, y, d)= output(x, y, d);
            for(int d= 0; d < input.dimension; d++)
                tmp(x, y, output.dimension + d)= input(x, y, d);
        }
        
        output= tmp;
        printf("output: %dx%d %dd\n", output.size, output.size, output.dimension);
    }
    
    printf("writing tilemap '%s'...\n", filename);
    output.saveTile(filename);
    return 0;
}
