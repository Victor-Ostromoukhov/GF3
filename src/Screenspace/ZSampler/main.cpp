
#include <cstdint>
#include <cstdio>
#include <cassert>

#include <vector>
#include <algorithm>

#include "zhash.h"
#include "z.h"

const int nspp= 4;
const int bounds= 4;
const int ndims= 1; // paires de dimensions

int main( )
{
    int samplesPerPixel = RoundUpPow2(nspp);
    int resolution = RoundUpPow2(bounds);
    int log2Resolution = Log2Int(resolution);
    
    // Z hash
    {
        FILE *out= fopen("zhash.dat", "wt");
        assert(out);
        for(int y= 0; y < bounds; y++)
        for(int x= 0; x < bounds; x++)
        {
            std::vector<Float> samples2D(2*nspp*ndims);
            for(int d = 0; d < ndims; ++d) 
            {
                Z2DHash(x, y, 1, samplesPerPixel, log2Resolution, &samples2D[2*d*nspp], d);
                
                printf("pixel (%d, %d) dim %d\n", x, y, d);
                for(int s= 0; s < nspp; s++)
                {
                    printf("  %f %f\n", float(samples2D[2*(d*nspp +s)]), float(samples2D[2*(d*nspp +s)+1]));
                    
                    fprintf(out, "%f %f\n", float(samples2D[2*(d*nspp +s)] + x), float(samples2D[2*(d*nspp +s)+1] + y));
                }
            }
        }
        fclose(out);
    }
    
    // Z permutations
    {
        ZTable table;
        initZTable(table, 4);
        
        FILE *out= fopen("z.dat", "wt");
        assert(out);
        for(int y= 0; y < bounds; y++)
        for(int x= 0; x < bounds; x++)
        {
            std::vector<Float> samples2D(2*nspp*ndims);
            for(int d = 0; d < ndims; ++d) 
            {
                Z2D(table, x, y, 1, samplesPerPixel, log2Resolution, &samples2D[2*d*nspp], d);
                
                printf("pixel (%d, %d) dim %d\n", x, y, d);
                for(int s= 0; s < nspp; s++)
                {
                    printf("  %f %f\n", float(samples2D[2*(d*nspp +s)]), float(samples2D[2*(d*nspp +s)+1]));
                    
                    fprintf(out, "%f %f\n", float(samples2D[2*(d*nspp +s)] + x), float(samples2D[2*(d*nspp +s)+1] + y));
                }
            }
        }
        fclose(out);
    }
    
    return 0;
}
