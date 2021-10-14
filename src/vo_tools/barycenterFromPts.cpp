#include <iostream>
#include <fstream>
#include <string>
#include "CLI11.hpp"

int main(int argc, char **argv)
{
  std::string filename;
  
  CLI::App app { "barycenterFromPts" };
  app.add_option("-i,--input", filename, "Filename of the (ascii) point set")->required()->check(CLI::ExistingFile);
  unsigned int scale = 1;
  app.add_option("-s,--scale", scale, "Scaling factor to rescale according to the number of points");
  CLI11_PARSE(app, argc, argv)

  std::ifstream in(filename, std::ifstream::in);
  
  double x,y,bx,by;
  unsigned int cpt=0;
  
  while (in.good())
  {
    in >> x;
    in >> y;
    if ( cpt==0 )
    {
      bx=x;
      by=y;
    }
    else
    {
      bx += x;
      by += y;
    }
  if (in.good())  cpt++;
  }

  bx /= (double)cpt;
  by /= (double)cpt;
  bx = (bx-0.5)*scale+0.5;
  by = (by-0.5)*scale+0.5;
  
  std::cout<<cpt<<" "<<bx<<" "<<by<<" "<< (bx-0.5)*(bx-0.5)+(by-0.5)*(by-0.5) <<std::endl;
  
  return 0;
}
  
