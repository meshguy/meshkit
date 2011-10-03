#include <iostream>
#include <fstream>

#include <time.h>
#include <stdlib.h>
#include <cstring>

#include "meshkit/MKCore.hpp"
#include "meshkit/MeshOp.hpp"
#include "meshkit/ModelEnt.hpp"
#include "meshkit/MBSplitOp.hpp"


using namespace MeshKit;

#include "TestUtil.hpp"

MKCore *mk;

std::string usage_string =
" <input file>      input file.\n"
" <polyline file> file with the direction and polyline for cropping\n"
" <output file >   output moab database \n";

std::string filenameS;
std::string filename_out;
std::string polygon_file_name;


int main(int argc, char* argv[])
{
  // check command line arg
  const char *filename = 0;
  const char *file_poly = 0;
  const char *outfile = 0;

  filenameS = TestDir + "/PB.h5m";
  polygon_file_name = TestDir + "/polyPB.txt";
  filename_out = "PB_new.h5m";

  if (argc<=1)
  {
    std::cout<<usage_string;
    std::cout << "\n\n";
    std::cout<< "default arguments: ../../data/PB.h5m ../../data/polyPB.txt PB_new.h5m \n";
    filename = filenameS.c_str();
    file_poly = polygon_file_name.c_str();
    outfile = 0; // do not output if default, do not save the db
  }
  else if (argc==4)
  {
    filename = argv[1];
    file_poly = argv[2];
    outfile = argv[3];
  }
  else
  {
    std::cout << usage_string << " abort.\n";
    return 1;
  }
  // initialize everything

  mk = new MKCore();
  mk->load_mesh(filename);
  MEntVector  selection;
  mk->get_entities_by_dimension(2, selection);
  //selection.push_back(*dum.rbegin());// push just the last one retrieved from core

  MBSplitOp *splitOp = (MBSplitOp*) mk->construct_meshop("MBSplitOp", selection);

  // get the direction, and the polygon/ polyline points
  std::ifstream datafile(file_poly, std::ifstream::in);
  if (!datafile) {
    std::cout << "can't read polyline file\n";
    return 1;
  }
    //
  char temp[100];
  double direction[3];// normalized
  double gridSize;
  datafile.getline(temp, 100);// first line

  // get direction and mesh size along polygon segments, from file
  sscanf(temp, " %lf %lf %lf %lf ", direction, direction + 1, direction + 2,
      &gridSize);


  std::vector<double> xyz;
  while (!datafile.eof()) {
    datafile.getline(temp, 100);
    //int id = 0;
    double x, y, z;
    int nr = sscanf(temp, "%lf %lf %lf", &x, &y, &z);
    if (nr == 3) {
      xyz.push_back(x);
      xyz.push_back(y);
      xyz.push_back(z);
    }
  }
  int sizePolygon = (int)xyz.size()/3;
  if (sizePolygon < 3) {
    std::cerr << " Not enough points in the polygon" << std::endl;
    return 1;
  }



  splitOp->set_options( /* int globalId*/ 1, direction[0], direction[1],
      direction[2], /* int closed*/ 1);

  for (int k=0 ; k<sizePolygon; k++)
    splitOp->add_points(xyz[3*k], xyz[3*k+1], xyz[3*k+2]);

  mk->setup_and_execute();

  if(outfile)
  {
    std::cout << "writing output to " << outfile << std::endl;
    mk->moab_instance()->write_mesh(outfile);// write everything left
  }

  return 0;
}
