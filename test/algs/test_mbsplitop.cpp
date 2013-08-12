/*! \file test_MBSplitOp.cpp \test
 *
 * test_MBSplitOp.cpp
 *
 */

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
#include "ReadPolyLine.hpp"

MKCore *mk;

std::string usage_string =
" <input file>      input file.\n"
" <polyline file> file with the direction and polyline for cropping\n"
" <closed loop>   1 or 0 for closed or not \n"
" <surface_id>    global id of the surface to split \n"
" <min_dot>       min dot for max angle allowed betwen consecutive edges"
" <output file >   output moab database \n";

std::string filenameS;
std::string filename_out;
std::string polygon_file_name;
int closed = 0;
int surfId = 1;// surface id to split
double min_dot = 0.8;

int main(int argc, char* argv[])
{
  // check command line arg
  const char *filename = 0;
  const char *file_poly = 0;
  const char *outfile = 0;
  closed = 1;// to split or not
  surfId = 1;
  min_dot = 0;


  filenameS = TestDir + "/PB.h5m";
  polygon_file_name = TestDir + "/polyPB.txt";
  filename_out = "PB_new.h5m";

  if (argc<=1)
  {
    std::cout<<usage_string;
    std::cout << "\n\n";
    std::cout<< "default arguments: ../../data/PB.h5m ../../data/polyPB.txt 1 1 0.8 PB_new.h5m \n";
    filename = filenameS.c_str();
    file_poly = polygon_file_name.c_str();
    outfile = 0; // do not output if default, do not save the db
  }
  else if (argc==7)
  {
    filename = argv[1];
    file_poly = argv[2];
    closed = atoi(argv[3]);
    surfId = atoi(argv[4]);
    min_dot = atof(argv[5]);
    outfile = argv[6];
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

  std::vector<double> xyz;
  double direction[3];

  int rc = ReadPolyLineFromFile(file_poly, direction, xyz);
  if (rc !=0)
  {
    std::cout<<" can't read from polyline file\n";
    return rc;

  }
  int sizePolygon = (int)xyz.size()/3;
  if ((closed && sizePolygon < 3 ) || (!closed&& sizePolygon<2)) {
    std::cerr << " Not enough points in the polygon" << std::endl;
    return 1;
  }

  splitOp->set_options( /* int globalId*/ surfId, direction[0], direction[1],
      direction[2], /* int closed*/ closed, /*min_dot*/ min_dot);

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
