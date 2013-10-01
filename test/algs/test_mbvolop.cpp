/**
 * \file test_MBVolOp.cpp \test
 *
 *  Started on: Jan 13, 2012
 *  test for volume creation operation
 *  volumes are created with 2 initial surfaces, top and bottom, one cropping line, and
 *  several splitting lines; also, a main direction for all operations, including weaving
 *  top and bottom are first cropped, then split several times (in this example, 2 times)
 */
#include <iostream>
#include <fstream>

#include <time.h>
#include <stdlib.h>
#include <cstring>

#include "meshkit/MKCore.hpp"
#include "meshkit/MeshOp.hpp"
#include "meshkit/ModelEnt.hpp"
#include "meshkit/MBVolOp.hpp"
#include "meshkit/MBSplitOp.hpp"


using namespace MeshKit;

#include "TestUtil.hpp"
#include "meshkit/ReadPolyLine.hpp"

MKCore *mk;

std::string usage_string =
" <input file bottom> <input file top>    input files.\n"
"< direction > 3 doubles for sweeping direction\n"
" <output file >   output moab database \n";

std::string filename_top;
std::string filename_bot;
double direction[3]={0., 0., 1.};
std::string output_file_name;


int main(int argc, char* argv[])
{
  // check command line arg

  filename_bot = TestDir + "/BedCropL3.h5m";
  filename_top = TestDir + "/SurfCropL3.h5m";//"/polyPB.txt";

  output_file_name = "volumesIce.h5m"; // output

  //number_tests_successful = 0;

  if (argc<=1)
  {
    std::cout<<usage_string;
    std::cout << "\n\n";
    std::cout<< "default arguments:" << filename_bot << " "<<filename_top << " "
        << " "<<direction[0] << " " << direction[1] << " "<< direction[2]<<" " << output_file_name <<"\n";
  }
  else if (argc==7)
  {
    filename_bot = argv[1];
    filename_top = argv[2];
    direction[0] = atof(argv[3]);
    direction[1] = atof(argv[4]);
    direction[2] = atof(argv[5]);
    output_file_name= argv[6];
  }
  else
  {
    std::cout << usage_string << " abort.\n";
    return 1;
  }
  // initialize everything

  mk = new MKCore();
  mk->load_mesh(filename_bot.c_str());
  MEntVector  botFaces;
  // we should have 2 faces only, for top and bottom
  mk->get_entities_by_dimension(2, botFaces);
  // load a second one
  mk->load_mesh(filename_top.c_str());
  MEntVector  allFaces;
  mk->get_entities_by_dimension(2, allFaces);

  MBVolOp *eVolOp = (MBVolOp*) mk->construct_meshop("MBVolOp", allFaces);

  //selection.push_back(*dum.rbegin());// push just the last one retrieved from core
  eVolOp->set_direction(direction[0], direction[1], direction[2]);

  std::cout<<"Total number of faces: " << allFaces.size() << "\n";
  mk->setup_and_execute();

  if(output_file_name.length()>0)
  {
    std::cout << "writing output to " << output_file_name << std::endl;
    // output only the root set with the volume
    //
    MEntSelection & mentset = eVolOp->me_selection();
    moab::Range rangeResult = (*(mentset.begin()) ).second;
    moab::EntityHandle resultSet = rangeResult[0]; // the first set
    mk->moab_instance()->write_mesh(output_file_name.c_str(), &resultSet, 1);// write the result set only
  }

  return 0;
}

