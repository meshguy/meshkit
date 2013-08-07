/*!
\example TriangleMesher.cpp

\section TriangleMesher_cpp_title <pretty-name-of-this-file>

\subsection TriangleMesher_cpp_in Input
\image html TriangleMesher.in.jpg
There is no input.

\subsection TriangleMesher_cpp_out Output
\image html TriangleMesher.out.jpg

\subsection TriangleMesher_cpp_inf Misc. Information
\author <your-name-here>
\date 7-15-2013
\bug <placeholder>
\warning <placeholder>

\subsection TriangleMesher_cpp_src Source Code
*/

#include <iostream>

#include <time.h>
#include <stdlib.h>
#include <cstring>
#include <cstring>

#include "meshkit/MKCore.hpp"
#include "meshkit/MeshOp.hpp"
#include "meshkit/ModelEnt.hpp"
#include "meshkit/TriangleMesher.hpp"

using namespace MeshKit;



MKCore *mk;

std::string usage_string =
"-o <file>      Output final model to given file.\n"
"-i <file>      input file \n"
"-O <options>   passed to Triangle mesher \n"
"-d  <1, 2, 3>  (create the mesh in x, y, z normal planes direction, default z(3)) "
"\n";


int main(int argc, char* argv[])
{
  // check command line arg
  const char *filename = 0;
  const char *outfile = 0;

  //process options
  // process all options, as given with regular qslim
  // the last argument must be the input file
  // one argument must be -o <the output file>, otherwise it is the standard output
  // in the process, a mesh file is read (instanced) and the
  // set of triangles is passed
  // it could be the root set for the mesh
  std::string fstr;
  char * opts  = (char *)("pc");
  int direction = 3; // default
  double fretting = 1.e+38;// huge, do not eliminate anything in general
  if (argc<=1)
  {
    std::cout<<usage_string;
    std::cout << "\n\n";
    std::cout<< "default arguments: -i TriangleInput.h5m -o mesh.h5m -O pc  -f 2. \n";
    //opts = "pc";

    fstr=TestDir + "/TriangleInput.h5m";
    filename = fstr.c_str();
    //opts = "pc"; // convex hull, poly input

    // ostr = "out.smf";
    // outfile = ostr.c_str();
  }
  else
  {
    int i=1;// will loop through arguments, and process them
    for (i=1; i<argc ; i++)
    {
      if (argv[i][0]=='-')
      {
        switch (argv[i][1])
        {
          case 'i':
          {
            filename = argv[i+1];
            i+=1;
            break;
          }
          case 'o':
          {
            outfile = argv[i+1];
            i+=1;
            break;
          }

          case 'O':
          {
            opts = argv[i+1];
            i+=1;
            break;
          }

          case 'd':
          {
            direction = atoi(argv[i+1]);
            i+=1;
            break;
          }

          case 'f':
          {
            fretting = atof(argv[i+1]);
            i+=1;
            break;
          }
          default :
          {
           std::cout << "unsupported wrong input argument " << argv[i] <<std::endl;
           return 1;
          }
        }
      }
    }
  }
  // initialize everything

  mk = new MKCore();
  mk->load_mesh(filename);
  MEntVector selection, dum;
  mk->get_entities_by_dimension(2, dum);
  selection.push_back(*dum.rbegin());// push just the last one retrieved from core

  TriangleMesher *qm = (TriangleMesher*) mk->construct_meshop("TriangleMesher", selection);

  qm->set_options(opts, direction, fretting);

  mk->setup_and_execute();

  if(outfile)
  {
    std::cout << "writing output to " << outfile << std::endl;
    mk->moab_instance()->write_mesh(outfile);// write everything left
  }

  return 0;
}
