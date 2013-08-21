/*!
\example example_onetooneswept.cpp

\section example_onetooneswept_cpp_title <pretty-name-of-this-file>

\subsection example_onetooneswept_cpp_in Input
\image html example_onetooneswept.in.jpg "(description of image)"
There is no input.

\subsection example_onetooneswept_cpp_out Output
\image html example_onetooneswept.out.jpg "(description of image)"

\subsection example_onetooneswept_cpp_inf Misc. Information
\author <your-name-here>
\date 7-15-2013
\bug <placeholder>
\warning <placeholder>

\subsection example_onetooneswept_cpp_src Source Code
*/

#include "meshkit/MKCore.hpp"
#include "meshkit/OneToOneSwept.hpp"
#include "meshkit/SizingFunction.hpp"
#include "meshkit/ModelEnt.hpp"
#include "meshkit/Matrix.hpp"
#include "meshkit/VertexMesher.hpp"


using namespace MeshKit;



std::string gfile_name;
std::string mfile_name;

// these are indices in the volume for source and target surfaces
int index_source, index_target;
int mesh_intervals;

MKCore *mk = NULL;

void test_brick();


int main(int argc, char **argv) 
{
  if (argc==1)
  {
#if HAVE_OCC
    gfile_name = TestDir + "/brick.stp";
    mfile_name = TestDir + "/sf.h5m";
    index_source = 0; index_target = 1; mesh_intervals = 6;
#else
    gfile_name = TestDir + "/BrickWithSrcMeshed.cub";
    mfile_name = gfile_name;
    index_source = 1; index_target = 0; mesh_intervals = 6;
#endif
    std::cout<<"using default arguments: ";
  }
  else if (argc==6)
  {
    gfile_name = argv[1];
    mfile_name = argv[2];
    index_source = atoi(argv[3]);
    index_target = atoi(argv[4]);
    mesh_intervals = atoi(argv[5]);
    std::cout << "Using arguments: ";
  }
  else
  {
    std::cout<<"Usage: " << argv[0] << " <geo_file> <mesh_file> <index_src> <index_tar> <mesh_intervals>\n";
    return 1; // error
  }
  std::cout << argv[0] << " " << gfile_name << " "
         << mfile_name << " " << index_source << " " << index_target << " " << mesh_intervals <<
         "\n";
  // start up MK and load the geometry
  mk = new MKCore();

  int num_fail = 0;
  
  test_brick();

  return num_fail;
  
}

void test_brick()
{

  mk->load_geometry_mesh(gfile_name.c_str(), mfile_name.c_str());

	// get the volumes
	MEntVector vols, surfs, curves, vertices;
	mk->get_entities_by_dimension(3, vols);
	std::cout << "Volume size = " << vols.size() << std::endl;	

	ModelEnt *this_vol = (*vols.rbegin());

	//make a one-to-one sweeping
	OneToOneSwept *sw = (OneToOneSwept*) mk->construct_meshop("OneToOneSwept", vols);

	sw->SetSourceSurface(index_source);
	sw->SetTargetSurface(index_target);

	//set up the size
	SizingFunction swSize(mk, mesh_intervals, -1);
	this_vol->sizing_function_index(swSize.core_index());

	//set up for the sweeping and sweep
	mk->setup_and_execute();

	//check the number of cells after OneToOneSwept
	moab::Range hex;
  moab::ErrorCode rval = mk->moab_instance()->get_entities_by_dimension(0, 3, hex);
  CHECK_EQUAL(moab::MB_SUCCESS, rval);
  std::cout << " generated " << hex.size() << " hexahedrons\n";

	mk->save_mesh("OneToOneSwept.h5m");

	delete sw;
	delete mk->vertex_mesher();
	mk->clear_graph();
		
}

