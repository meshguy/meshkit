/*!
\example example_parallelmesher.cpp

\section example_parallelmesher_cpp_title <pretty-name-of-this-file>

\subsection example_parallelmesher_cpp_in Input
\image html example_parallelmesher.in.jpg
There is no input.

\subsection example_parallelmesher_cpp_out Output
\image html example_parallelmesher.out.jpg

\subsection example_parallelmesher_cpp_inf Misc. Information
\author <your-name-here>
\date 7-15-2013
\bug <placeholder>
\warning <placeholder>

\subsection example_parallelmesher_cpp_src Source Code
*/

// Test file for parallel mesh generation
// It reads geometry in parallel with CGM and meshed by camal in parallel.
#include <iostream>
#include <sstream>

#include "meshkit/MKCore.hpp"
#include "meshkit/ParallelMesher.hpp"
#include "meshkit/SizingFunction.hpp"
#include "meshkit/ModelEnt.hpp"

using namespace MeshKit;



#ifdef HAVE_ACIS
#define DEFAULT_TEST_FILE "bricks_3.sat"
#elif defined (HAVE_OCC)
#define DEFAULT_TEST_FILE "bricks_3.occ"
#endif

int load_and_mesh(const char *geom_filename,
                  const char *mesh_filename,
                  const char *options, const double interval_size,
                  const int n_interval, const int rank);

int main( int argc, char *argv[] )
{
  int rank, size;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  std::cout << "Hello" << std::endl;
 
  // Check command line arg
  std::string geom_filename;
  const char *mesh_filename = NULL;
  int send_method = 1; // read and delete
  int part_method = 0; // round_robin
  double mesh_size = 0.1;
  int mesh_interval = -1;
  bool force_intervals = false;
  std::string options;

  if (argc < 2) {
    if (rank == 0) {
      std::cout << "Usage: " << argv[0] << " <input_geom_filename> <output_mesh_filename> [<send_methond>] [<partition_method>] [-s <uniform_size>] [-i <uniform_int>] [-f] "
                << std::endl
                << "  <send_method> = 0:read, 1:read_delete, 2:broadcast"
                << " , 3:bcast_delete, 4:scatter,"
                << " , 5:read_parallel" << std::endl;
      std::cout << " <partition_method> = 0:round-robin, 1:static" << std::endl;
      std::cout << "  -s <uniform_size> = mesh with this size" << std::endl;
      std::cout << "  -i <uniform_int> = mesh curves with this # intervals" << std::endl;
      std::cout << "  -f = force these size/interval settings even if geometry has interval settings" << std::endl;
      std::cout << "No file specified.  Defaulting to: " << DEFAULT_TEST_FILE << std::endl;
      std::cout << "No send method specified.  Defaulting to: read_delete" << std::endl;
    }
    std::string file_name = TestDir + "/" + DEFAULT_TEST_FILE;
    geom_filename += TestDir;
    geom_filename += "/";
    geom_filename += DEFAULT_TEST_FILE;
  }
  else {
    geom_filename = argv[1];
    if (argc > 2) mesh_filename = argv[2];
    if (argc > 3) send_method = atoi(argv[3]);
    if (argc > 4) part_method = atoi(argv[4]);
    if (argc > 5) {
      int argno = 5;
      while (argno < argc) {
        if (!strcmp(argv[argno], "-s")) {
          argno++;
          sscanf(argv[argno], "%lf", &mesh_size);
          argno++;
        }
        else if (!strcmp(argv[argno], "-i")) {
          argno++;
          sscanf(argv[argno], "%d", &mesh_interval);
          argno++;
        }
        else if (!strcmp(argv[argno], "-f")) {
          argno++;
          force_intervals = true;
        }
        else {
          std::cerr << "Unrecognized option: " << argv[argno] << std::endl;
          return 1;
        }
      }
    }
#ifdef HAVE_ACIS
    if (send_method != 1) {
      std::cerr << "Other send methods than read_and_delete for ACIS geometry are not supported." << std::endl;
      return 0;
    }
#endif
  }

  // set geometry send method
  if (send_method == 0) options += "PARALLEL=READ;";
  else if (send_method == 1) options += "PARALLEL=READ_DELETE;";
  else if (send_method == 2) options += "PARALLEL=BCAST;";
  else if (send_method == 3) options += "PARALLEL=BCAST_DELETE;";
  else if (send_method == 4) options += "PARALLEL=SCATTER;";
  else if (send_method == 5) options += "PARALLEL=READ_PARALLEL;";
  else {
    std::cout << "Send method " << send_method
              << " is not supported. Defaulting to: read_delete" << std::endl;
    options = "PARALLEL=READ_DELETE;";
  }

  // set partition method
  if (part_method == 0) options += "PARTITION=GEOM_DIMENSION;PARTITION_VAL=3;PARTITION_DISTRIBUTE;";
  else if (part_method == 1) options += "PARTITION=PAR_PARTITION_STATIC;";
  else {
    std::cout << "Partition method " << part_method
              << " is not supported. Defaulting to: round_robin" << std::endl;
    options = "PARTITION=GEOM_DIMENSION;PARTITION_VAL=3;PARTITION_DISTRIBUTE;";
  }

  if (load_and_mesh(geom_filename.c_str(), mesh_filename,
                    options.c_str(), mesh_size, mesh_interval, rank)) return 1;

  return 0;
}

int load_and_mesh(const char *geom_filename,
                  const char *mesh_filename,
                  const char *options, const double interval_size,
                  const int n_interval, const int rank)
{
  // start up MK and load the geometry
  double t1 = MPI_Wtime();
  MKCore *mk = new MKCore;
  mk->load_geometry(geom_filename, options);
  MPI_Barrier(MPI_COMM_WORLD);
  double t2 = MPI_Wtime();

  // get the volumes
  MEntVector dum, vols, part_vols;
  mk->get_entities_by_dimension(3, vols);

  // make a sizing function and set it on the surface
  SizingFunction esize(mk, n_interval, interval_size);
  unsigned int i_sf = esize.core_index();
  for (int i = 0; i < vols.size(); i++) vols[i]->sizing_function_index(i_sf);
  
  // do parallel mesh
  mk->construct_meshop("ParallelMesher", vols);
  mk->setup_and_execute();
  double t3 = MPI_Wtime();

  if (mesh_filename != NULL) {
    std::string out_name;
    std::stringstream ss;
    ss << mesh_filename << ".h5m";
    ss >> out_name;
    iMesh::Error err = mk->imesh_instance()->save(0, out_name.c_str(),
                                                  "moab:PARALLEL=WRITE_PART");
    IBERRCHK(err, "Couldn't save mesh.");
    double t4 = MPI_Wtime();
    std::cout << "Export_time="
              << t4 - t3 << std::endl;
  }

  MPI_Barrier(MPI_COMM_WORLD);
  double t5 = MPI_Wtime();

  // report the number of tets
  moab::Range tets;
  moab::ErrorCode rval = mk->moab_instance()->get_entities_by_dimension(0, 3, tets);
  CHECK_EQUAL(moab::MB_SUCCESS, rval);
  std::cout << tets.size() << " tets generated." << std::endl;

  if (rank == 0) {
    std::cout << "Geometry_loading_time=" << t2 - t1
              << " secs, Meshing_time=" << t3 - t2
              << " secs, Total_time=" << t5 - t1 << std::endl;
  }
  
  mk->clear_graph();
  MPI_Finalize();
  
  return 0;
}
