// Test file for parallel mesh generation
// It reads geometry in parallel with CGM and meshed by camal in parallel.
#include <iostream>
#include <sstream>

#include "meshkit/MKCore.hpp"
#include "meshkit/ParallelMesher.hpp"
#include "meshkit/SizingFunction.hpp"
#include "meshkit/ModelEnt.hpp"

using namespace MeshKit;

#include "TestUtil.hpp"

#define DEFAULT_TEST_FILE "bricks_3.occ"

int load_and_mesh(const char *geom_filename,
                  const char *mesh_filename,
                  const char *options, const double interval_size,
                  const int n_interval, const int rank);

int main( int argc, char *argv[] )
{
#ifdef HAVE_OCC
  int rank, size;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  std::cout << "Hello" << std::endl;
 
  // Check command line arg
  std::string geom_filename;
  const char *mesh_filename = NULL;
  int send_method = 3; // broadcst and delete
  double mesh_size = 0.1;
  int mesh_interval = -1;
  bool force_intervals = false;
  std::string options;

  if (argc < 3) {
    if (rank == 0) {
      std::cout << "Usage: " << argv[0] << " <input_geom_filename> <output_mesh_filename> [<send_methond>] [-s <uniform_size>] [-i <uniform_int>] [-f] "
                << std::endl
                << "  <send_method> = 0:read, 1:read_delete, 2:broadcast"
                << " , 3:bcast_delete, 4:scatter,"
                << " , 5:read_parallel" << std::endl;
      std::cout << "  -s <uniform_size> = mesh with this size" << std::endl;
      std::cout << "  -i <uniform_int> = mesh curves with this # intervals" << std::endl;
      std::cout << "  -f = force these size/interval settings even if geometry has interval settings" << std::endl;
      std::cout << "No file specified.  Defaulting to: merge_test.OCC" << std::endl;
      std::cout << "No send method specified.  Defaulting to: bcast_delete" << std::endl;
    }

    std::string file_name = TestDir + "/" + DEFAULT_TEST_FILE;
    geom_filename += TestDir;
    geom_filename += "/";
    geom_filename += DEFAULT_TEST_FILE;
  }
  else {
    geom_filename = argv[1];
    mesh_filename = argv[2];
    if (argc > 3) send_method = atoi(argv[3]);
    if (argc > 4) {
      int argno = 4;
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
              << " is not supported. Defaulting to: broadcast_delete" << std::endl;
    options = "PARALLEL=BCAST_DELETE;";
  }

  // do body partitioning with round robin distribution
  options += "PARTITION=GEOM_DIMENSION;PARTITION_VAL=3;PARTITION_DISTRIBUTE;";

  if (load_and_mesh(geom_filename.c_str(), mesh_filename,
                    options.c_str(), mesh_size, mesh_interval, rank)) return 1;
#else
  std::cout << "Parallel meshing is only supported OCC geometry now." << std::endl;
#endif
  return 0;
}

int load_and_mesh(const char *geom_filename,
                  const char *mesh_filename,
                  const char *options, const double interval_size,
                  const int n_interval, const int rank)
{
  // start up MK and load the geometry
  double t_start = MPI_Wtime();
  MKCore *mk = new MKCore;
  mk->load_geometry(geom_filename, options);
  MPI_Barrier(MPI_COMM_WORLD);
  double t_load = MPI_Wtime();
  std::cout << "Geometry is loaded in "
            << (double) (t_load - t_start)/CLOCKS_PER_SEC
            << " seconds." << std::endl;

  // get the volumes
  MEntVector dum, vols, part_vols;
  mk->get_entities_by_dimension(3, vols);

  // make a sizing function and set it on the surface
  SizingFunction esize(mk, n_interval, interval_size);
  vols[0]->sizing_function_index(esize.core_index());
  
  // do parallel mesh
  mk->construct_meshop("ParallelMesher", vols);
  mk->setup_and_execute();

  // report the number of tets
  moab::Range tets;
  moab::ErrorCode rval = mk->moab_instance()->get_entities_by_dimension(0, 3, tets);
  CHECK_EQUAL(moab::MB_SUCCESS, rval);
  std::cout << tets.size() << " tets generated." << std::endl;

  if (mesh_filename != NULL) {
    std::string out_name;
    std::stringstream ss;
    ss << mesh_filename << rank << ".vtk";
    ss >> out_name;
    iMesh::Error err = mk->imesh_instance()->save(0, out_name.c_str());
    IBERRCHK(err, "Couldn't save mesh.");
  }

  mk->clear_graph();
  MPI_Finalize();
  
  return 0;
}
