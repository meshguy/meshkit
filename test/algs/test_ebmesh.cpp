/** \file test_ebmesh.cpp
 *
 * Test the EBMesher
 *
 */

#include "meshkit/MKCore.hpp"
#include "meshkit/MeshOp.hpp"
#include "meshkit/EBMesher.hpp"
#include "meshkit/SCDMesh.hpp"
#include "meshkit/ModelEnt.hpp"

using namespace MeshKit;

#include "TestUtil.hpp"

#ifdef HAVE_ACIS
#define DEFAULT_TEST_FILE "sphere.sat"
#elif defined(HAVE_OCC)
#define DEFAULT_TEST_FILE "sphere.stp"
#endif

const bool debug_ebmesher = true;

int load_and_mesh(const char *input_filename,
		  const char *output_filename,
		  int* n_intervals, double box_increase,
		  int vol_frac_res);

int main(int argc, char **argv) 
{
#ifdef HAVE_OCC
  // check command line arg
  std::string input_filename;
  const char *output_filename = NULL;
  int n_interval[3] = {10, 10, 10};
  double box_increase = .2;
  int vol_frac_res = 0;

  if (argc > 2 && argc < 9) {
    input_filename += argv[1];
    if (argc > 2) n_interval[0] = atoi(argv[2]);
    if (argc > 3) n_interval[1] = atoi(argv[3]);
    if (argc > 4) n_interval[2] = atoi(argv[4]);
    if (argc > 5) output_filename = argv[5];
    if (argc > 6) box_increase = atof(argv[6]);
    if (argc > 7) vol_frac_res = atoi(argv[7]);
  }
  else {
    std::cout << "Usage: " << argv[0] << "<input_geom_filename> {x: # of intervals} {y: # of intervals} {z: # of intervals} {output_mesh_filename} {#_add_layer} {vol_frac_res}" << std::endl;
    std::cout << "<input_geom_filename> : input geometry file name" << std::endl;
    std::cout << "{x/y/z: # ofintervals} : optional argument. # of intervals. if it is not set, set to 10." << std::endl;
    std::cout << "{output_mesh_filename} : optional argument. if it is not set, dosn't export. can output mesh file (e.g. output.vtk.)" << std::endl;
    std::cout << "{box size increase} : optional argument. Cartesian mesh box increase form geometry. default 0.03" << std::endl;
    std::cout << "{vol_frac_res} : optional argument, volume fraction resolution of boundary cells for each material, you can specify it as # of divisions (e.g. 4)." << std::endl;
    std::cout << std::endl;
    if (argc != 1) return 1;
    std::cout << "No file specified.  Defaulting to: " << DEFAULT_TEST_FILE << std::endl;
    std::string file_name = TestDir + "/" + DEFAULT_TEST_FILE;
    input_filename += TestDir;
    input_filename += "/";
    input_filename += DEFAULT_TEST_FILE;
  }
  
  if (load_and_mesh(input_filename.c_str(), output_filename,
		    n_interval, box_increase, vol_frac_res)) return 1;
  
  return 0;
}

int load_and_mesh(const char *input_filename,
		  const char *output_filename,
		  int* n_interval, double box_increase,
		  int vol_frac_res)
{
  bool result;
  time_t start_time, load_time, mesh_time, vol_frac_time,
    export_time, query_time_techX, query_time;

  // start up MK and load the geometry
  MKCore mk;
  time(&start_time);
  mk.load_mesh(input_filename);
  time(&load_time);

  if (debug_ebmesher) {
    mk.save_mesh("input.vtk");
  }

  // populate mesh to relate geometry entities and mesh sets
  mk.populate_mesh();

  // get the volumes
  MEntVector vols;
  mk.get_entities_by_dimension(3, vols);

  // make and set input for structured mesher
  double box_min[3], box_max[3];
  SCDMesh *sm = (SCDMesh*) mk.construct_meshop("SCDMesh", vols);
  sm->set_name("structured_mesh");
  sm->set_interface_scheme(SCDMesh::full);
  sm->set_grid_scheme(SCDMesh::cfMesh);
  sm->set_axis_scheme(SCDMesh::cartesian);
  sm->set_box_increase_ratio(box_increase); // add some extra layer to box
  sm->set_box_dimension(); // set box dimension
  sm->get_box_dimension(box_min, box_max); // get box dimension

  // set # of intervals for 3 directions
  std::vector<int> fine_i (n_interval[0], 1);
  sm->set_coarse_i_grid(n_interval[0]);
  sm->set_fine_i_grid(fine_i);
  std::vector<int> fine_j (n_interval[1], 1);
  sm->set_coarse_j_grid(n_interval[1]);
  sm->set_fine_j_grid(fine_j);
  std::vector<int> fine_k (n_interval[2], 1);
  sm->set_coarse_k_grid(n_interval[2]);
  sm->set_fine_k_grid(fine_k);

  // make EBMesher
  EBMesher *ebm = (EBMesher*) mk.construct_meshop("EBMesher", vols);
  ebm->set_name("embedded_boundary_mesh");
  ebm->set_division(box_min, box_max, n_interval);

  // put them in the graph
  mk.get_graph().addArc(mk.root_node()->get_node(), sm->get_node());
  mk.get_graph().addArc(sm->get_node(), ebm->get_node());
  mk.get_graph().addArc(ebm->get_node(), mk.leaf_node()->get_node());

  // mesh embedded boundary mesh, by calling execute
  mk.setup_and_execute();
  time(&mesh_time);

  // caculate volume fraction, only for geometry input
  if (vol_frac_res > 0) {
    result = ebm->get_volume_fraction(vol_frac_res);
    if (!result) {
      std::cerr << "Couldn't get volume fraction." << std::endl;
      return 1;
    }
  }
  time(&vol_frac_time);

  // export mesh
  if (output_filename != NULL) {
    ebm->export_mesh(output_filename);
  }
  time(&export_time);

  // techX query function test
  double boxMin[3], boxMax[3];
  int nDiv[3];
  std::map< CutCellSurfEdgeKey, std::vector<double>, LessThan > mdCutCellSurfEdge;
  std::vector<int> vnInsideCellTechX;

  ebm->get_grid_and_edges_techX(boxMin, boxMax, nDiv,
					 mdCutCellSurfEdge, vnInsideCellTechX);
  time(&query_time_techX);

  // multiple intersection fraction query test
  std::map< CutCellSurfEdgeKey, std::vector<double>, LessThan > mdCutCellEdge;
  std::vector<int> vnInsideCell;
  result = ebm->get_grid_and_edges(boxMin, boxMax, nDiv,
				   mdCutCellEdge, vnInsideCell);
  if (!result) {
    std::cerr << "Couldn't get mesh information." << std::endl;
    return 1;
  }
  time(&query_time);

  std::cout << "EBMesh is succesfully finished." << std::endl;
  std::cout << "# of TechX cut-cell surfaces: " << mdCutCellSurfEdge.size() 
	    << ", # of nInsideCell: " << vnInsideCell.size()/3 << std::endl;
  std::cout << "Time including loading: "
	    << difftime(mesh_time, start_time)
	    << " secs, Time excluding loading: "
	    << difftime(mesh_time, load_time)
	    << " secs, Time volume fraction: "
	    << difftime(vol_frac_time, mesh_time)
	    << " secs, Time export mesh: "
	    << difftime(export_time, vol_frac_time)
	    << " secs, TechX query time: "
	    << difftime(query_time_techX, export_time)
	    << " secs, normal query time(elems, edge-cut fractions): "
	    << difftime(query_time, query_time_techX)
	    << " secs." << std::endl;
#else
  std::cout << "Current EBMesher is not supported for ACIS geometry. Only works with OCC." << std::endl;
#endif
  return 0;
}
