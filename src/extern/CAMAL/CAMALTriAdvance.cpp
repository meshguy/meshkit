#include "meshkit/CAMALTriAdvance.hpp"
#include "meshkit/MeshScheme.hpp"
#include "meshkit/ModelEnt.hpp"
#include "moab/Interface.hpp"
#include "moab/ReadUtilIface.hpp"
#include "moab/Range.hpp"
#include "meshkit/RegisterMeshOp.hpp"
#include "CAMALSurfEval.hpp"
#include "CAMALSizeEval.hpp"
#include "CMLTriAdvance.hpp"
#include "RefEntity.hpp"

#ifdef PARALLEL
#ifdef HAVE_PARALLEL_MOAB
#include "moab/ParallelComm.hpp"
#endif
#endif

#include <vector>

const bool debug_camaltriadv = false;

namespace MeshKit
{

moab::EntityType CAMALTriAdvance::meshTps[] = {moab::MBVERTEX, moab::MBTRI, moab::MBMAXTYPE};

CAMALTriAdvance::CAMALTriAdvance(MKCore *mk_core, const MEntVector &me_vec)
        : MeshScheme(mk_core, me_vec)
{
}

CAMALTriAdvance::~CAMALTriAdvance()
{
}

void CAMALTriAdvance::setup_this()
{
    // just call setup_boundary, since that's the only constraint we have
  setup_boundary();
}

void CAMALTriAdvance::execute_this()
{
  for (MEntSelection::iterator sit = mentSelection.begin(); sit != mentSelection.end(); sit++) {
      // make a me, for convenience
    ModelEnt *me = (*sit).first;

      // create a surface evaluator for this modelent, and a size evaluator
    CAMALSurfEval cse(me);

    SizingFunction * sz = mk_core()->sizing_function(me->sizing_function_index());
    if (!sz->variable())
      sz = NULL; // no function variable

    CAMALSizeEval mesize(me->mesh_interval_size(), sz);
      // make sure the size isn't negative
    if (mesize.get_size() == -1) mesize.set_size(1.0);
    
    // assemble bounding mesh
    std::vector<moab::EntityHandle> bdy;
    std::vector<int> group_sizes, bdy_ids;
    me->boundary(0, bdy, NULL, &group_sizes);

      // convert the handles to integers for input to TriAdv
    moab::Range bdy_vrange;
    std::vector<double> coords;
    me->get_indexed_connect_coords(bdy, NULL, NULL, bdy_ids, coords, &bdy_vrange);

    // now construct the CAMAL mesher, and pass it initial conditions
    CMLTriAdvance triadv(&cse, &mesize);
    bool success = triadv.set_boundary_mesh(bdy_vrange.size(), &coords[0], group_sizes.size(), &group_sizes[0], &bdy_ids[0]);
    if (!success)
      ECERRCHK(MK_FAILURE, "Trouble setting boundary mesh.");

      // ok, now generate the mesh
    int num_pts, num_tris;
    success = triadv.generate_mesh(num_pts, num_tris);
    if (!success) {
      if (debug_camaltriadv) print_debug(me, coords, bdy_vrange, group_sizes, bdy_ids);
      ECERRCHK(MK_FAILURE, "Trouble generating tri mesh.");
    }

    moab::Range &new_ents = (*sit).second;
    moab::ErrorCode rval;

      // resize the coords array, then get the coords of the new points
    assert(num_pts >= (int)bdy_vrange.size());
    if (num_pts > (int)bdy_vrange.size()) {
      coords.resize(3*(num_pts-bdy_vrange.size()));
      unsigned int pts_returned = triadv.get_points_buf(coords.size(), &coords[0], bdy_vrange.size());
      if (pts_returned != num_pts-bdy_vrange.size()) 
        ECERRCHK(MK_FAILURE, "Number of new points returned from TriAdv doesn't agree with previous value output.");
    
        // create the new vertices' entities on the face
      rval = mk_core()->moab_instance()->create_vertices(&coords[0], pts_returned, new_ents);
      MBERRCHK(rval, mk_core()->moab_instance());
    }

      // for tris, pre-allocate connectivity
    moab::ReadUtilIface *iface;
    rval = mk_core()-> moab_instance() -> query_interface("ReadUtilIface", (void**)&iface);
    MBERRCHK(rval, mk_core()->moab_instance());		

      //create the tris, get a direct ptr to connectivity
    moab::EntityHandle starth, *connect;
    rval = iface->get_element_connect(num_tris, 3, moab::MBTRI, 1, starth, connect);
    MBERRCHK(rval, mk_core()->moab_instance());

      // read connectivity directly into that array, as int's
    int *connecti = (int*) connect;
    int tris_returned = triadv.get_tris_buf(3*num_tris, connecti);
    if (tris_returned != num_tris)
      ECERRCHK(MK_FAILURE, "Number of new tris returned from TriAdv doesn't agree with previous value output.");

      // put vertex handles into an indexible array
    std::vector<moab::EntityHandle> bdy_vvec;
    std::copy(bdy_vrange.begin(), bdy_vrange.end(), std::back_inserter(bdy_vvec));
    std::copy(new_ents.begin(), new_ents.end(), std::back_inserter(bdy_vvec));
    
      // now convert vertex indices into handles in-place, working from the back
    for (int i = 3*num_tris-1; i >= 0; i--) {
      assert(connecti[i] >= 0 && connecti[i] < (int)bdy_vvec.size());
      connect[i] = bdy_vvec[connecti[i]];
    }
    
      // put new tris into new entity range, then commit the mesh
    new_ents.insert(starth, starth+num_tris-1);
    me->commit_mesh(new_ents, COMPLETE_MESH);
  }
}

void CAMALTriAdvance::print_debug(ModelEnt *me, std::vector<double> &coords,
                                  moab::Range &bdy_vrange, std::vector<int> &group_sizes,
                                  std::vector<int> &bdy_ids)
{
  std::cout << "Surface_bounadry_mesh: mesh_size = "
            << me->mesh_interval_size() << std::endl;
  
  for (int i = 0; i < (int) bdy_vrange.size(); i++) {
    std::cout << coords[3 * i] << "  " << coords[3 * i + 1]
              << "  " << coords[3 * i + 2] << std::endl;
  }
  
  std::cout << "bdy_vertex_size:" << bdy_vrange.size()
            << ", group_size:" << group_sizes.size() << std::endl;
  
  int index = 0;
  for (int i = 0; i < (int) group_sizes.size(); i++) {
    int g_size = group_sizes[i];
    std::cout << "boundary_order_group" << i + 1 << ", group_size="
              << g_size << std::endl;
    for (int j = 0; j < g_size; j++) {
          std::cout << bdy_ids[index + j] << " ";
    }
    std::cout << std::endl;
    index += g_size;
  }
  
  moab::ErrorCode rval;
  moab::EntityHandle outset;
  std::string outfile;
  std::stringstream ss;
  
  RefEntity* entity = reinterpret_cast<RefEntity*> (me->geom_handle());
  ss << "CAMALTri_boundary_surf";
  ss << entity->id();
  ss << "_";
#ifdef PARALLEL
#ifdef HAVE_PARALLEL_MOAB
  moab::ParallelComm* pcomm = moab::ParallelComm::get_pcomm(mk_core()->moab_instance(), 0);
  ss << "proc";
  ss << pcomm->proc_config().proc_rank();
#endif
#endif
  ss >> outfile;
  outfile += ".vtk";
  rval = mk_core()->moab_instance()->create_meshset(0, outset);
  MBERRCHK(rval, mk_core()->moab_instance());
  rval = mk_core()->moab_instance()->add_entities(outset, bdy_vrange);
  MBERRCHK(rval, mk_core()->moab_instance());
  rval = mk_core()->moab_instance()->write_mesh(outfile.c_str(), &outset, 1);
  MBERRCHK(rval, mk_core()->moab_instance());
  std::cout << outfile.c_str() << " is saved." << std::endl;
}
} // namespace MeshKit
