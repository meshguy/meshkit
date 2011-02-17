#include "meshkit/CAMALTetMesher.hpp"
#include "meshkit/MeshScheme.hpp"
#include "meshkit/ModelEnt.hpp"
#include "moab/Interface.hpp"
#include "moab/ReadUtilIface.hpp"
#include "moab/Range.hpp"
#include "meshkit/RegisterMeshOp.hpp"
#include "CAMALSurfEval.hpp"
#include "CAMALSizeEval.hpp"
#include "CMLTetMesher.hpp"
#include <vector>

namespace MeshKit
{

moab::EntityType CAMALTetMesher::meshTps[] = {moab::MBVERTEX, moab::MBTET, moab::MBMAXTYPE};

CAMALTetMesher::CAMALTetMesher(MKCore *mk_core, const MEntVector &me_vec)
        : MeshScheme(mk_core, me_vec)
{
}

CAMALTetMesher::~CAMALTetMesher()
{
}

MeshOp *CAMALTetMesher::get_tri_mesher() 
{
  std::vector<MeshOpProxy *> proxies;
  mk_core()->meshop_by_mesh_type(moab::MBTRI, proxies);
  if (proxies.empty()) throw Error(MK_FAILURE, "Couldn't find a MeshOp capable of producing triangles.");
  return mk_core()->construct_meshop(*proxies.begin());
}

void CAMALTetMesher::setup_this()
{
  MeshOp *tri_mesher = NULL;
  std::vector<MeshOp*> meshops;
  MEntVector surfs;
  
  for (MEntSelection::iterator sit = mentSelection.begin(); sit != mentSelection.end(); sit++) {
    ModelEnt *me = (*sit).first;
    if (me->dimension() != 3) throw Error(MK_BAD_INPUT, "Tet mesher assigned to an entity with dimension != 3.");
    
      // get the bounding faces
    surfs.clear();
    me->get_adjacencies(2, surfs);
    
      // check the mesh scheme on each one; if there is one, verify it can generate tris; if there
      // isn't one, make one
    bool inserted = false;
    
    for (MEntVector::iterator vit = surfs.begin(); vit != surfs.end(); vit++) {
      meshops.clear();
      (*vit)->get_meshops(meshops);
      if (meshops.empty()) {
          // get a tri mesher if we haven't already
        if (!tri_mesher) tri_mesher = get_tri_mesher();
        
          // add this surface to it, and if first for the volume, make sure it's added upstream
        tri_mesher->add_modelent(*vit);
        if (!inserted) {
          mk_core()->insert_node(tri_mesher, this);
          inserted = true;
        }
      } // if no meshops
    } // over surfs
  } // over vols
}

void CAMALTetMesher::execute_this()
{
  for (MEntSelection::iterator sit = mentSelection.begin(); sit != mentSelection.end(); sit++) {
      // make a me, for convenience
    ModelEnt *me = (*sit).first;

    // assemble bounding mesh
    std::vector<moab::EntityHandle> bdy;
    std::vector<int> group_sizes, senses, bdy_ids;
    me->boundary(2, bdy, &senses, &group_sizes);
    
      // get connectivity 

      // convert the handles to integers for input to CMLTetMesher
    moab::Range bdy_vrange;
    std::vector<double> coords;
    me->get_indexed_connect_coords(bdy, &senses, NULL, bdy_ids, coords, &bdy_vrange);

    CMLTetMesher tet_mesher;
    bool success = tet_mesher.set_boundary_mesh(bdy_vrange.size(), &coords[0], bdy.size(), &bdy_ids[0]);
    if (!success)
      ECERRCHK(MK_FAILURE, "Failed setting boundary mesh");

      // generate the mesh
    int num_tets, num_pts;
    success = tet_mesher.generate_mesh(num_pts, num_tets);
    if (!success)
      ECERRCHK(MK_FAILURE, "Trouble generating tet mesh.");

    moab::Range &new_ents = (*sit).second;
    moab::ErrorCode rval;

      // resize the coords array, then get the coords of the new points
    assert(num_pts >= (int)bdy_vrange.size());
    if (num_pts > (int)bdy_vrange.size()) {
      coords.resize(3*(num_pts-bdy_vrange.size()));
      int pts_returned = tet_mesher.get_points_buf(coords.size(), &coords[0], bdy_vrange.size());
      if (pts_returned != num_pts-bdy_vrange.size()) 
        ECERRCHK(MK_FAILURE, "Number of new points returned from TetMesher doesn't agree with previous value output.");
    
        // create the new vertices' entities 
      rval = mk_core()->moab_instance()->create_vertices(&coords[0], pts_returned, new_ents);
      MBERRCHK(rval, mk_core()->moab_instance());
    }

      // for tets, pre-allocate connectivity
    moab::ReadUtilIface *iface;
    rval = mk_core()-> moab_instance() -> query_interface("ReadUtilIface", (void**)&iface);
    MBERRCHK(rval, mk_core()->moab_instance());		

      //create the tris, get a direct ptr to connectivity
    moab::EntityHandle starth, *connect, *tmp_connect;
    rval = iface->get_element_connect(num_tets, 4, moab::MBTET, 1, starth, connect);
    MBERRCHK(rval, mk_core()->moab_instance());

      // read connectivity directly into that array, as int's
    int *connecti = (int*) connect;
    int tets_returned = tet_mesher.get_tets_buf(4*num_tets, connecti);
    if (tets_returned != num_tets)
      ECERRCHK(MK_FAILURE, "Number of new tets returned from TetMesher doesn't agree with previous value output.");

      // put vertex handles into an indexible array
    std::vector<moab::EntityHandle> bdy_vvec;
    std::copy(bdy_vrange.begin(), bdy_vrange.end(), std::back_inserter(bdy_vvec));
    std::copy(new_ents.begin(), new_ents.end(), std::back_inserter(bdy_vvec));
    
      // now convert vertex indices into handles in-place, working from the back
    for (int i = 4*num_tets-1; i >= 0; i--) {
      assert(connecti[i] >= 0 && connecti[i] < (int)bdy_vvec.size());
      connect[i] = bdy_vvec[connecti[i]];
    }
    
      // put new tris into new entity range, then commit the mesh
    new_ents.insert(starth, starth+num_tets-1);
    me->commit_mesh(new_ents, COMPLETE_MESH);
  }
}

} // namespace MeshKit
