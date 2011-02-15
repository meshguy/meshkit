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
#include <vector>

namespace MeshKit
{

moab::EntityType CAMALTriAdvance::meshTps[] = {moab::MBVERTEX, moab::MBTRI};
iBase_EntityType CAMALTriAdvance::geomTps[] = {iBase_FACE};

//---------------------------------------------------------------------------//
    static RegisterMeshOp<CAMALTriAdvance,true> INIT("CAMALTriAdvance", CAMALTriAdvance::geomTps, 2,
                                                     CAMALTriAdvance::meshTps, 2);
//---------------------------------------------------------------------------//

CAMALTriAdvance::CAMALTriAdvance(MKCore *mk_core, const MEntVector &me_vec)
        : MeshScheme(mk_core, me_vec)
{
}

CAMALTriAdvance::~CAMALTriAdvance()
{
}

void CAMALTriAdvance::mesh_types(std::vector<moab::EntityType> &tps)
{
  tps.push_back(moab::MBVERTEX);
  tps.push_back(moab::MBTRI);
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
    CAMALSizeEval mesize(me->mesh_interval_size());
      // make sure the size isn't negative
    if (mesize.get_size() == -1) mesize.set_size(1.0);
    
    // assemble bounding mesh
    std::vector<moab::EntityHandle> bdy;
    std::vector<int> group_sizes, bdy_ids;
    me->boundary(0, bdy, NULL, &group_sizes);

      // convert the handles to integers for input to TriAdv
    moab::Range bdy_vrange;
    me->handles_to_ids(bdy, bdy_ids, NULL, &bdy_vrange);

      // get the vertex coordinates
    std::vector<double> coords(3*bdy_vrange.size());
    moab::ErrorCode rval = mk_core()->moab_instance()->get_coords(bdy_vrange, &coords[0]);
    MBERRCHK(rval, mk_core()->moab_instance());

    if (coords.size() != 3*bdy_vrange.size())
      ECERRCHK(MK_FAILURE, "Trouble getting vertex coordinates.");
    
      // now construct the CAMAL mesher, and pass it initial conditions
    CMLTriAdvance triadv(&cse, &mesize);
    bool success = triadv.set_boundary_mesh(bdy_vrange.size(), &coords[0], group_sizes.size(), &group_sizes[0], &bdy_ids[0]);
    if (!success)
      ECERRCHK(MK_FAILURE, "Trouble setting boundary mesh.");

      // ok, now generate the mesh
    int num_pts, num_tris;
    success = triadv.generate_mesh(num_pts, num_tris);
    if (!success)
      ECERRCHK(MK_FAILURE, "Trouble generating tri mesh.");

    moab::Range &new_ents = (*sit).second;

      // resize the coords array, then get the coords of the new points
    assert(num_pts >= (int)bdy_vrange.size());
    if (num_pts > (int)bdy_vrange.size()) {
      coords.resize(3*(num_pts-bdy_vrange.size()));
      int pts_returned = triadv.get_points_buf(coords.size(), &coords[0], bdy_vrange.size());
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
    moab::EntityHandle starth, *connect, *tmp_connect;
    rval = iface->get_element_connect(num_tris, 3, moab::MBTRI, 1, starth, connect);
    MBERRCHK(rval, mk_core()->moab_instance());

      // read connectivity directly into that array, as int's
    int *connecti = (int*) connect;
    int tris_returned = triadv.get_tris_buf(3*num_tris, connecti);
    if (tris_returned != num_tris)
      ECERRCHK(MK_FAILURE, "Number of new tris returned from TriAdv doesn't agree with previous value output.");

      // put vertex handles into an indexible array
    std::vector<moab::EntityHandle> verts;
    verts.reserve(bdy_vrange.size() + new_ents.size());
    std::copy(bdy_vrange.begin(), bdy_vrange.end(), std::back_inserter(verts));
    std::copy(new_ents.begin(), new_ents.end(), std::back_inserter(verts));
    
      // now convert vertex indices into handles in-place, working from the back
    for (int i = 3*num_tris-1; i > 0; i--) {
      assert(connecti[i] >= 0 && connecti[i] < (int)verts.size());
      connect[i] = verts[connecti[i]];
    }
    
      // put new tris into new entity range, then commit the mesh
    new_ents.insert(starth+num_tris-1);
    me->commit_mesh(new_ents, COMPLETE_MESH);
  }
}

MeshOp *CAMALTriAdvance::factory(MKCore *mkcore, const MEntVector &me_vec)
{
    // construct a CAMALTriAdvance and pass it back
  return new CAMALTriAdvance(mkcore, me_vec);
}

} // namespace MeshKit
