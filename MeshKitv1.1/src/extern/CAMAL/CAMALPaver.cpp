#include "meshkit/CAMALPaver.hpp"
#include "meshkit/MeshScheme.hpp"
#include "meshkit/ModelEnt.hpp"
#include "moab/Interface.hpp"
#include "moab/ReadUtilIface.hpp"
#include "moab/Range.hpp"
#include "meshkit/RegisterMeshOp.hpp"
#include "CAMALSurfEval.hpp"
#include "CAMALSizeEval.hpp"
#include "CMLPaver.hpp"
#include <vector>

namespace MeshKit
{

moab::EntityType CAMALPaver::meshTps[] = {moab::MBVERTEX, moab::MBQUAD, moab::MBMAXTYPE};

CAMALPaver::CAMALPaver(MKCore *mk_core, const MEntVector &me_vec)
        : MeshScheme(mk_core, me_vec)
{
}

CAMALPaver::~CAMALPaver()
{
}

void CAMALPaver::setup_this()
{
    // need to constrain all edges to be even
  constrain_even();
  
    // then call setup_boundary, to set up edge meshers
  setup_boundary();
}

void CAMALPaver::execute_this()
{

#ifdef HAVE_FBIGEOM
  if (mentSelection.empty())
  {
    // create model ents from previous op
    create_model_ents_from_previous_ops();
    // now look at the latest SizingFunction, and set it or each model ent
    int latestIndexSF = 0; // maybe we would need to set it right
    for (MEntSelection::iterator sit = mentSelection.begin();
        sit != mentSelection.end(); sit++) {
        // make a me, for convenience
        ModelEnt *me = (*sit).first;
        me->sizing_function_index(latestIndexSF); // need to work on this one; how do we know?
    }
    // now, force setup of this node again, as we have added model entities to it
    setup_called(false);
    mk_core()->setup(false);
    // debug
    mk_core()->print_graph();
    // it may not be enough, we may have to execute the previous ops, that were just
    // created during setup... not very clean code;
    mk_core()->execute_before((GraphNode *) this);
  }
#endif
  for (MEntSelection::iterator sit = mentSelection.begin(); sit != mentSelection.end(); sit++) {
      // make a me, for convenience
    ModelEnt *me = (*sit).first;

    if(me->get_meshed_state()>=COMPLETE_MESH)
      continue;
      // create a surface evaluator for this modelent, and a size evaluator
    CAMALSurfEval cse(me);

    SizingFunction * sz = mk_core()->sizing_function(me->sizing_function_index());
    if (!sz->variable())
      sz = NULL; // no function variable

    CAMALSizeEval mesize(me->mesh_interval_size(), sz);
          // make sure the size isn't negative
      // make sure the size isn't negative
    if (mesize.get_size() == -1) mesize.set_size(1.0);
    
    // assemble bounding mesh
    std::vector<moab::EntityHandle> bdy;
    std::vector<int> group_sizes, bdy_ids;
    me->boundary(0, bdy, NULL, &group_sizes);

      // convert the handles to integers for input to Paver
    moab::Range bdy_vrange;
    std::vector<double> coords;
    me->get_indexed_connect_coords(bdy, NULL, NULL, bdy_ids, coords, &bdy_vrange);

      // now construct the CAMAL mesher, and pass it initial conditions
    CMLPaver paver(&cse, &mesize);
    bool success = paver.set_boundary_mesh(bdy_vrange.size(), &coords[0], group_sizes.size(), &group_sizes[0], &bdy_ids[0]);
    if (!success)
      ECERRCHK(MK_FAILURE, "Trouble setting boundary mesh.");

      // ok, now generate the mesh
    int num_pts, num_quads;
    success = paver.generate_mesh(num_pts, num_quads);
    if (!success)
      ECERRCHK(MK_FAILURE, "Trouble generating quad mesh.");

    moab::Range &new_ents = (*sit).second;
    moab::ErrorCode rval;

      // resize the coords array, then get the coords of the new points
    assert(num_pts >= (int)bdy_vrange.size());
    if (num_pts > (int)bdy_vrange.size()) {
      coords.resize(3*(num_pts-bdy_vrange.size()));
      int pts_returned = paver.get_points_buf(coords.size(), &coords[0], bdy_vrange.size());
      if (pts_returned != num_pts-(int)bdy_vrange.size())
        ECERRCHK(MK_FAILURE, "Number of new points returned from Paver doesn't agree with previous value output.");
    
        // create the new vertices' entities on the face
      rval = mk_core()->moab_instance()->create_vertices(&coords[0], pts_returned, new_ents);
      MBERRCHK(rval, mk_core()->moab_instance());
    }

      // for quads, pre-allocate connectivity
    moab::ReadUtilIface *iface;
    rval = mk_core()-> moab_instance() -> query_interface(iface);
    MBERRCHK(rval, mk_core()->moab_instance());		

      //create the quads, get a direct ptr to connectivity
    moab::EntityHandle starth, *connect;
    rval = iface->get_element_connect(num_quads, 4, moab::MBQUAD, 1, starth, connect);
    MBERRCHK(rval, mk_core()->moab_instance());

      // read connectivity directly into that array, as int's
    int *connecti = (int*) connect;
    int quads_returned = paver.get_quads_buf(4*num_quads, connecti);
    if (quads_returned != num_quads)
      ECERRCHK(MK_FAILURE, "Number of new quads returned from Paver doesn't agree with previous value output.");

      // put vertex handles into an indexible array
    std::vector<moab::EntityHandle> bdy_vvec;
    std::copy(bdy_vrange.begin(), bdy_vrange.end(), std::back_inserter(bdy_vvec));
    std::copy(new_ents.begin(), new_ents.end(), std::back_inserter(bdy_vvec));
    
      // now convert vertex indices into handles in-place, working from the back
    for (int i = 4*num_quads-1; i >= 0; i--) {
      assert(connecti[i] >= 0 && connecti[i] < (int)bdy_vvec.size());
      connect[i] = bdy_vvec[connecti[i]];
    }
    
      // put new quads into new entity range, then commit the mesh
    new_ents.insert(starth, starth+num_quads-1);
    me->commit_mesh(new_ents, COMPLETE_MESH);
  }
}

} // namespace MeshKit
