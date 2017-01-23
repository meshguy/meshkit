#include "meshkit/NGTetMesher.hpp"
#include "meshkit/MeshScheme.hpp"
#include "meshkit/ModelEnt.hpp"
#include "moab/Interface.hpp"
#include "moab/ReadUtilIface.hpp"
#include "moab/Range.hpp"
#include "meshkit/RegisterMeshOp.hpp"
namespace nglib 
{
#include "nglib.h"
}

#include <vector>

namespace MeshKit
{

moab::EntityType NGTetMesher::meshTps[] = {moab::MBVERTEX, moab::MBTET, moab::MBMAXTYPE};

NGTetMesher::NGTetMesher(MKCore *mk_core, const MEntVector &me_vec)
        : MeshScheme(mk_core, me_vec)
{
}

NGTetMesher::~NGTetMesher()
{
}

MeshOp *NGTetMesher::get_tri_mesher() 
{
  std::vector<MeshOpProxy *> proxies;
    mk_core()->meshop_by_mesh_type(moab::MBTRI, proxies);
   if (proxies.empty()) throw Error(MK_FAILURE, "Couldn't find a MeshOp capable of producing triangles.");
  return mk_core()->construct_meshop(*proxies.begin());
  return mk_core()->construct_meshop("CAMALTriAdvance");
}

void NGTetMesher::setup_this()
{
  MeshOp *tri_mesher = NULL;
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
      if ((*vit)->is_meshops_list_empty()) {
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

void NGTetMesher::execute_this()
{
  nglib::Ng_Init();

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
    me->get_indexed_connect_coords(bdy, &senses, NULL, bdy_ids, coords, &bdy_vrange, 1);

    nglib::Ng_Mesh *ngmesh = nglib::Ng_NewMesh();
    unsigned int numvs = bdy_vrange.size();
    
      // add the points
    for (unsigned int i = 0; i < numvs; i++)
      nglib::Ng_AddPoint(ngmesh, &coords[3*i]);

      // add the tris
    unsigned int numtris = bdy.size();
    for (unsigned int i = 0; i < numtris; i++)
      nglib::Ng_AddSurfaceElement(ngmesh, nglib::NG_TRIG, &bdy_ids[3*i]);
    
    double my_size = me->mesh_interval_size();
    if (0 > my_size) my_size = 1.0;
    nglib::Ng_RestrictMeshSizeGlobal(ngmesh, my_size);

    nglib::Ng_Meshing_Parameters ngp;
    ngp.maxh = my_size;
    ngp.fineness = 0.5;
    
    // Use variable name second_order if using newer version of NetGen
    //ngp.secondorder = 0;
        ngp.second_order = 0;
    //    nglib.h Rev. 663 http://netgen-mesher.svn.sourceforge.net/viewvc/netgen-mesher/netgen/nglib/nglib.h?revision=663&view=markup
    //uses: second_order instead of secondorder(4.9.13) 

    
    nglib::Ng_Result result = nglib::Ng_GenerateVolumeMesh(ngmesh, &ngp);
    if (nglib::NG_OK != result) ECERRCHK(MK_FAILURE, "Netgen mesher returned !ok.");
    
    moab::ErrorCode rval;
    moab::Range new_ents;
    int num_pts = nglib::Ng_GetNP(ngmesh);
    assert(num_pts >= (int)numvs);
    if (num_pts > (int)numvs) {
      coords.resize(3*(num_pts-numvs));
      
      for (unsigned int i = (unsigned int) numvs; i < (unsigned int) num_pts; i++)
        nglib::Ng_GetPoint(ngmesh, i, &coords[3*(i-numvs)]);
      
        // create the new vertices' entities 
      rval = mk_core()->moab_instance()->create_vertices(&coords[0], num_pts-numvs, new_ents);
      MBERRCHK(rval, mk_core()->moab_instance());
    }

      // for tets, pre-allocate connectivity
    moab::ReadUtilIface *iface;
    rval = mk_core()-> moab_instance() -> query_interface(iface);
    MBERRCHK(rval, mk_core()->moab_instance());		

      //create the tris, get a direct ptr to connectivity
    int num_tets = (int) nglib::Ng_GetNE(ngmesh);
    moab::EntityHandle starth, *connect;
    rval = iface->get_element_connect(num_tets, 4, moab::MBTET, 1, starth, connect);
    MBERRCHK(rval, mk_core()->moab_instance());

      // read connectivity directly into that array, as int's
    int *connecti = (int*) connect;
    for (unsigned int i = 0; i < (unsigned int) num_tets; i++)
      nglib::Ng_GetVolumeElement(ngmesh, i+1, connecti+4*i);

      // put vertex handles into an indexible array
    std::vector<moab::EntityHandle> bdy_vvec;
    std::copy(bdy_vrange.begin(), bdy_vrange.end(), std::back_inserter(bdy_vvec));
    std::copy(new_ents.begin(), new_ents.end(), std::back_inserter(bdy_vvec));
    
      // now convert vertex indices into handles in-place, working from the back
    for (int i = 4*num_tets-1; i >= 0; i--) {
      assert(connecti[i] > 0 && connecti[i] <= (int)bdy_vvec.size());
      connect[i] = bdy_vvec[connecti[i]-1];
    }
    
      // put new tris into new entity range, then commit the mesh
    new_ents.insert(starth, starth+num_tets-1);
    me->commit_mesh(new_ents, COMPLETE_MESH);

    nglib::Ng_DeleteMesh(ngmesh);
  }

  nglib::Ng_Exit();
}

} // namespace MeshKit
