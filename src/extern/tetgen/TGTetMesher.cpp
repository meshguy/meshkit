#include "meshkit/TGTetMesher.hpp"
#include "meshkit/MeshScheme.hpp"
#include "meshkit/ModelEnt.hpp"
#include "moab/Interface.hpp"
#include "moab/ReadUtilIface.hpp"
#include "moab/Range.hpp"
#include "meshkit/RegisterMeshOp.hpp"

namespace tetgen {
  #include "tetgen.h"
}


#include <vector>

namespace MeshKit
{

moab::EntityType TGTetMesher::meshTps[] = {moab::MBVERTEX, moab::MBTET, moab::MBMAXTYPE};

TGTetMesher::TGTetMesher(MKCore *mk_core, const MEntVector &me_vec)
        : MeshScheme(mk_core, me_vec)
{
  tet_switch = (char*)"pq1.414Y";
}

TGTetMesher::~TGTetMesher()
{
}

MeshOp *TGTetMesher::get_tri_mesher()
{
  std::vector<MeshOpProxy *> proxies;
  //  mk_core()->meshop_by_mesh_type(moab::MBTRI, proxies);
  // if (proxies.empty()) throw Error(MK_FAILURE, "Couldn't find a MeshOp capable of producing triangles.");
  //return mk_core()->construct_meshop(*proxies.begin());
  return mk_core()->construct_meshop("CAMALTriAdvance");
}

void TGTetMesher::setup_this()
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

void TGTetMesher::execute_this()
{
  tetgen::tetgenio in, out;
  tetgen::tetgenio::facet *f;
  tetgen::tetgenio::polygon *p;
  in.firstnumber = 0;

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

    // add points
    in.numberofpoints = bdy_vrange.size();
    in.pointlist = new REAL[3 * bdy_vrange.size()];
    for (int i = 0; i < (int)bdy_vrange.size(); i++)
    {
        in.pointlist[3 * i + 0] = coords[3*i];
        in.pointlist[3 * i + 1] = coords[3*i+1];
        in.pointlist[3 * i + 2] = coords[3*i+2];
      }

     // add the tris
     unsigned int numtris = bdy.size();

     in.numberoffacets = numtris;
     in.facetlist = new tetgen::tetgenio::facet[numtris];
     in.facetmarkerlist = NULL;

    for (unsigned int i = 0; i < numtris; i++){
        f = &in.facetlist[i];
        f->numberofpolygons = 1;
        f->polygonlist = new tetgen::tetgenio::polygon[1];
        f->numberofholes = 0;
        f->holelist = NULL;
        p = &f->polygonlist[0];
        p->numberofvertices = 3;
        p->vertexlist = new int[3];
        p->vertexlist[0] = bdy_ids[3*i];
        p->vertexlist[1] = bdy_ids[3*i + 1];
        p->vertexlist[2] = bdy_ids[3*i + 2];
      }

    // CALL TETGEN
    tetgen::tetrahedralize(tet_switch, &in, &out);

   // TODO: Bring back mesh to MeshKit database
    }
}

} // namespace MeshKit
