
#include "meshkit/MKCore.hpp"
#include "meshkit/CurveMesher.hpp"
#include "meshkit/ModelEnt.hpp"
#include <iostream>
#include <iGeom.h>

#include "MBTagConventions.hpp"
#include "MBCore.hpp"
#include "moab/GeomTopoTool.hpp"

namespace MeshKit
{
// Construction Function for CurveMesher

moab::EntityType CurveMesher_types[] = {moab::MBVERTEX, moab::MBEDGE, moab::MBMAXTYPE};
const moab::EntityType* CurveMesher::output_types()
  { return CurveMesher_types; }


CurveMesher::CurveMesher(MKCore *mk_core, const MEntVector &ments)
  : MeshScheme(mk_core, ments)
{
  mk = mk_core;
  model_ents = ments;
}

// Destructor Function for Curvemesher
CurveMesher::~CurveMesher()
{
}

void CurveMesher::setup_this()
{

  MEntVector::iterator i;
  for ( i = model_ents.begin(); i != model_ents.end(); i++)
    { 
      ModelEnt *me = *i;
      //do an initial check that the model entity is of the correct dimension
      if( me->dimension() != 1)
	{
	  std::cout << "Found an entity that is of an incorrect dimension" << std::endl;
          model_ents.erase(i);
          continue;
        } 

    }

  //set the faceting_tolerance
  facet_tol = 1e-3;

}

void CurveMesher::execute_this()
{
 
  //setup the category tag value  
  char category[CATEGORY_TAG_SIZE] = "Curve\0";


  //create iMesh tag for categories
  iBase_TagHandle category_tag;
  mk_core()->imesh_instance()->createTag(CATEGORY_TAG_NAME,CATEGORY_TAG_SIZE,iBase_BYTES, category_tag);

  MEntVector::iterator i;
  for( i = model_ents.begin(); i != model_ents.end(); i++)
    {

      ModelEnt *me = *i;
      
      //get the mesh set handle from the ModelEnt
      iMesh::EntitySetHandle msh = IBSH(me->mesh_handle());

      //get the geom handle from the ModelEnt
      iGeom::EntityHandle gh = me->geom_handle();

      //---------------------------------//
      //           SET TAGS              //
      //---------------------------------//
      //set the category tag
      mk_core()->imesh_instance()->setEntSetData(msh, category_tag, &category);
    
      //---------------------------------//
      //            FACET                //
      //---------------------------------//
   
      facet(gh,msh);
      //set_senses(me);
      me->set_meshed_state(COMPLETE_MESH);

    }

}

 void CurveMesher::facet( iGeom::EntityHandle h, iMesh::EntitySetHandle sh)
{

      
      //get the facets for this curve/ref_edge
      std::vector<double> pnts;
      std::vector<int> conn;
      mk_core()->igeom_instance()->getFacets(h,facet_tol,pnts,conn);
      std::cout << "Points returned from getFacets: " << pnts.size()/3 << std::endl;
      std::cout << "Facets returned from getFacets: " << conn.size() << std::endl;
      
      //create vector for keeping track of the vertices
      std::vector<iBase_EntityHandle> verts;
 
      // loop over the facet points
      for(unsigned int j=0; j<pnts.size(); j+=3)
	{
	  //create a new vertex handle 
          iBase_EntityHandle v;
          mk_core()->imesh_instance()->createVtx(pnts[j],pnts[j+1],pnts[j+2],v);
          verts.push_back(v);
	}

      //vector to contain the edges
      std::vector<iBase_EntityHandle> edges;
      //loop over vertices and create edges
      for (unsigned int j=0; j<verts.size()-1; j++)
	{
	  //create edges
          iBase_EntityHandle e; 
          mk_core()->imesh_instance()->createEnt(iMesh_LINE_SEGMENT, &verts[j],2,e);
	  edges.push_back(e);
	}

      //add vertices and edges to the entity set
      mk_core()->imesh_instance()->addEntArrToSet(&verts[0],verts.size(),sh);
      mk_core()->imesh_instance()->addEntArrToSet(&edges[0],edges.size(),sh);


}

  void CurveMesher::set_senses(ModelEnt *ent)
{
  
  //get the geom_handle for this curve
  iGeom::EntityHandle gh = ent->geom_handle();

  //get all surfaces adjacent to this curve
  MEntVector adj_surfs;
  ent->get_adjacencies(2,adj_surfs);

  MEntVector::iterator i; 

  std::vector<int> senses;
  std::vector<moab::EntityHandle> meshsets;
  for(i = adj_surfs.begin(); i !=adj_surfs.end(); i++)
    {
      int sense;
      ModelEnt *adj_ent = *i;
      mk_core()->igeom_instance()->getEgFcSense(gh, adj_ent->geom_handle(),sense);
      std::cout << "Sense: " << sense << std::endl;
      senses.push_back(sense);

      meshsets.push_back(adj_ent->mesh_handle());

    }

   
  moab::GeomTopoTool gt(mk_core()->moab_instance());
  moab::ErrorCode rval = gt.set_senses(ent->mesh_handle(),meshsets,senses);
  if(rval != MB_SUCCESS) std::cout << "Error setting curve senses!" << std::endl;

}

}
