
#include "meshkit/MKCore.hpp"
#include "meshkit/CurveMesher.hpp"
#include "meshkit/ModelEnt.hpp"
#include <iostream>
#include <iGeom.h>

#include "MBTagConventions.hpp"

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

  //setup the geom dimension tag value
  int dim = 1;

  //create iMesh tag for categories
  iBase_TagHandle category_tag;
  mk_core()->imesh_instance()->createTag(CATEGORY_TAG_NAME,CATEGORY_TAG_SIZE,iBase_BYTES, category_tag);

  //create iMesh tag for geom_dim
  iBase_TagHandle geom_tag; 
  mk_core()->imesh_instance()->createTag(GEOM_DIMENSION_TAG_NAME,1,iBase_INTEGER, geom_tag);

  MEntVector::iterator i;
  for( i = model_ents.begin(); i != model_ents.end(); i++)
    {

      ModelEnt *me = *i;
      
      //get the mesh set handle from the ModelEnt
      iMesh::EntitySetHandle msh = IBSH(me->mesh_handle());

      //---------------------------------//
      //           SET TAGS              //
      //---------------------------------//
      //set the geom_dim tag
      //mk_core()->imesh_instance()->setEntSetIntData(msh,geom_tag,dim);

      //set the category tag
      mk_core()->imesh_instance()->setEntSetData(msh, category_tag, &category);

      //get the geom entity handle from the ModelEnt
      iGeom::EntityHandle gh = me->geom_handle();

      //get the facets for this curve/ref_edge
      std::vector<double> pnts;
      std::vector<int> conn;
      mk_core()->igeom_instance()->getFacets(gh,facet_tol,pnts,conn);
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

      //add verticess and edges to the entity set
      mk_core()->imesh_instance()->addEntArrToSet(&verts[0],verts.size(),msh);
      mk_core()->imesh_instance()->addEntArrToSet(&edges[0],edges.size(),msh);

    }

}


}
