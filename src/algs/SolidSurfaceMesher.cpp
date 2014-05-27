

#include "meshkit/MKCore.hpp"
#include "meshkit/SolidSurfaceMesher.hpp"
#include "meshkit/ModelEnt.hpp"
#include <iostream>
#include <iGeom.h>

#include "MBTagConventions.hpp"
#include "MBCore.hpp"
#include "moab/GeomTopoTool.hpp"


namespace MeshKit
{

// Output mesh types for this class
moab::EntityType SolidSurfaceMesher_types[] = {moab::MBVERTEX, moab::MBEDGE, moab::MBTRI, moab::MBMAXTYPE};
const moab::EntityType* SolidSurfaceMesher::output_types()
{ return SolidSurfaceMesher_types; }

SolidSurfaceMesher::SolidSurfaceMesher(MKCore *mk_core, const MEntVector &ments)
  : MeshScheme(mk_core, ments)
{
  mk = mk_core; 
  model_ents = ments;
}

SolidSurfaceMesher::~SolidSurfaceMesher()
{
}

void SolidSurfaceMesher::setup_this()
{

  //create a solid curve mesher 
  MeshOp *scm = (MeshOp*) mk_core()->construct_meshop("CurveMesher");
 
  MEntVector::iterator i; 
  for(i = model_ents.begin(); i != model_ents.end(); i++)
    {
      ModelEnt *me = *i;
      //do an initial check that the model entity is of the correct dimension
      if(me->dimension() != 2)
	{
	  std::cout << "Found a model entity with an incorrect dimension for this mesher" << std::endl;
          model_ents.erase(i); 
          continue; 
	}

      //make sure that the model entity's children have been meshed
      MEntVector children; 
      me->get_adjacencies(1, children); 
 
      for(MEntVector::iterator j = children.begin(); j != children.end(); j++)
	{
	  ModelEnt *child = *j; 
          if(child->is_meshops_list_empty()) scm-> add_modelent(child); 
	}

    }

  //have the children be meshed first
  mk_core()->insert_node(scm, this, mk_core()->root_node());

  //set the faceting parameters to default values if they are not already set
  
  //set the faceting tolerance
  if(!facet_tol) facet_tol = 1e-4; 
  
  //set the geom_resabs value
  if(!geom_res) geom_res = 1e-6;

}

void SolidSurfaceMesher::execute_this()
{

  //setup the category tag value
  char category[CATEGORY_TAG_SIZE] = "Surface\0";
 
  //create iMesh tag for the category
  iBase_TagHandle category_tag; 
  mk_core()->imesh_instance()->createTag(CATEGORY_TAG_NAME, CATEGORY_TAG_SIZE, iBase_BYTES, category_tag); 

  MEntVector::iterator i; 
  for(i = model_ents.begin(); i !=model_ents.end(); i++)
    {
      ModelEnt *me = *i; 
      
      //get the mesh set handle from the ModelEnt
      iMesh::EntitySetHandle msh = IBSH(me->mesh_handle());

      //---------------------------------//
      //           SET TAGS              //
      //---------------------------------//
      //set the category tag
      mk_core()->imesh_instance()->setEntSetData(msh, category_tag, &category);
    
      //---------------------------------//
      //            FACET                //
      //---------------------------------//
   
      //facet(me);

      me->set_meshed_state(COMPLETE_MESH);
    }
}

void SolidSurfaceMesher::facet(ModelEnt *surf)
{

  iGeom::EntityHandle h = surf->geom_handle(); 
  iMesh::EntitySetHandle sh = IBSH(surf->mesh_handle()); 

  //get the facets for this surface/ ref_face
  std::vector<double> pnts;
  std::vector<int> conn;
  mk->igeom_instance()->getFacets(h,facet_tol,pnts,conn);
  std::cout << "Triangles returned from getFacets: " << pnts.size()/3 << std::endl;
  std::cout << "Facets returned from getFacets: " << conn.size()/3 << std::endl;
  
  //create vector for keeping track of the vertices
  std::vector<iBase_EntityHandle> verts;

  // loop over the facet points
  for(unsigned int j=0; j<pnts.size(); j+=3)
    {
      //create a new vertex handle 
      iBase_EntityHandle v;
      mk->imesh_instance()->createVtx(pnts[j],pnts[j+1],pnts[j+2],v);
      verts.push_back(v);
    }
      
  //vector for keeping track of the triangles
  std::vector<iBase_EntityHandle> tris;

  //loop over the connectivity
  for(unsigned int j=0; j<conn.size()-2; j+=3)
    {
      //get the appropriate points for a triangle and add them to a vector
      std::vector<iBase_EntityHandle> tri_verts; 
      tri_verts.push_back(verts[conn[j]]);
      tri_verts.push_back(verts[conn[j+1]]);
      tri_verts.push_back(verts[conn[j+2]]);

      //create a new triangle handle
      iBase_EntityHandle t; 
      mk->imesh_instance()->createEnt(iMesh_TRIANGLE, &tri_verts[0],3,t);
      tri_verts.clear();
      tris.push_back(t);
    }
  std::cout << "Created " << tris.size() << " triangles" << std::endl;

  //add verticess and edges to the entity set
  mk->imesh_instance()->addEntArrToSet(&verts[0],verts.size(),sh);
  mk->imesh_instance()->addEntArrToSet(&tris[0],tris.size(),sh);

}

}
