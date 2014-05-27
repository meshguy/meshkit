
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

  //create a vertex mesher
  int meshop_dim = 0;
  MeshOp *vm = (MeshOp*) mk_core()->construct_meshop(meshop_dim);

  MEntVector::iterator i;
  for (i = model_ents.begin(); i != model_ents.end(); i++)
    { 
      ModelEnt *me = *i;
      //do an initial check that the model entity is of the correct dimension
      if(me->dimension() != 1)
	{
	  std::cout << "Found an entity that is of an incorrect dimension" << std::endl;
          model_ents.erase(i);
          continue;
        } 

     
      //make sure that its children (vertices) have been meshed
      MEntVector children; 
      me->get_adjacencies(0, children);

      for(MEntVector::iterator j = children.begin(); j != children.end(); j++)
	{
          ModelEnt *child = *j;
          if(child->is_meshops_list_empty()) vm->add_modelent(child);
	}

    }

  //make sure we're meshing the veritces first
  mk_core()->insert_node(vm, this, mk_core()->root_node());
  
  //set faceting parameters to default values if the are not already set. 

  //set the faceting_tolerance
  if(!facet_tol) facet_tol = 1e-4;

  //set the geom_reabs value
  if(!geom_res) geom_res = 1e-6;

}

void CurveMesher::execute_this()
{
 
  //setup the category tag value  
  char category[CATEGORY_TAG_SIZE] = "Curve\0";


  //create iMesh tag for categories
  iBase_TagHandle category_tag;
  mk_core()->imesh_instance()->createTag(CATEGORY_TAG_NAME,CATEGORY_TAG_SIZE,iBase_BYTES, category_tag);

  MEntVector::iterator i;
  for(i = model_ents.begin(); i != model_ents.end(); i++)
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
   
      facet(me);
      // Not sure if this needs to be done yet
      //set_senses(me);
      me->set_meshed_state(COMPLETE_MESH);
    }
}

 void CurveMesher::facet(ModelEnt *curve)
{
  iGeom::EntityHandle h = curve->geom_handle();
  iMesh::EntitySetHandle sh = IBSH(curve->mesh_handle());

  //get the facets for this curve/ref_edge
  std::vector<double> pnts;
  std::vector<int> conn;
  mk_core()->igeom_instance()->getFacets(h,facet_tol,pnts,conn);
  std::cout << "Points returned from getFacets: " << pnts.size()/3 << std::endl;
  std::cout << "Facets returned from getFacets: " << conn.size() << std::endl;


  //create vector for keeping track of the vertices
  std::vector<iBase_EntityHandle> verts;
          

  // loop over the facet points
  for(unsigned int j = 0; j < pnts.size(); j+=3)
    {
      //create a new vertex handle 
      iBase_EntityHandle v;
      mk_core()->imesh_instance()->createVtx(pnts[j],pnts[j+1],pnts[j+2],v);
      verts.push_back(v);
    }
  //std::cout << "Number of verts created: " << verts.size() << std::endl;
  //--------------------------------------------------//
  //check that the start and end vertices are within  //
  //the geom resabs of the facet verts                //
  //--------------------------------------------------//
     
  //get the curve vertices
  MEntVector end_verts;
  curve->get_adjacencies(0, end_verts);
  std::cout << "Number of end_verts: " << end_verts.size() << std::endl;

  //check for a closed curve
  bool closed = (end_verts.size() == 1);

  //Special case for a point curve
  if (verts.size() < 2 && (!closed || curve->measure() > geom_res))
    {
      std::cerr << "Warning: No facetting for curve" << curve->id() << std::endl;
      return;
    }


  if(closed)
    {
      //check the proximity of the single vertex to both the start and end points
      iGeom::EntityHandle end_vert_hand = end_verts[0]->geom_handle(); 
      if(vtx2vtx_dist(end_vert_hand,verts.front()) > geom_res || vtx2vtx_dist(end_vert_hand,verts.back()) > geom_res)
        std::cerr << "Warning: closed curve vertex not at the end of curve facets" << std::endl;

      //capture the beginning and end vertex handles for deletion 
      iMesh::EntityHandle front_vert_to_del = verts.front();
      iMesh::EntityHandle back_vert_to_del = verts.back();

      //insert the single vertex in at the beginning and the end of the curve
      std::vector<moab::EntityHandle> dum_handles;
      end_verts[0]->get_mesh(0,dum_handles);
      verts.front() = IBEH(dum_handles[0]);
      verts.back() = IBEH(dum_handles[0]);
     
      //delete the old vertices
      mk_core()->imesh_instance()->deleteEnt(front_vert_to_del);
      mk_core()->imesh_instance()->deleteEnt(back_vert_to_del);

    }
  else
    {
      //check the proximity of the front and end vertices to the start and end points, respectively
      if(vtx2vtx_dist(end_verts[0]->geom_handle(), verts.front()) > geom_res ||
         vtx2vtx_dist(end_verts[1]->geom_handle(), verts.back()) > geom_res)
        {
          //try reversing the points
          std::reverse(verts.begin(), verts.end());
        }
          //check again, if this time it fails, give a warning
          if(vtx2vtx_dist(end_verts[0]->geom_handle(), verts.front()) > geom_res ||
             vtx2vtx_dist(end_verts[1]->geom_handle(), verts.back()) > geom_res)
            {
              std::cerr << "Warning: vertices not at the ends of the curve" << std::endl;
            }
      //now replace start and end facet points with the child vert

      //capture the beginning and end vertex handles for deletion 
      iMesh::EntityHandle front_vert_to_del = verts.front();
      iMesh::EntityHandle back_vert_to_del = verts.back();


      //this is all necessary to keep withing imesh interface
      std::vector<moab::EntityHandle> dum_handles;
      end_verts[0]->get_mesh(0,dum_handles);
      verts.front() = IBEH(dum_handles[0]);
      dum_handles.clear();
      end_verts[1]->get_mesh(0,dum_handles);
      verts.back() = IBEH(dum_handles[0]);

      //delete the old vertices
      mk_core()->imesh_instance()->deleteEnt(front_vert_to_del);
      mk_core()->imesh_instance()->deleteEnt(back_vert_to_del);

    }
    
  //now create the edges with the vertices we've chosen

  //vector to contain the edges
  std::vector<iBase_EntityHandle> edges;
  //loop over vertices and create edges
  for (unsigned int j = 0; j < verts.size()-1; j++)
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

  void CurveMesher::set_senses(ModelEnt *curve)
{
 
  //get the geom_handle for this curve
  iGeom::EntityHandle gh = curve->geom_handle();

  //get all surfaces adjacent to this curve
  MEntVector adj_surfs;
  curve->get_adjacencies(2, adj_surfs);

  MEntVector::iterator i; 
  std::vector<int> senses;
  std::vector<moab::EntityHandle> meshsets;
  for(i = adj_surfs.begin(); i !=adj_surfs.end(); i++)
    {
      int sense;
      //create a model ent for the current adjacent surface
      ModelEnt *adj_ent = *i;
      //get the sense of the curve wrt the surface
      mk_core()->igeom_instance()->getEgFcSense(gh, adj_ent->geom_handle(),sense);
      //std::cout << "Sense: " << sense << std::endl;
      // add the sense and the entityhandle to the appropriate vectors
      senses.push_back(sense);
      meshsets.push_back(adj_ent->mesh_handle());
    }

  // use moab to set the senses of the entities
  // using moab to do this allows the use of variable length tags
  // which are not currently available in iGeom 
  moab::GeomTopoTool gt(mk_core()->moab_instance());
  moab::ErrorCode rval = gt.set_senses(curve->mesh_handle(),meshsets,senses);
  if(rval != MB_SUCCESS) std::cout << "Error setting curve senses!" << std::endl;

}

// should I be using an iMesh entity to do this comparison??
double CurveMesher::vtx2vtx_dist( iGeom::EntityHandle vtx1, iMesh::EntityHandle vtx2)
{

  double x1,y1,z1;
  double x2,y2,z2;

  mk_core()->igeom_instance()->getVtxCoord( vtx1, x1, y1, z1);
  mk_core()->imesh_instance()->getVtxCoord( vtx2, x2, y2, z2);

  double dist = pow((x1-x2),2) + pow((y1-y2),2) + pow((z1-z2),2);

  return sqrt(dist);
}

void CurveMesher::set_mesh_params(double faceting_tolerance, double geom_resabs)
{
  //Assign the faceting values if they are passed into the function
  if(faceting_tolerance) facet_tol = faceting_tolerance; 
  if(geom_resabs) geom_res = geom_resabs; 
}

}
