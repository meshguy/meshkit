
#include "meshkit/MKCore.hpp"
#include "meshkit/CurveFacetMeshReader.hpp"
#include "meshkit/ModelEnt.hpp"
#include <iostream>
#include <iGeom.h>
#include <cmath>
#include "MBTagConventions.hpp"
#include "MBCore.hpp"
#include "moab/GeomTopoTool.hpp"

namespace MeshKit
{
// Construction Function for CurveFacetMeshReader

moab::EntityType CurveFacetMeshReader_types[] = {moab::MBVERTEX, moab::MBEDGE, moab::MBMAXTYPE};
const moab::EntityType* CurveFacetMeshReader::output_types()
  { return CurveFacetMeshReader_types; }


CurveFacetMeshReader::CurveFacetMeshReader(MKCore *mk_core, const MEntVector &ments)
  : MeshScheme(mk_core, ments)
{
  mk = mk_core;
}

// Destructor Function for CurveFacetMeshReader
CurveFacetMeshReader::~CurveFacetMeshReader()
{
}

void CurveFacetMeshReader::setup_this()
{

  //create a vertex mesher
  int meshop_dim = 0;
  MeshOp *vm = (MeshOp*) mk_core()->construct_meshop(meshop_dim);

  for (MEntSelection::iterator mit = mentSelection.begin(); mit != mentSelection.end(); mit++)
    { 
      ModelEnt *me = mit->first;
      //do an initial check that the model entity is of the correct dimension
      if(me->dimension() != 1)
	{
	  std::cout << "Found an entity that is of an incorrect dimension" << std::endl;
          mentSelection.erase(mit);
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

void CurveFacetMeshReader::execute_this()
{
 
  for(MEntSelection::iterator mit = mentSelection.begin(); mit != mentSelection.end(); mit++)
    {

      ModelEnt *me = mit->first;
      
    
      //---------------------------------//
      //            FACET                //
      //---------------------------------//
      facet(me);

      // Not sure if this needs to be done yet
      //set_senses(me);

      me->set_meshed_state(COMPLETE_MESH);
    }
}

void CurveFacetMeshReader::facet(ModelEnt *curve)
{
  iGeom::EntityHandle h = curve->geom_handle();
  iMesh::EntitySetHandle sh = IBSH(curve->mesh_handle());

  //get the facets for this curve/ref_edge
  std::vector<double> pnts;
  std::vector<int> conn;
  iGeom::Error errval = mk_core()->igeom_instance()->getFacets(h,facet_tol,pnts,conn);
  if (errval != iBase_SUCCESS) throw Error(MK_FAILURE, "Unable get the facets of the curve.");
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

  //--------------------------------------------------//
  //check that the start and end vertices are within  //
  //the geom resabs of the facet verts                //
  //--------------------------------------------------//
     
  //get the curve vertices
  MEntVector end_verts;
  curve->get_adjacencies(0, end_verts);
  
  //get the entity handle of the start vertex
  std::vector<moab::EntityHandle> dum_handles;
  end_verts.front()->get_mesh(0,dum_handles);
  assert( 1 == dum_handles.size() );   // should only be one vertex in a vertex entity set
  iMesh::EntityHandle start_vtx = IBEH(dum_handles.front());

  dum_handles.clear();

  // do the same for the end vertex
  end_verts.back()->get_mesh(0, dum_handles);
  assert( 1 == dum_handles.size() );   // should only be one vertex in a vertex entity set
  iMesh::EntityHandle end_vtx = IBEH(dum_handles.front());

  dum_handles.clear();

  //std::cout << "Number of end_verts: " << end_verts.size() << std::endl;

  //special case for a point curve
  if(verts.size() < 2)
    {
      if( start_vtx != end_vtx || curve->measure() > geom_res)
	{
	  std::cout << "Warning: no facetting for curve " << curve->id() << std::endl;
	  return;
	}
      else
	{
	  mk_core()->imesh_instance()->addEntToSet(start_vtx, sh);
	  return;
	}
    }

  //check for a closed curve
  bool closed = (mvtx2mvtx_dist(verts.front(),verts.back()) < geom_res);

  // check that geometry and facetting agree that the curve is closed
  if ( closed != (start_vtx == end_vtx))
    std::cout << "Warning: topology and geometry inconsistent for possibly closed curve "
	      << curve->id() << std::endl;

  //check the proximity of the front and end vertices to the start and end points, respectively
  if(vtx2vtx_dist(end_verts.front()->geom_handle(), verts.front()) > geom_res ||
     vtx2vtx_dist(end_verts.back()->geom_handle(), verts.back()) > geom_res)
    {
      //try reversing the points
      //std::reverse(verts.begin(), verts.end());
      
      //check again, if this time it fails, give a warning
      if(vtx2vtx_dist(end_verts.front()->geom_handle(), verts.front()) > geom_res ||
	 vtx2vtx_dist(end_verts.back()->geom_handle(), verts.back()) > geom_res)
	{
	  std::cerr << "Warning: vertices not at the ends of the curve " << curve->id() << std::endl;
	}
    }


  //now replace start and end facet points with the curve's child vert(s)
  
  //capture the beginning and end vertex handles for deletion 
  iMesh::EntityHandle front_vert_to_del = verts.front();
  iMesh::EntityHandle back_vert_to_del = verts.back();
  
  verts.front() = start_vtx;
  verts.back() = end_vtx;
  
  //delete the old vertices
  mk_core()->imesh_instance()->deleteEnt(front_vert_to_del);
  mk_core()->imesh_instance()->deleteEnt(back_vert_to_del);
    
  //now create the edges with the vertices we've chosen

  //vector to contain the edges
  std::vector<iBase_EntityHandle> edges;
  //loop over vertices and create edges
  for (unsigned int j = 0; j < verts.size()-1; j++)
    {
      //create edges
      iBase_EntityHandle e; 
      mk_core()->imesh_instance()->createEnt(iMesh_LINE_SEGMENT, &verts[j], 2, e);
      edges.push_back(e);
    }

  //remove duplicate vertex if the curve is closed
  if(closed && start_vtx == end_vtx) verts.pop_back(); 

  //add vertices and edges to the entity set
  mk_core()->imesh_instance()->addEntArrToSet(&verts[0], verts.size(), sh);
  mk_core()->imesh_instance()->addEntArrToSet(&edges[0], edges.size(), sh);
}

void CurveFacetMeshReader::set_senses(ModelEnt *curve)
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
double CurveFacetMeshReader::vtx2vtx_dist(iGeom::EntityHandle vtx1, iMesh::EntityHandle vtx2)
{

  double x1,y1,z1;
  double x2,y2,z2;

  mk_core()->igeom_instance()->getVtxCoord(vtx1, x1, y1, z1);
  mk_core()->imesh_instance()->getVtxCoord(vtx2, x2, y2, z2);

  double dist = pow((x1-x2),2) + pow((y1-y2),2) + pow((z1-z2),2);

  return sqrt(dist);
}

double CurveFacetMeshReader::mvtx2mvtx_dist(iMesh::EntityHandle vtx1, iMesh::EntityHandle vtx2)
{

  double x1,y1,z1;
  double x2,y2,z2;

  mk_core()->imesh_instance()->getVtxCoord(vtx1, x1, y1, z1);
  mk_core()->imesh_instance()->getVtxCoord(vtx2, x2, y2, z2);

  double dist = pow((x1-x2),2) + pow((y1-y2),2) + pow((z1-z2),2);

  return sqrt(dist);
}

void CurveFacetMeshReader::set_mesh_params(double faceting_tolerance, double geom_resabs)
{
  //Assign the faceting values if they are passed into the function
  if(faceting_tolerance) facet_tol = faceting_tolerance; 
  if(geom_resabs) geom_res = geom_resabs; 
}

}
