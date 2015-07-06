

#include "meshkit/MKCore.hpp"
#include "meshkit/SurfaceFacetMeshReader.hpp"
#include "meshkit/CurveFacetMeshReader.hpp"
#include "meshkit/ModelEnt.hpp"
#include <iostream>
#include <iGeom.h>

#include "MBCore.hpp"


namespace MeshKit
{

// Output mesh types for this class
moab::EntityType SurfaceFacetMeshReader_types[] = {moab::MBVERTEX, moab::MBEDGE, moab::MBTRI, moab::MBMAXTYPE};
const moab::EntityType* SurfaceFacetMeshReader::output_types()
{ return SurfaceFacetMeshReader_types; }

SurfaceFacetMeshReader::SurfaceFacetMeshReader(MKCore *mk_core, const MEntVector &ments)
  : MeshScheme(mk_core, ments)
{
  mk = mk_core; 
}

SurfaceFacetMeshReader::~SurfaceFacetMeshReader()
{
}

void SurfaceFacetMeshReader::setup_this()
{
  //set the faceting parameters to default values if they are not already set
  
  //set the faceting tolerance
  if(!facet_tol) facet_tol = 1e-4; 
  
  //set the geom_resabs value
  if(!geom_res) geom_res = 1e-6;

  //create a solid curve mesher 
  CurveFacetMeshReader *cfmr = (CurveFacetMeshReader*) mk_core()->construct_meshop("CurveFacetMeshReader");

  //set the parameters of the curvemesher to match those of the surface mesher
  cfmr->set_mesh_params(facet_tol,geom_res);

  for(MEntSelection::iterator mit = mentSelection.begin(); mit != mentSelection.end(); mit++)
    {
      ModelEnt *me = mit->first;
      //do an initial check that the model entity is of the correct dimension
      if(me->dimension() != 2)
	{
	  std::cout << "Found a model entity with an incorrect dimension for this mesher" << std::endl;
          mentSelection.erase(mit); 
          continue; 
	}

      //make sure that the model entity's children have been meshed
      MEntVector children; 
      me->get_adjacencies(1, children); 
      //std::cout << "Number of children found: " << children.size() << std::endl;
 
      for(MEntVector::iterator j = children.begin(); j != children.end(); j++)
	{
	  ModelEnt *child = *j; 
          if(child->is_meshops_list_empty()) cfmr-> add_modelent(child); 
	}

    }

  //have the children be meshed first
  mk_core()->insert_node(cfmr, this, mk_core()->root_node());

}

void SurfaceFacetMeshReader::execute_this()
{

  for(MEntSelection::iterator mit = mentSelection.begin(); mit !=mentSelection.end(); mit++)
    {
      ModelEnt *me = mit->first; 
      
      //---------------------------------//
      //            FACET                //
      //---------------------------------//
      facet(me);

      me->set_meshed_state(COMPLETE_MESH);
    }
}

void SurfaceFacetMeshReader::facet(ModelEnt *surf)
{

  iGeom::EntityHandle h = surf->geom_handle(); 
  iMesh::EntitySetHandle sh = IBSH(surf->mesh_handle()); 

  //get the facets for this surface/ ref_face
  std::vector<double> pnts;
  std::vector<int> conn;
  iGeom::Error errval = mk->igeom_instance()->getFacets(h,facet_tol,pnts,conn);
  if (errval != iBase_SUCCESS) throw Error(MK_FAILURE, "Unable get the facets of the surface.");
  
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
      

  // get the geometric vertices of the surface and merge them into the mesh by 
  // a proximity check

  //get the geometric vertices
  MEntVector geom_verts;
  surf->get_adjacencies(0,geom_verts);
  unsigned int matches = 0;
 
  //loop over all the vertices we've created
  for(unsigned int j=0; j<verts.size() ; j++)
    {
      //loop over all geometric vertices adjacent to this surface
      for(unsigned int k=0; k<geom_verts.size(); k++)
	{
          // if thhe current vert is within the geom_res proximity of a geometric vert
          if(vtx2vtx_dist(geom_verts[k]->geom_handle(),verts[j]) < geom_res)
	    {
              //replace the vertex with the geometric vertex

	      //capture vert to remove from the list
	      iMesh::EntityHandle vert_to_del = verts[j];
                
              //replace the vertex with the geom vertex in the vector
	      std::vector<moab::EntityHandle> dum_handle;
	      geom_verts[k]->get_mesh(0, dum_handle);
              verts[j] = IBEH(dum_handle[0]);

              //delete the former vertex
              mk_core()->imesh_instance()->deleteEnt(vert_to_del);

              matches++;
	    }
	}
    }

  //now check that we have matched the correct number of vertices
  if (matches != geom_verts.size())
    {
      if(matches > geom_verts.size()) std::cout << "Warning: coincident vertices in surface " << surf->id() << std::endl;
      if(matches < geom_verts.size()) std::cout << "Warning: one or more geom vertices could not be matched in surface " << surf->id() << std::endl;
    }

  //vector for keeping track of the triangles
  std::vector<iBase_EntityHandle> tris;

  //loop over the connectivity
  for(int j = 0; j < (int)conn.size()-2; j+=3)
    {
      //get the appropriate points for a triangle and add them to a vector
      std::vector<iBase_EntityHandle> tri_verts; 
      tri_verts.push_back(verts[conn[j]]);
      tri_verts.push_back(verts[conn[j+1]]);
      tri_verts.push_back(verts[conn[j+2]]);

      //create a new triangle handle
      iBase_EntityHandle t; 
      mk->imesh_instance()->createEnt(iMesh_TRIANGLE, &tri_verts[0], 3, t);
      tri_verts.clear();
      tris.push_back(t);
    }
  //std::cout << "Created " << tris.size() << " triangles" << std::endl;

  //add verticess and edges to the entity set
  mk->imesh_instance()->addEntArrToSet(&verts[0], verts.size(), sh);
  mk->imesh_instance()->addEntArrToSet(&tris[0], tris.size(), sh);

}

void SurfaceFacetMeshReader::set_mesh_params(double faceting_tolerance, double geom_resabs)
{
  //Assign the faceting values if they are passed into the function
  if(faceting_tolerance) facet_tol = faceting_tolerance; 
  if(geom_resabs) geom_res = geom_resabs; 
}

double SurfaceFacetMeshReader::vtx2vtx_dist(iGeom::EntityHandle vtx1, iMesh::EntityHandle vtx2)
{

  double x1,y1,z1;
  double x2,y2,z2;

  mk_core()->igeom_instance()->getVtxCoord(vtx1, x1, y1, z1);
  mk_core()->imesh_instance()->getVtxCoord(vtx2, x2, y2, z2);

  double dist = pow((x1-x2),2) + pow((y1-y2),2) + pow((z1-z2),2);

  return sqrt(dist);
}

}
