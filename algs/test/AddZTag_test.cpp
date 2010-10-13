#include <iostream>
#include "iMesh.h"
#include "vec_utils.hpp"
#include "string.h"
#include "stdlib.h"
#define ERROR(a) {if (iBase_SUCCESS != err) std::cerr << a << std::endl;}
#define ERRORR(a,b) {if (iBase_SUCCESS != err) {std::cerr << a << std::endl; return b;}}

// add a Z tag to a mesh read from a file, and save it
int main(int argc, char* argv[] )
{
  // check command line arg
  const char *filename_mesh = 0;
  
  const char *newFile = 0;
  if (argc == 3)
    {
      filename_mesh = argv[1];
      newFile = argv[2];
    }
  else {
    std:: cout<< "Usage : addZTag <smf_file> <outvtk> \n";
    return 0;
  }
  int err;
  iMesh_Instance mesh;
  iMesh_newMesh("", &mesh, &err, 0);
  ERRORR("Couldn't create mesh.", 1);

  iBase_EntitySetHandle root_set;
  iMesh_getRootSet(mesh, &root_set, &err);
  ERRORR("Couldn't get root set.", 1);

  // read  mesh
 
  iMesh_load(mesh, root_set, filename_mesh, NULL, &err, strlen(filename_mesh), 0);
  ERRORR("Couldn't load mesh file.", 1);
  
  // read coordinates of nodes, and add a tag equal to the z of each node
 
  iBase_EntityHandle *verts = NULL;
  int verts_alloc = 0;
  int numNodes = 0;
  iMesh_getEntities(mesh, root_set, iBase_VERTEX, 
                    iMesh_POINT, &verts, &verts_alloc, &numNodes, &err);
  ERRORR("failed to get vertices.", 1);
  

  /* get the coordinates in one array */
   
  int vert_coords_alloc = 0;
  double * xyz = 0; // not allocated 
  int vertex_coord_size = 0;

  iMesh_getVtxArrCoords(mesh, verts, numNodes, iBase_INTERLEAVED, 
                        &xyz, &vert_coords_alloc, &vertex_coord_size, &err);
  ERRORR("failed to get vertex coordinates of entities in getMeshData.", 1);

  iBase_TagHandle elev_tag_handle;
  const char * tagName1 = "ZTag";
  iMesh_createTag       ( mesh,
                        tagName1,
                        /*  size ? */ 1,
                        iBase_DOUBLE ,
                        &elev_tag_handle,
                        &err,
                        strlen(tagName1) ) ;
  ERRORR("failed to create tag.", 1);
  double * dArr = new double[numNodes];
  int j=0;
  for (j=0; j<numNodes; j++)
   {
     dArr[j] = xyz[j*3+2];
   }    

  iMesh_setDblArrData(mesh,
       /*in const iBase_EntityHandle* */ verts,
       /*in const int */  numNodes,
       /*in const iBase_TagHandle*/ elev_tag_handle,
       /*in const double* */  dArr,
       /*in const int */ numNodes,
       /*out int * */ &err);

   ERRORR("failed to set tag values.", 1);

  // get the triangles and the vertices in one shot
  iBase_EntityHandle *triangles = NULL;
  int triangles_alloc = 0;
  iBase_EntityHandle *vert_adj = NULL;
  int  vert_adj_alloc = 0, vert_adj_size;
  int * offsets = NULL, offsets_alloc = 0, indices_size;
  int * indices = NULL, indices_alloc = 0, offsets_size;
  int numTriangles;
  iMesh_getAdjEntIndices( mesh, root_set, 
			  iBase_FACE, iMesh_TRIANGLE, iBase_VERTEX,
			  &triangles, &triangles_alloc, &numTriangles,
			  &vert_adj, &vert_adj_alloc, &vert_adj_size,
			  &indices, &indices_alloc, &indices_size,
			  &offsets, &offsets_alloc, &offsets_size,
			  &err );
  ERRORR("failed to get triangles and vertices.", 1);

  /* get the coordinates in one array */
   
  vert_coords_alloc = 0;
  free (xyz);
  xyz = NULL;

  iMesh_getVtxArrCoords(mesh, vert_adj, vert_adj_size, iBase_INTERLEAVED, 
                        &xyz, &vert_coords_alloc, &vertex_coord_size, &err);
  ERRORR("failed to get vertex coordinates.", 1);

  iBase_TagHandle orientation_tag_handle;
  const char * tagName = "Orientation";
  iMesh_createTag       ( mesh,
                        tagName,
                        /*  size ? */ 1,
                        iBase_INTEGER ,
                        &orientation_tag_handle,
                        &err,
                        strlen(tagName) ) ;
  ERRORR("failed to create orientation tag.", 1);
  int * orientation = new int[numTriangles];
  
  iBase_TagHandle normal_tag_handle;
  const char * tagNameNormal = "normalTag";
  iMesh_createTag       ( mesh,
                        tagNameNormal,
                        /*  size ? */ 1,
                        iBase_DOUBLE ,
                        &normal_tag_handle,
                        &err,
                        strlen(tagNameNormal) ) ;
  ERRORR("failed to create normal tag.", 1);
  double * normtag = new double[numTriangles];
  
  // now, a loop over each triangle, to compute the normal
  // start copy
  int numberNegative = 0;
  int numberZeroArea = 0;
  for (j=0; j<numTriangles; j++)
    {
      int v[3];
      for (int ii=0; ii<3; ii++)
         v[ii]=indices[offsets[j]+ii];
      // now, the triangle has coordinates x, y of interest
      double area = (xyz[3*v[1]] - xyz[3*v[0]]) * (xyz[3*v[2]+1] - xyz[3*v[0]+1]) -
		     (xyz[3*v[2]] - xyz[3*v[0]]) * (xyz[3*v[1]+1] - xyz[3*v[0]+1]);
      orientation[j] = 0;
      if (area<0)
      {
	orientation[j] = -1;
        numberNegative++;
      }
      else if (area >0)
	orientation[j] =1;
      else 
	numberZeroArea++;
      // repeat some calculations here :((
      double normal3[3];
      normal3D(normal3, &(xyz[3*v[0]]), &(xyz[3*v[1]]), &(xyz[3*v[2]]));
      normtag[j] = normal3[2];
    }

  // positive or negative triangles
  iMesh_setIntArrData  	(  mesh,
			   triangles,
			   numTriangles,
			   orientation_tag_handle,
			   orientation,
			   numTriangles,
			   &err); 
  ERRORR("failed to set orientation tag.", 1);

  iMesh_setDblArrData(mesh,
       /*in const iBase_EntityHandle* */ triangles,
       /*in const int */  numTriangles,
       /*in const iBase_TagHandle*/ normal_tag_handle,
       /*in const double* */  normtag,
       /*in const int */ numTriangles,
       /*out int * */ &err);

   ERRORR("failed to set normal tag values.", 1);


  ///
  iMesh_save(mesh, root_set, newFile, NULL, &err, strlen(newFile), 0);
  ERRORR("Couldn't save mesh.", 1);
  
  iMesh_dtor(mesh, &err);

  ERRORR("Couldn't destroy mesh.", 1);
  std::cout<<"numberNegative triangles : " << numberNegative << " number zero area triangles "<< numberZeroArea << std::endl;
  delete [] dArr;
  delete [] normtag;
  delete [] orientation;
  free (xyz);
  free (verts);
  return 0;
}
