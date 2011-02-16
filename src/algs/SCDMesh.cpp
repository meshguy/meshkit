//-----------------------------------C++-------------------------------------//
// File: algs/SCDMesh.cpp
// Author: Stuart R. Slattery
// Wednesday February 2 16:15:16 2011
// Brief: SCDMesh member definitions
//---------------------------------------------------------------------------//

#include "meshkit/SCDMesh.hpp"
#include "meshkit/MKCore.hpp"
#include "meshkit/ModelEnt.hpp"
#include "meshkit/SizingFunction.hpp"
#include "meshkit/RegisterMeshOp.hpp"

#include <vector>

#include "moab/Core.hpp"
#include "moab/Range.hpp"
#include "moab/Interface.hpp"
#include "moab/ScdInterface.hpp"

namespace MeshKit
{
//---------------------------------------------------------------------------//
// Static registration of this mesh scheme   
  moab::EntityType SCDMesh_tps[] = {moab::MBVERTEX, moab::MBHEX, moab::MBMAXTYPE};
  const moab::EntityType* SCDMesh::output_types() 
    { return SCDMesh_tps; }

//---------------------------------------------------------------------------//
// setup function
void SCDMesh::setup_this()
{
  // check the grid entries to make sure they're correct
  if (gridType == 0) {
    if (coarse_i != fine_i.size()) {
      throw Error(MK_INCOMPLETE_MESH_SPECIFICATION, "Number of coarse i divisions not equal to the number of fine i divisions.");
    }
    if (coarse_j != fine_j.size()) {
      throw Error(MK_INCOMPLETE_MESH_SPECIFICATION, "Number of coarse j divisions not equal to the number of fine j divisions.");
      }
    if (coarse_k != fine_k.size()) {
      throw Error(MK_INCOMPLETE_MESH_SPECIFICATION, "Number of coarse k divisions not equal to the number of fine k divisions.");
    }
  }
}

//---------------------------------------------------------------------------//
// Structured cartesian mesh generation execution
  void SCDMesh::execute_this()
  {
      if (axisType == 0) {
        create_cart_box();
      }

      /* Generate the mesh */
 
      // full representation case //
      if (interfaceType == 0) {
        create_full_mesh();
      }
  }

//---------------------------------------------------------------------------//
// export the mesh to a file
  void SCDMesh::export_mesh(const char* file_name)
  {
    rval = mk_core()->moab_instance()->write_mesh(file_name);
    MBERRCHK(rval, mk_core()->moab_instance());
  }


void SCDMesh::set_box_dimension()
{
  gerr = mk_core()->igeom_instance()->getBoundBox(minCoord[0], minCoord[1], minCoord[2],
                                                  maxCoord[0], maxCoord[1], maxCoord[2]);
  IBERRCHK(gerr, "ERROR: couldn't get geometry bounding box");

  // increase box size for 3 directions
  for (int i = 0; i < 3; i++) {
    minCoord[i] -= minCoord[i]*boxIncrease;
    maxCoord[i] += maxCoord[i]*boxIncrease;
  }
}

//---------------------------------------------------------------------------//
// get the axis cartesian bounding box
// for now this just assumes that gridType has been set to cfmesh
void SCDMesh::create_cart_box()
{
  // generate the vector arrays for cartesian grid.
  // at some point this should be generalized to make 
  // the grid for a axis aligned box from and obb tree too
  i_arr.push_back(minCoord[0]);
  double icrs_size = (maxCoord[0] - minCoord[0]) / coarse_i;
  double ifn_size;
  for (int crsi = 0; crsi != coarse_i; crsi++) {
    ifn_size = icrs_size / fine_i[crsi];
    for (int fni = 1; fni != fine_i[crsi] + 1; fni++) {
      i_arr.push_back(minCoord[0] + crsi*icrs_size + fni*ifn_size);
    }
  }
  
  j_arr.push_back(minCoord[1]);
  double jcrs_size = (maxCoord[1] - minCoord[1]) / coarse_j;
  double jfn_size;
  for (int crsj = 0; crsj != coarse_j; crsj++) {
    jfn_size = jcrs_size / fine_j[crsj];
    for (int fnj = 1; fnj != fine_j[crsj] + 1; fnj++) {
      j_arr.push_back(minCoord[1] + crsj*jcrs_size + fnj*jfn_size);
    }
  }
  
  k_arr.push_back(minCoord[2]);
  double kcrs_size = (maxCoord[2] - minCoord[2]) / coarse_k;
  double kfn_size;
  for (int crsk = 0; crsk != coarse_k; crsk++) {
    kfn_size = kcrs_size / fine_k[crsk];
    for (int fnk = 1; fnk != fine_k[crsk] + 1; fnk++) {
      k_arr.push_back(minCoord[2] + crsk*kcrs_size + fnk*kfn_size);
    }
  }
}

//---------------------------------------------------------------------------//
// create the full mesh representaton
  void SCDMesh::create_full_mesh()
  {
    int num_i = i_arr.size();
    int num_j = j_arr.size();
    int num_k = k_arr.size();
    int num_verts = num_i*num_j*num_k;

    // generate the vertices
    std::vector<double> full_coords;
    for (int kval = 0; kval != num_k; kval++) {
      for (int jval = 0; jval != num_j; jval++) {
        for (int ival = 0; ival != num_i; ival++) {
          full_coords.push_back(i_arr[ival]);
          full_coords.push_back(j_arr[jval]);
          full_coords.push_back(k_arr[kval]);
        }
      }
    }

    moab::Range vtx_range;
    rval = mk_core()->moab_instance()->create_vertices(&full_coords[0],
                                                       num_verts,
                                                       vtx_range);
    MBERRCHK(rval, mk_core()->moab_instance());

    // generate hex entities from the vertices
    std::vector<moab::EntityHandle* > conn;
    moab::EntityHandle local_conn[8];
    int idx;
    int num_hex = (num_i - 1)*(num_j - 1)*(num_k - 1);
    for (int kv = 0; kv != num_k - 1; kv++) {
      for (int jv = 0; jv != num_j - 1; jv++) {
        for (int iv = 0; iv != num_i - 1; iv++) {
          idx = iv + jv*num_i + kv*num_i*num_j;
          local_conn[0] = vtx_range[ idx];
          local_conn[1] = vtx_range[ idx + 1];
          local_conn[2] = vtx_range[ idx + 1 + num_i];
          local_conn[3] = vtx_range[ idx +     num_i];
          local_conn[4] = vtx_range[ idx +             num_i*num_j];
          local_conn[5] = vtx_range[ idx + 1 +         num_i*num_j];
          local_conn[6] = vtx_range[ idx + 1 + num_i + num_i*num_j];
          local_conn[7] = vtx_range[ idx +     num_i + num_i*num_j];
          conn.push_back(local_conn);
        }
      }
    }

    moab::EntityHandle elements[num_hex];
    for (unsigned hx = 0; hx != num_hex; hx++) {
      rval = mk_core()->moab_instance()->create_element( moab::MBHEX,
                                                         conn[hx],
                                                         8,
                                                         elements[hx]);
      MBERRCHK(rval, mk_core()->moab_instance());
    }
  }

//---------------------------------------------------------------------------//

} // end namespace MeshKit

//---------------------------------------------------------------------------//
// end algs/SCDMesh.cpp 
//---------------------------------------------------------------------------//
