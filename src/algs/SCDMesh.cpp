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
// Initialize the entity types for SCDMesh
  moab::EntityType SCDMesh_tps[] = {moab::MBVERTEX, moab::MBHEX, moab::MBMAXTYPE};
  const moab::EntityType* SCDMesh::output_types() 
    { return SCDMesh_tps; }

//---------------------------------------------------------------------------//
// Constructor for SCDMesh
  SCDMesh::SCDMesh(MKCore *mk_core, const MEntVector &me_vec)
    : MeshScheme(mk_core, me_vec)
  {
    boxIncrease = .0;
  }

//---------------------------------------------------------------------------//
// setup function
void SCDMesh::setup_this()
{
  // if cfMesh chosen, check the grid entries to make sure they're correct
  if (gridType == 0) {
    if (coarse_i != fine_i.size()) {
      throw Error(MK_INCOMPLETE_MESH_SPECIFICATION, 
      "Number of coarse i divisions not equal to the number of fine i divisions.");
    }
    if (coarse_j != fine_j.size()) {
      throw Error(MK_INCOMPLETE_MESH_SPECIFICATION, 
      "Number of coarse j divisions not equal to the number of fine j divisions.");
      }
    if (coarse_k != fine_k.size()) {
      throw Error(MK_INCOMPLETE_MESH_SPECIFICATION, 
      "Number of coarse k divisions not equal to the number of fine k divisions.");
    }
  }

  // do setup for the model entities
  for (MEntSelection::iterator mit = mentSelection.begin(); mit != mentSelection.end(); mit++) {
    ModelEnt *me = mit->first;

    // check if the entity is already meshed
    if (me->get_meshed_state() >= COMPLETE_MESH || me->mesh_intervals() > 0) continue;
  }
}

//---------------------------------------------------------------------------//
// Structured cartesian mesh generation execution
  void SCDMesh::execute_this()
  {
    for (MEntSelection::iterator mit = mentSelection.begin(); mit != mentSelection.end(); mit++) {
      ModelEnt *me = mit->first;


      /* Generate the bounding box */

      // cartesian case
      if (axisType == 0) {
        set_box_dimension();
        create_cart_edges();
      }


      /* Generate the mesh */

      // create mesh vertices
      create_vertices();

      // full representation case
      if (interfaceType == 0) {
        create_full_mesh();
      }

      // lighweight ScdInterface case
      if (interfaceType == 1) {
        create_light_mesh();
      }


      /* Finish */

      // commit the mesh when finished
      me->commit_mesh(mit->second, COMPLETE_MESH);
    }
  }

//---------------------------------------------------------------------------//
// fix the box dimension for the bounding box
void SCDMesh::set_box_dimension()
{
  // this should be modified to get the bounding box on the geometry
  // that is only associated with the current model entity instead
  // of the entire igeom instance
  gerr = mk_core()->igeom_instance()->getBoundBox(minCoord[0], minCoord[1], minCoord[2],
                                                  maxCoord[0], maxCoord[1], maxCoord[2]);
  IBERRCHK(gerr, "ERROR: couldn't get geometry bounding box");

  // increase box size for 3 directions 
  for (int i = 0; i < 3; i++) {
    double box_size = maxCoord[i] - minCoord[i];
    minCoord[i] -= box_size*boxIncrease; 
    maxCoord[i] += box_size*boxIncrease; 
  } 
}

//---------------------------------------------------------------------------//
// get the axis cartesian bounding box edges if cfmesh has been chosen
void SCDMesh::create_cart_edges()
{
  // generate the vector arrays for cartesian grid.
  i_arr.push_back(minCoord[0]);
  double icrs_size = (maxCoord[0] - minCoord[0]) / coarse_i;
  double ifn_size;
  unsigned int crsi;
  for (crsi = 0; crsi != coarse_i; crsi++) {
    ifn_size = icrs_size / fine_i[crsi];
    for (int fni = 1; fni != fine_i[crsi] + 1; fni++) {
      i_arr.push_back(minCoord[0] + crsi*icrs_size + fni*ifn_size);
    }
  }
  
  j_arr.push_back(minCoord[1]);
  double jcrs_size = (maxCoord[1] - minCoord[1]) / coarse_j;
  double jfn_size;
  unsigned int crsj;
  for (crsj = 0; crsj != coarse_j; crsj++) {
    jfn_size = jcrs_size / fine_j[crsj];
    for (int fnj = 1; fnj != fine_j[crsj] + 1; fnj++) {
      j_arr.push_back(minCoord[1] + crsj*jcrs_size + fnj*jfn_size);
    }
  }
  
  k_arr.push_back(minCoord[2]);
  double kcrs_size = (maxCoord[2] - minCoord[2]) / coarse_k;
  double kfn_size;
  unsigned int crsk;
  for (crsk = 0; crsk != coarse_k; crsk++) {
    kfn_size = kcrs_size / fine_k[crsk];
    for (int fnk = 1; fnk != fine_k[crsk] + 1; fnk++) {
      k_arr.push_back(minCoord[2] + crsk*kcrs_size + fnk*kfn_size);
    }
  }
}

//---------------------------------------------------------------------------//
// create the mesh vertices
  void SCDMesh::create_vertices()
  {
    num_i = i_arr.size();
    num_j = j_arr.size();
    num_k = k_arr.size();
    num_verts = num_i*num_j*num_k;

    full_coords.resize(3*num_verts);
    unsigned int idx;

    // generate the vertices
    for (int kval = 0; kval != num_k; kval++) {
      for (int jval = 0; jval != num_j; jval++) {
        for (int ival = 0; ival != num_i; ival++) {
          idx = ival + num_i*jval + num_i*num_j*kval;
          full_coords[3*idx] = i_arr[ival];
          full_coords[3*idx + 1] = j_arr[jval];
          full_coords[3*idx + 2] = k_arr[kval];
        }
      }
    }

    rval = mk_core()->moab_instance()->create_vertices(&full_coords[0],
                                                       num_verts,
                                                       vtx_range);
    MBERRCHK(rval, mk_core()->moab_instance());
  }

//---------------------------------------------------------------------------//
// create the full mesh representation
  void SCDMesh::create_full_mesh()
  {
    // generate hex entities from the vertices
    moab::EntityHandle local_conn[8];
    unsigned int idx;
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

          moab::EntityHandle this_elem;
          rval = mk_core()->moab_instance()->create_element( moab::MBHEX,
                                                             local_conn,
                                                             8,
                                                             this_elem);
          MBERRCHK(rval, mk_core()->moab_instance());
        }
      }
    }
  }

//---------------------------------------------------------------------------//
// create mesh using light-weight MOAB ScdInterface
  void SCDMesh::create_light_mesh()
  {
    // create an instance of the ScdInterface
    rval = mk_core()->moab_instance()->query_interface(scdIface);
    MBERRCHK(rval, mk_core()->moab_instance());

    // create an instance of the ScdBox class
    // for now the coordinate indexing starts at 0
    // this can be changed if necessary
    moab::ScdBox *scd_box;
    rval = scdIface->construct_box(moab::HomCoord(0, 0, 0, 1),
                                   moab::HomCoord(num_i-1, num_j-1, num_k-1, 1),
                                   &full_coords[0], 
                                   num_verts,
                                   scd_box);
    if (moab::MB_SUCCESS != rval) mk_core()->moab_instance()->release_interface(scdIface);
    MBERRCHK(rval, mk_core()->moab_instance());

    // get rid of ScdInterface once we are done
    rval = mk_core()->moab_instance()->release_interface(scdIface);
    MBERRCHK(rval, mk_core()->moab_instance());
  }

//---------------------------------------------------------------------------//

} // end namespace MeshKit

//---------------------------------------------------------------------------//
// end algs/SCDMesh.cpp 
//---------------------------------------------------------------------------//
