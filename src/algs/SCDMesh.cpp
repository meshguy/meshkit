//-----------------------------------C++-------------------------------------//
// File: algs/SCDMesh.cpp
// Author: Stuart R. Slattery
// Wednesday February 2 16:15:16 2011
// Brief: SCDMesh member definitions
//---------------------------------------------------------------------------//

#include "meshkit/SCDMesh.hpp"
#include "meshkit/MKCore.hpp"
#include "meshkit/ModelEnt.hpp"
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
    boxIncrease = 0.0;

    // set bounding box size tag
    double bb_default[6] = { 0., 0., 0., 0., 0., 0. };
    rval = mk_core->moab_instance()->tag_create("BOUNDING_BOX_SIZE",
                                                6*sizeof(double), moab::MB_TAG_SPARSE,
                                                moab::MB_TYPE_DOUBLE, bb_tag, bb_default, true);
    if (moab::MB_SUCCESS != rval && moab::MB_ALREADY_ALLOCATED != rval) 
      MBERRCHK(rval, mk_core->moab_instance());
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

      // Cartesian case
      if (axisType == 0) {
        if (geomType == 0) {
          set_cart_box_all();
        }
        else if (geomType == 1) {
          set_cart_box_individual(me);
        }
        create_cart_edges();
      }

      /* Generate the mesh */

      // create mesh vertex coordinates
      create_vertex_coords();

      // full representation case
      if (interfaceType == 0) {
        create_full_mesh(mit->second);
      }

      // lighweight ScdInterface case
      if (interfaceType == 1) {
        create_light_mesh(mit->second);
      }

      // stop if mesh is created for whole geometry
      if (geomType == 0) break;

      // commit the mesh when finished
      me->commit_mesh(mit->second, COMPLETE_MESH);
    }
  }

//---------------------------------------------------------------------------//
// set the Cartesian bounding box dimension for the entire geometry
  void SCDMesh::set_cart_box_all()
  {
    gerr = mk_core()->igeom_instance()->getBoundBox(minCoord[0], minCoord[1], minCoord[2],
                                                  maxCoord[0], maxCoord[1], maxCoord[2]);
    IBERRCHK(gerr, "ERROR: couldn't get geometry bounding box");

    // increase box size
    double box_size;
    for (int i = 0; i < 3; i++) {
      box_size = maxCoord[i] - minCoord[i];
      minCoord[i] -= box_size*boxIncrease; 
      maxCoord[i] += box_size*boxIncrease; 
    }
  }

//---------------------------------------------------------------------------//
// set the Cartesian bounding box dimension for an individual volume
  void SCDMesh::set_cart_box_individual(ModelEnt *this_me)
  {
    gerr = this_me->igeom_instance()->getBoundBox(minCoord[0], minCoord[1], minCoord[2],
                                                  maxCoord[0], maxCoord[1], maxCoord[2]);
    IBERRCHK(gerr, "ERROR: couldn't get geometry bounding box");

    // increase box size
    double box_size;
    for (int i = 0; i < 3; i++) {
      box_size = maxCoord[i] - minCoord[i];
      minCoord[i] -= box_size*boxIncrease; 
      maxCoord[i] += box_size*boxIncrease; 
    }

    // set bounding box size as tag
    moab::EntityHandle meshset = this_me->mesh_handle();
    double box_min_max[6] = { minCoord[0], minCoord[1], minCoord[2],
                              maxCoord[0], maxCoord[1], maxCoord[2] };
    rval = mk_core()->moab_instance()->tag_set_data(bb_tag, &meshset, 1,
                                                    box_min_max);
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
// create the mesh vertex coordinates
  void SCDMesh::create_vertex_coords()
  {
    num_i = i_arr.size();
    num_j = j_arr.size();
    num_k = k_arr.size();
    num_verts = num_i*num_j*num_k;

    full_coords.resize(3*num_verts);
    unsigned int idx;

    // generate the vertex coordinate array in an interleaved fashion
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
  }

//---------------------------------------------------------------------------//
// create the full mesh representation
  void SCDMesh::create_full_mesh(moab::Range& me_range)
  {
    // create the vertices from the coordinates array
    moab::Range vtx_range;
    rval = mk_core()->moab_instance()->create_vertices(&full_coords[0],
                                                       num_verts,
                                                       vtx_range);
    MBERRCHK(rval, mk_core()->moab_instance());

    if (geomType == 0) me_range.merge(vtx_range);

    // generate hex entities from the vertices
    // the entity handles here should be contiguous with the i blocks 
    // moving the fastest and the k blocks moving the slowest
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
          if (geomType == 0) me_range.insert(this_elem);
        }
      }
    }
  }

//---------------------------------------------------------------------------//
// create mesh using light-weight MOAB ScdInterface
  void SCDMesh::create_light_mesh(moab::Range& me_range)
  {
    // create an instance of the ScdInterface
    rval = mk_core()->moab_instance()->query_interface(scdIface);
    MBERRCHK(rval, mk_core()->moab_instance());

    // create an instance of the ScdBox class
    // for now the coordinate indexing starts at 0
    // this should be changed to account for multiple geometric instances being meshed
    moab::ScdBox *scd_box;
    rval = scdIface->construct_box(moab::HomCoord(0, 0, 0, 1),
                                   moab::HomCoord(num_i-1, num_j-1, num_k-1, 1),
                                   &full_coords[0], 
                                   num_verts,
                                   scd_box);
    MBERRCHK(rval, mk_core()->moab_instance());

    // get rid of ScdInterface once we are done
    rval = mk_core()->moab_instance()->release_interface(scdIface);
    MBERRCHK(rval, mk_core()->moab_instance());

    // add created vertex and hexes
    if (geomType == 1) {
      moab::EntityHandle start;
      start = scd_box->start_vertex();
      moab::Range vert_range(start, start + scd_box->num_vertices()-1);
      me_range.merge(vert_range);
      start = scd_box->start_element();
      moab::Range hex_range(start, start + scd_box->num_elements()-1);
      me_range.merge(hex_range);
    }
  }

//---------------------------------------------------------------------------//

} // end namespace MeshKit

//---------------------------------------------------------------------------//
// end algs/SCDMesh.cpp 
//---------------------------------------------------------------------------//
