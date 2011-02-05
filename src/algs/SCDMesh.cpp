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

#include <vector>

#include "moab/Core.hpp"
#include "moab/Range.hpp"
#include "moab/Interface.hpp"
#include "moab/ScdInterface.hpp"

namespace MeshKit
{
//---------------------------------------------------------------------------//
// Static registration of this mesh scheme   
  moab::EntityType SCDMesh_tps[] = {moab::MBVERTEX, moab::MBHEX};
  iBase_EntityType SCDMesh_mtp = iBase_REGION;

  static int SCD_success= MKCore::register_meshop("SCDMesh",
                                                  &SCDMesh_mtp,
                                                  1,
                                                  SCDMesh_tps,
                                                  2,
                                                  SCDMesh::factory,
                                                  MeshOp::canmesh_region);

//---------------------------------------------------------------------------//
// make an instance of the SCDMesh class
  MeshOp *SCDMesh::factory(MKCore *mkcore, const MEntVector &me_vec)
  {
   return new SCDMesh(mkcore, me_vec);
  }

//---------------------------------------------------------------------------//
// return the type of entities this mesh creates
  void SCDMesh::mesh_types(std::vector<moab::EntityType> &tps)
  {
    tps.push_back(moab::MBVERTEX);
    tps.push_back(moab::MBHEX);
  }

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

    // not sure yet how to implement the sizing function for this meshop
    // but it will go here

    // call generic setup function to finish
    setup_boundary();
  }

//---------------------------------------------------------------------------//
// Structured cartesian mesh generation execution
  void SCDMesh::execute_this()
  {

    for (MEntSelection::iterator mit = mentSelection.begin();
         mit != mentSelection.end(); mit++) {

      ModelEnt *me = mit->first;

      /* get the bounding box on the geometry */

      // bounding box parameters
      double *xmn, *ymn, *zmn, *xmx, *ymx, *zmx;

      // cartesian case //
      if (axisType == 0) {
        create_cart_box(me, xmn, ymn, zmn, xmx, ymx, zmx);
      }


      /* Generate the mesh */
 
      // full representation case //
      if (interfaceType == 0) {
        create_full_mesh(me);
      }


      /* Commit the mesh */
      me->commit_mesh(mit->second, COMPLETE_MESH);
    }
  }

//---------------------------------------------------------------------------//
// export the mesh to a file
  void SCDMesh::export_mesh(const char* file_name)
  {
    rval = mk_core()->moab_instance()->write_mesh(file_name);
    MBERRCHK(rval, "ERROR: couldn't write SCDMesh");
  }

//---------------------------------------------------------------------------//
// get the axis cartesian bounding box
// for now this just assumes that gridType has been set to cfmesh
  void SCDMesh::create_cart_box(ModelEnt *ent, 
                                double *xmin, double *ymin, double *zmin,
                                double *xmax, double *ymax, double *zmax)
  {
      gerr = mk_core()->igeom_instance()->getBoundBox(*xmin,
                                                      *ymin,
                                                      *zmin,
                                                      *xmax,
                                                      *ymax,
                                                      *zmax);
      IBERRCHK(gerr, "ERROR: couldn't get geometry bounding box");

      // generate the vector arrays for cartesian grid.
      // at some point this should be generalized to make 
      // the grid for a axis aligned box from and obb tree too
      i_arr.push_back(*xmin);
      double icrs_size = (*xmax - *xmin) / coarse_i;
      double ifn_size;
      for (int crsi = 0; crsi != coarse_i; crsi++) {
        ifn_size = icrs_size / fine_i[crsi];
        for (int fni = 1; fni != fine_i[crsi] + 1; fni++) {
          i_arr.push_back(*xmin + crsi*icrs_size + fni*ifn_size);
        }
      }

      j_arr.push_back(*ymin);
      double jcrs_size = (*ymax - *ymin) / coarse_j;
      double jfn_size;
      for (int crsj = 0; crsj != coarse_j; crsj++) {
        jfn_size = jcrs_size / fine_j[crsj];
        for (int fnj = 1; fnj != fine_j[crsj] + 1; fnj++) {
          j_arr.push_back(*ymin + crsj*jcrs_size + fnj*jfn_size);
        }
      }

      k_arr.push_back(*zmin);
      double kcrs_size = (*zmax - *zmin) / coarse_k;
      double kfn_size;
      for (int crsk = 0; crsk != coarse_k; crsk++) {
        kfn_size = kcrs_size / fine_k[crsk];
        for (int fnk = 1; fnk != fine_k[crsk] + 1; fnk++) {
          k_arr.push_back(*zmin + crsk*kcrs_size + fnk*kfn_size);
        }
      }
  }

//---------------------------------------------------------------------------//
// create the full mesh representaton
  void SCDMesh::create_full_mesh(ModelEnt *ent)
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
    MBERRCHK(rval, "ERROR: couldn't create mesh nodes");

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
      MBERRCHK(rval, "ERROR: couldn't create hex elements");
    }
  }

//---------------------------------------------------------------------------//

} // end namespace MeshKit

//---------------------------------------------------------------------------//
// end algs/SCDMesh.cpp 
//---------------------------------------------------------------------------//
