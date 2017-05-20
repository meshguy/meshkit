
#include "MKVersion.h"
#include "meshkit/MeshImprove.hpp"
#include <iostream>
#include <math.h>
#include <map>

//#include "SweepWrapper.hpp"
#include "moab/mesquite/UntangleWrapper.hpp"
#include "moab/mesquite/ShapeImprover.hpp"
#include "moab/mesquite/SmartLaplacianSmoother.hpp"
#include "SweepWrapper.hpp"
#include "SmartLaplaceWrapper.hpp"
#include "moab/mesquite/ShapeImprovementWrapper.hpp"

#include "moab/mesquite/IdealWeightInverseMeanRatio.hpp"
#include "moab/mesquite/QualityAssessor.hpp"
#include "moab/mesquite/InstructionQueue.hpp"
#include "moab/mesquite/TerminationCriterion.hpp"

#include "moab/mesquite/LaplaceWrapper.hpp"
#include "moab/mesquite/SizeAdaptShapeWrapper.hpp"
#include "moab/mesquite/PaverMinEdgeLengthWrapper.hpp"
#include "moab/mesquite/DeformingDomainWrapper.hpp"
#include <moab/mesquite/MsqError.hpp>
#include <moab/mesquite/ShapeImprovementWrapper.hpp>
#include "moab/mesquite/MsqIMesh.hpp"
#include "moab/mesquite/MsqIGeom.hpp"
#ifdef HAVE_FBIGEOM
#include "meshkit/MsqFBiGeom.hpp"
#endif
using namespace MESQUITE_NS;

namespace MeshKit {

MeshImprove::MeshImprove(MKCore* core, bool isLaplacian, bool isUntangle,
    bool isShapeImprove, bool isSizeAdapt, iGeom * ig_inst)
{
  mk_core = core;
  IsLaplacian = isLaplacian;
  IsUntangle = isUntangle;
  IsShapeImprove = isShapeImprove;
  IsSizeAdapt = isSizeAdapt;
  if (ig_inst)
    igeom_inst = ig_inst;
  else
    igeom_inst = mk_core->igeom_instance();// it will pick the first one
}

void MeshImprove::SurfMeshImprove(iBase_EntityHandle surface,
    iBase_EntitySetHandle surfMesh, iBase_EntityType entity_type)
{

  MsqError mError;
  const char* VERTEX_FIXED_TAG_NAME = "MesquiteVertexFixed";

  iBase_TagHandle fixed_tag = 0;
  iMesh::Error m_err = mk_core->imesh_instance()->getTagHandle(
      VERTEX_FIXED_TAG_NAME, fixed_tag);
  if (m_err) {
    m_err = mk_core->imesh_instance()->createTag(VERTEX_FIXED_TAG_NAME, 1,
        iBase_INTEGER, fixed_tag);
    IBERRCHK(m_err, "Trouble create the tag handle.");
  }

  MsqIMesh mesh_adapter(mk_core->imesh_instance()->instance(), surfMesh,
      entity_type, mError, &fixed_tag);
  cout << "error =" << mError << endl;
  if (mError)
    throw mError;

  //get all the vertices in surface mesh: entity_handles_out---quads     adj_entity_handles_out---vertices
  std::vector<iBase_EntityHandle> entity_handles_out, adj_entity_handles_out;
  std::vector<int> offsets_out, adj_entity_indices_out;

  m_err = mk_core->imesh_instance()->getAdjEntIndices(surfMesh, entity_type,
      iMesh_ALL_TOPOLOGIES, iBase_VERTEX, entity_handles_out,
      adj_entity_handles_out, adj_entity_indices_out, offsets_out);
  IBERRCHK(m_err, "Trouble get the adjacent entity indices.");

  cout << "number of faces is " << entity_handles_out.size() << endl;
  cout << "number of vertices is " << adj_entity_handles_out.size() << endl;

  //set fixed flag on all vertices
  std::vector<int> tag_data(adj_entity_handles_out.size(), 1);
  m_err = mk_core->imesh_instance()->setIntArrData(&adj_entity_handles_out[0],
      adj_entity_handles_out.size(), fixed_tag, &tag_data[0]);
  IBERRCHK(m_err, "Trouble set an array of int data for a list of vertices.");

  //clear fixed flag for vertices contained directly in set
  int count = -1;
  m_err
      = mk_core->imesh_instance()->getNumOfType(surfMesh, iBase_VERTEX, count);
  IBERRCHK(m_err, "Trouble get the number of vertices in the set.");

  adj_entity_handles_out.clear();

  cout << "Num of Vertices on the target surface is " << count << endl;

  m_err = mk_core->imesh_instance()->getEntities(surfMesh, iBase_VERTEX,
      iMesh_ALL_TOPOLOGIES, adj_entity_handles_out);
  IBERRCHK(m_err, "Trouble get the nodes from the mesh entity set.");

  tag_data.clear();
  tag_data.resize(adj_entity_handles_out.size(), 0);

  m_err = mk_core->imesh_instance()->setIntArrData(&adj_entity_handles_out[0],
      adj_entity_handles_out.size(), fixed_tag, &tag_data[0]);
  IBERRCHK(m_err, "Trouble set an array of int data for mesh nodes.");

  //Finally smooth the mesh

  //SweepWrapper smoother( 1e-6,  "COORDINATES_MAP");
  if (IsUntangle) {
    UntangleWrapper smoother;
    //smoother.set_cpu_time_limit(300);
    smoother.set_vertex_movement_limit_factor(0.001);

    if (surface) {
#ifdef HAVE_FBIGEOM
      if (this->igeom_inst->isFBiGeom())
      {
        MsqFBiGeom fbgeom_adapter((FBiGeom *)igeom_inst, surface);
        MeshDomainAssoc mesh_and_domain(&mesh_adapter, &fbgeom_adapter);
        smoother.run_instructions(&mesh_and_domain, mError);
        cout << "Mesquite error in fb surface mesh smoothing=" << mError << endl;
      }
      else
      {
#endif
      MsqIGeom geom_adapter(mk_core->igeom_instance()->instance(), surface);
      MeshDomainAssoc mesh_and_domain(&mesh_adapter, &geom_adapter);
      smoother.run_instructions(&mesh_and_domain, mError);
      cout << "Mesquite error in surface mesh smoothing=" << mError << endl;
#ifdef HAVE_FBIGEOM
    }
#endif
    } else {
      smoother.run_instructions(&mesh_adapter, mError);
      cout << "Mesquite error in surface mesh smoothing=" << mError << endl;
    }
  }

  //use the SmartLaplaceWrapper class the smooth the target surface mesh
  if (IsLaplacian) {
    LaplaceWrapper sl_smooth;
    TerminationCriterion terminate;
    terminate.write_iterations("mesquite.gpt", mError);
    if (surface) {
#ifdef HAVE_FBIGEOM
      if (this->igeom_inst->isFBiGeom())
      {
        MsqFBiGeom fbgeom_adapter((FBiGeom *)igeom_inst, surface);
        MeshDomainAssoc mesh_and_domain(&mesh_adapter, &fbgeom_adapter);
        sl_smooth.run_instructions(&mesh_and_domain, mError);
        cout << "Mesquite error in fb smart Laplacian surface mesh smoothing with the geometry domain=" << mError << endl;
      }
      else
      {
#endif
      MsqIGeom geom_adapter(mk_core->igeom_instance()->instance(), surface);
      MeshDomainAssoc mesh_and_domain(&mesh_adapter, &geom_adapter);
      sl_smooth.run_instructions(&mesh_and_domain, mError);
      cout
          << "Mesquite error in the smart Laplacian surface mesh smoothing with the geometry domain="
          << mError << endl;
#ifdef HAVE_FBIGEOM
    }
#endif
    } else {
      MeshDomainAssoc mesh_and_domain(&mesh_adapter, 0);
      sl_smooth.run_instructions(&mesh_and_domain, mError);
      cout
          << "Mesquite error in the smart Laplacian surface mesh smoothing without the geometry domain="
          << mError << endl;
    }
  }

  //use the ShapeImprover class to smooth the target surface mesh
  if (IsShapeImprove) {
    IdealWeightInverseMeanRatio extra_metric;
    ShapeImprovementWrapper smoother1;
    smoother1.quality_assessor().add_quality_assessment(&extra_metric);
    //smoother1.set_vertex_movement_limit_factor(0.01);
    //smoother1.set_cpu_time_limit(300);
    if (surface) {
#ifdef HAVE_FBIGEOM
      if (this->igeom_inst->isFBiGeom())
      {
        MsqFBiGeom fbgeom_adapter((FBiGeom *)igeom_inst, surface);
        MeshDomainAssoc mesh_and_domain(&mesh_adapter, &fbgeom_adapter);
        smoother1.run_instructions(&mesh_and_domain, mError);
        cout << "Mesquite error in ShapeImprover fb surface mesh smoothing=" << mError << endl;
      }
      else
      {
#endif
      MsqIGeom geom_adapter(mk_core->igeom_instance()->instance(), surface);
      MeshDomainAssoc mesh_and_domain(&mesh_adapter, &geom_adapter);
      smoother1.run_instructions(&mesh_and_domain, mError);
      cout << "Mesquite error in the ShapeImprover surface mesh smoothing="
          << mError << endl;
#ifdef HAVE_FBIGEOM
    }
#endif
    } else {
      MeshDomainAssoc mesh_and_domain(&mesh_adapter, 0);
      smoother1.run_instructions(&mesh_and_domain, mError);
      cout << "Mesquite error in the ShapeImprover surface mesh smoothing="
          << mError << endl;
    }
  }

  //Use the Minimum Edge-Length Improvement

  if (IsSizeAdapt) {
    SizeAdaptShapeWrapper Smooth_SizeAdapt(1.0e-2);
    if (surface) {
#ifdef HAVE_FBIGEOM
      if (this->igeom_inst->isFBiGeom())
      {
        MsqFBiGeom fbgeom_adapter((FBiGeom *)igeom_inst, surface);
        MeshDomainAssoc mesh_and_domain(&mesh_adapter, &fbgeom_adapter);
        Smooth_SizeAdapt.run_instructions(&mesh_and_domain, mError);
        cout << "Mesquite error in SizeAdaptShape fb surface mesh smoothing=" << mError << endl;
      }
      else
      {
#endif
      MsqIGeom geom_adapter(mk_core->igeom_instance()->instance(), surface);
      MeshDomainAssoc mesh_and_domain(&mesh_adapter, &geom_adapter);
      Smooth_SizeAdapt.run_instructions(&mesh_and_domain, mError);
      cout << "Mesquite error in the SizeAdaptShape surface mesh smoothing="
          << mError << endl;
#ifdef HAVE_FBIGEOM
    }
#endif

    } else {
      Smooth_SizeAdapt.run_instructions(&mesh_adapter, mError);
      cout << "Mesquite error in the SizeAdaptShape surface mesh smoothing="
          << mError << endl;
    }
  }

}

void MeshImprove::VolumeMeshImprove(iBase_EntitySetHandle volMesh,
    iBase_EntityType entity_type)
{

  cout << "Volume smoothing is starting..." << endl;

  MsqError mError;
  const char* VERTEX_FIXED_TAG_NAME = "MesquiteVertexFixed";

  iBase_TagHandle fixed_tag = 0;

  iMesh::Error m_err = mk_core->imesh_instance()->getTagHandle(
      VERTEX_FIXED_TAG_NAME, fixed_tag);
  if (m_err) {
    m_err = mk_core->imesh_instance()->createTag(VERTEX_FIXED_TAG_NAME, 1,
        iBase_INTEGER, fixed_tag);
    IBERRCHK(m_err, "Trouble create a taghandle.");
  }

  MsqIMesh mesh_adapter(mk_core->imesh_instance()->instance(), volMesh,
      entity_type, mError, &fixed_tag);
  cout << "error =" << mError << endl;
  if (mError)
    throw mError;

  //get all the vertices in volume mesh
  int num_vtx, count;
  std::vector<iBase_EntityHandle> faces, verts;
  std::vector<int> indices, offsets;

  m_err = mk_core->imesh_instance()->getAdjEntIndices(volMesh, entity_type,
      iMesh_ALL_TOPOLOGIES, iBase_VERTEX, faces, verts, indices, offsets);
  IBERRCHK(m_err, "Trouble get the quads and nodes on the target surface.");
  num_vtx = verts.size();

  cout << "number of faces is " << faces.size() << endl;
  cout << "number of vertices is " << verts.size() << endl;

  //set fixed flag on all vertices
  vector<int> tag_data(num_vtx, 1);
  m_err = mk_core->imesh_instance()->setIntArrData(&verts[0], verts.size(),
      fixed_tag, &tag_data[0]);
  IBERRCHK(m_err, "Trouble set an array of int data for nodes on the target surface.");

  //clear fixed flag for vertices contained directly in set
  m_err = mk_core->imesh_instance()->getNumOfType(volMesh, iBase_VERTEX, count);
  IBERRCHK(m_err, "Trouble get the number of interior nodes on the target surface.");

  //get the interior mesh nodes on the target surface
  verts.clear();
  cout << "Num of Vertices on the target surface is " << count << endl;

  m_err = mk_core->imesh_instance()->getEntities(volMesh, iBase_VERTEX,
      iMesh_ALL_TOPOLOGIES, verts);
  IBERRCHK(m_err, "Trouble get the number of interior nodes on the target surface.");

  tag_data.clear();
  tag_data.resize(verts.size(), 0);
  m_err = mk_core->imesh_instance()->setIntArrData(&verts[0], verts.size(),
      fixed_tag, &tag_data[0]);
  IBERRCHK(m_err, "Trouble set the int data for interior nodes on the target surface.");

  //using the UntangleWrapper class to smooth the mesh with the inverted elements.
  UntangleWrapper smoother1;

  smoother1.set_cpu_time_limit(1000);
  smoother1.set_vertex_movement_limit_factor(0.001);

  smoother1.run_instructions(&mesh_adapter, mError);

  //Finally smooth the mesh
  //Use the ShapeImprover class to smooth the target surface mesh
  ShapeImprover smoother2;
  smoother2.set_cpu_time_limit(1000);
  smoother2.set_vertex_movement_limit_factor(0.001);
  smoother2.run_instructions(&mesh_adapter, mError);
  if (mError)
    cout << "Mesquite error in volume mesh smoothing is as follows\n" << mError
        << std::endl;
  /*

   SweepWrapper smoother3( 1e-6,  "COORDINATES_MAP");
   ShapeImprovementWrapper smoother2(mError);


   smoother2.run_instructions(&mesh_adapter, mError);

   if (mError)
   cout << "Mesquite error in volume mesh smoothing=" << mError << endl;
   */

}

MeshImprove::~MeshImprove()
{
  cout << "It is over now in smoothing" << endl;
}

}
