//-----------------------------------C++-------------------------------------//
// File: src/algs/OneToOneSwept.hpp
// Wednesday February 11 10:50 2011
// Brief: OneToOneSwept class definition: generate the all-hexahedral mesh by
//        general sweeping 
//---------------------------------------------------------------------------//


#ifndef MESHKIT_ONETOONESWEPT_HPP
#define MESHKIT_ONETOONESWEPT_HPP

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string>
#include <iostream>
#include <fstream>
#include <string.h>
#include <limits.h>

#include "Global.hpp"
#include "meshkit/MKCore.hpp"
#include "meshkit/ModelEnt.hpp"
#include "meshkit/SizingFunction.hpp"

#include <iMesh.h>
#include <iGeom.h>
#include <set>
#include <iRel.h>
#include <vector>
#include <set>
#include <list>

#include "meshkit/MeshScheme.hpp"

using namespace std;

namespace MeshKit {
//===========================================================================//
/*!
 * \class OneToOneSwept
 * \brief Generate the all-hexahedral mesh by sweeping
 *
 * OneToOneSwept generates the all-hexahedral mesh by sweeping the source mesh
 * to the target surface
 */
//===========================================================================//

class OneToOneSwept: public MeshScheme {
public:

  OneToOneSwept(MKCore *mk_core, const MEntVector &me_vec);
  ~OneToOneSwept();

  virtual void setup_this();
  virtual void execute_this();

  /**\brief Get class name */
  static const char* name()
  {
    return "OneToOneSwept";
  }

  /**\brief Function returning whether this scheme can mesh entities of t
   *        the specified dimension.
   *\param dim entity dimension
   */
  static bool can_mesh(iBase_EntityType dim)
  {
    return iBase_REGION == dim;
  }

  /** \brief Function returning whether this scheme can mesh the specified entity
   *
   * Used by MeshOpFactory to find scheme for an entity.
   * \param me ModelEnt being queried
   * \return If true, this scheme can mesh the specified ModelEnt
   */
  static bool can_mesh(ModelEnt *me)
  {
    return canmesh_region(me);
  }

  /**\brief Get list of mesh entity types that can be generated.
   *\return array terminated with \c moab::MBMAXTYPE
   */
  static const moab::EntityType* output_types();

  /** \brief Return the mesh entity types operated on by this scheme
   * \return array terminated with \c moab::MBMAXTYPE
   */
  virtual const moab::EntityType* mesh_types_arr() const
  {
    return output_types();
  }

  //specify the source surface for OneToOneSwept class
  void SetSourceSurface(int index);

  //specify the target surface for OneToOneSwept class
  void SetTargetSurface(int index);

private:
#ifdef HAVE_ARMADILLO
  void SurfMeshHarmonic(iBase_EntityHandle vol);
#endif

  void BuildLateralBoundaryLayers(ModelEnt * me, std::vector<moab::EntityHandle> & layers);

  moab::ErrorCode NodeAbove(moab::EntityHandle node1, moab::EntityHandle node2, moab::Range & quadsOnLateralSurfaces,
      moab::EntityHandle & node3, moab::EntityHandle & node4);

  moab::ErrorCode FourthNodeInQuad(moab::EntityHandle node1, moab::EntityHandle node2, moab::EntityHandle node3,
      moab::Range & quadsOnLateralSurfaces, moab::EntityHandle & node4);

  //function for all-quad meshing on the target surface
  int TargetSurfProjection(std::vector<moab::EntityHandle> & boundLayers);

  // to compute node positions in the interior of volume, numLayers-1 times
  // use similar code to TargetSurfProjection, but do not project on surface...
  int ProjectInteriorLayers(std::vector<moab::EntityHandle> & boundLayers, vector<vector<Vertex> > &linkVertexList);

  //create the hexahedral elements between the source surface and target surface
  int CreateElements(vector<vector<Vertex> > &linkVertexList);

  // input: list of nodes on source, boundary center, list of nodes on target, target center
  // output: 3x3 matrix A such that
  //  target= A * ( source - 2*sc + tc) + sc
  void computeTransformationFromSourceToTarget(std::vector<Vector3D> & sNodes, Vector3D & sc,
      std::vector<Vector3D> & tNodes, Vector3D & tc, Matrix3D & transMatrix);

#ifdef HAVE_MESQUITE
  //target surface mesh by Mesquite
  void SurfMeshOptimization();
#endif

private:
  //private member variable
  iBase_EntityHandle sourceSurface;
  iBase_EntityHandle targetSurface;
  iBase_EntityHandle volume;
  moab::Tag markTag;
  moab::Interface *mb;
  int numLayers;
  int numLoops;
  int sizeBLayer;
  iGeom * igeom_inst;
  int index_src, index_tar;
  iBase_TagHandle geom_id_tag;

  std::vector<Vertex> NodeList; //mesh nodes on the source surface
  std::vector<Vertex> TVertexList; //mesh nodes on the target surface

  std::vector<Face> FaceList;
  iBase_EntitySetHandle volumeSet;

};

}

#endif

