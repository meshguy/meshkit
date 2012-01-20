/*
 * MBVolOp.h
 *
 *  Created on: Jan 13, 2012
 */

#ifndef MBVOLOP_H_
#define MBVOLOP_H_

#include "meshkit/MeshScheme.hpp"
#include "moab/GeomTopoTool.hpp"
#include "moab/FBEngine.hpp"

namespace MeshKit {

class MBVolOp: public MeshKit::MeshScheme {
public:
  MBVolOp(MKCore *mk_core, const MEntVector &me_vec);
  virtual ~MBVolOp();
  //  globalIds: to identify / single out the gfaces to be weaved
    // the gface will be found among gfaces in the mentSelection
  // gface sources, targets
  void add_pair(int source, int target);// these are global ids in input ment vectors
  // or maybe orders in some lists?

  // this is the main direction of sweeping / weaving
  void set_direction(double x, double y, double z)
      {_direction[0]=x; _direction[1]=y; _direction[2]=z;}
  //set up the weaving, needed for volume creation
  virtual void setup_this();

  // construct the mesh
  virtual void execute_this();

  /**\brief Get class name */
  static const char* name()
  {
    return "MBVolOp";
  }

  static bool can_mesh(iBase_EntityType dim)
  {
    return false;
  }

  /** \brief Function returning whether this scheme can mesh the specified entity
   *
   * Used by MeshOpFactory to find scheme for an entity.
   * \param me ModelEnt being queried
   * \return If true, this scheme can mesh the specified ModelEnt
   */
  static bool can_mesh(ModelEnt *me)
  {
    return false;// this is not used for meshing, but for geometry creation
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


private:
  /**
     * \brief find correspondence between vertices, edges, surfaces, for later weaving
     * used later to weave between source and target
     */
  void establish_mapping();

  double _direction[3];
  moab::GeomTopoTool * _pGTT ;
  moab::EntityHandle _rootSet;// this is the root set for gtt result
  std::map<moab::EntityHandle, moab::EntityHandle> vertexMap;
  std::map<moab::EntityHandle, moab::EntityHandle> edgeMap;
  std::map<moab::EntityHandle, moab::EntityHandle> faceMap;
  moab::FBEngine * _fbe;
};

}// namespace MeshKit

#endif /* MBVOLOP_H_ */
