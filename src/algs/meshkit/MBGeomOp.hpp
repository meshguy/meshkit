/*
 * MBGeomOp.hpp
 *
 *  Created on: Sep 30, 2011
 *      Author: iulian
 */

#ifndef MBGEOMOP_HPP_
#define MBGEOMOP_HPP_

#include "meshkit/MeshScheme.hpp"

namespace MeshKit {

class MBGeomOp: public virtual MeshKit::MeshScheme
{
public:
  MBGeomOp(MKCore *mk_core, const MEntVector &me_vec);
  virtual ~MBGeomOp();

  //set up the geometrization of a model ent (face) for mesh-based geometry
  virtual void setup_this();

  // construct the mesh
  virtual void execute_this();

  /**\brief Get class name */
  static const char* name()
  {
    return "MBGeomOp";
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

};

}

#endif /* MBGEOMOP_HPP_ */
