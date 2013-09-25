/*
 * TriangleMesher.h
 *
 *  Created on: Sep 23, 2011
 *      Author: iulian
 */

#ifndef TRIANGLEMESHER_H_
#define TRIANGLEMESHER_H_


#include <string>
#include "meshkit/MeshScheme.hpp"

namespace MeshKit {

class TriangleMesher: public MeshScheme
{
public:
  virtual ~TriangleMesher();

  //construction function for Triangle mesher
  TriangleMesher(MKCore *mk_core, const MEntVector &me_vec);

  //set up the parameters for triangle meshing
  virtual void setup_this();

  // construct the mesh
  virtual void execute_this();

  /**\brief Get class name */
  static const char* name()
  {
    return "TriangleMesher";
  }

  // pass
  void set_options(char * opts, int direction, double fretting)
  {
    _opts = opts;
    _dir = direction;
    _fretting = fretting;
  }

  static bool can_mesh(iBase_EntityType dim)
  {
    return iBase_FACE == dim;
  }

  /** \brief Function returning whether this scheme can mesh the specified entity
   *
   * Used by MeshOpFactory to find scheme for an entity.
   * \param me ModelEnt being queried
   * \return If true, this scheme can mesh the specified ModelEnt
   */
  static bool can_mesh(ModelEnt *me)
  {
    return canmesh_face(me);
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

  char * _opts;
  int _dir; // do it in x, y, z normal planes
  double _fretting;
};

}

#endif /* TRIANGLEMESHER_H_ */
