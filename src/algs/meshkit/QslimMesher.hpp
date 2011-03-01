//-----------------------------------C++-------------------------------------//
// File: src/algs/meshkit/QslimMesher.hpp
//
// Brief: QslimMesher class definition:
//        Operates on a set of triangles, options passed with another
//         class, QslimOptions
//---------------------------------------------------------------------------//

#ifndef MESHKIT_QSLIM_MESHER_HPP
#define MESHKIT_QSLIM_MESHER_HPP

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string>
#include <iostream>
#include <fstream>
#include <string.h>
#include <limits.h>

//#include "meshkit/iGeom.hpp"
#include <set>
#include <vector>

#include "meshkit/MeshScheme.hpp"
#include "meshkit/QslimOptions.hpp"

namespace MeshKit {
//===========================================================================//
/*!
 * \class QslimMesher
 * \brief Decimate a set of triangles
 *
 * QslimMesher decimates a set of triangles that form a 3d surface
 * It uses edge collapse sequentially, while keeping the error in
 * the quadric sense minimal at each step.
 */
//===========================================================================//

using namespace std;

class QslimDecimation;

class QslimMesher: public MeshScheme {

public:
  //construction function for Qslim mesher
  QslimMesher(MKCore *mk_core, const MEntVector &me_vec);

  //set up the parameters for decimation meshing
  // these will be passed with QslimOptions
  virtual void setup_this();

  //Decimate the mesh
  virtual void execute_this();

  /**\brief Get class name */
   static const char* name()
     { return "QslimMesher"; }

   // pass
  void set_options(QslimOptions & opts)
  {
    _opts = opts;
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

  ~QslimMesher();

private:

  QslimOptions _opts;
  QslimDecimation * _worker;

};

}

#endif
