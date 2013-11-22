/*
 * MBSplitOp.hpp
 *
 *  Created on: Oct 2, 2011
 *      Author: iulian
 */

#ifndef MBSPLITOP_HPP_
#define MBSPLITOP_HPP_

#include "meshkit/MeshScheme.hpp"

namespace MeshKit {

class MBSplitOp: public MeshKit::MeshScheme
{
public:
  MBSplitOp(MKCore *mk_core, const MEntVector &me_vec);
  virtual ~MBSplitOp();

  // set options polyline[3*nPoints], closed (1 or 0, closed loop or open polyline,
  //                                           intersect  the boundary)
  //  globalId: to identify / single out the gface to be split
  // the gface will be found among gfaces in the mentSelection
  void set_options(int globalId,
      double dirx, double diry, double dirz, int closed, double min_dot);

  void add_points(double x, double y, double z);
  //set up the splitting of a gface, with a polyline or loop
  virtual void setup_this();

  // construct the mesh
  virtual void execute_this();

  /**\brief Get class name */
  static const char* name()
  {
    return "MBSplitOp";
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
  int _globalId;
  std::vector<double>_polyline;
  double _direction[3];
  int _closed;
  double _min_dot;
};

}

#endif /* MBSPLITOP_HPP_ */
