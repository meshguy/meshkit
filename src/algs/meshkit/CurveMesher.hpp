

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string>
#include <iostream>
#include <fstream>
#include <string.h>
#include <limits.h>

#include <iGeom.h>
#include <iMesh.h>
#include <set>
#include <iRel.h>
#include <vector>
#include "meshkit/MKCore.hpp"
#include "meshkit/MeshScheme.hpp"

namespace MeshKit
{
 
using namespace std;

class CurveMesher : public MeshScheme
{
public: 
CurveMesher(MKCore *mk, const MEntVector &ments);

~CurveMesher();
private:
  double facet_tol;
  double geom_res;
  MKCore *mk;
  MEntVector model_ents;

  /** \brief Returns the distance between an iGeom vertex and an iMesh vertex.
   * \param vtx1 iGeom vertex handle
   * \param vtx2 iMesh vertex handle
   */
  virtual double vtx2vtx_dist(iGeom::EntityHandle vtx1, iMesh::EntityHandle vtx2);

  /** \brief Returns the distance between an iGeom vertex and an iMesh vertex.
   * \param vtx1 iMesh vertex handle
   * \param vtx2 iMesh vertex handle
   */

  virtual void facet(ModelEnt *curve);

  /** \brief Sets the senses wrt all surfaces adjacent to the curve
   * \param curve Pointer to the ModelEnt to be meshed
   */
  virtual void set_senses( ModelEnt *curve);

public:  
  virtual void setup_this();
  virtual void execute_this();

  /** \brief Sets the faceting tolerance and geom_reabs values. If
       this function is not run before mk->setup(). Default values
       for these parameters will be used. 
   * \param faceting_tolerance value to be set for the faceting tolerance
   * \param geom_resabs value to be set for the geom_resabs (used for vertex proximity checks)
   */
  void set_mesh_params(double faceting_tolerance = 0, double geom_resabs = 0);

  static bool can_mesh(iBase_EntityType dim)
  { return iBase_EDGE == dim; }

  static bool can_mesh(ModelEnt *me)
  { return canmesh_edge(me); }
 
  static const char* name()
  {return "CurveMesher";}

  static const moab::EntityType* output_types();

  virtual const moab::EntityType* mesh_types_arr() const
  { return output_types(); }

};


}
