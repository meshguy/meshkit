

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

public:  
  virtual void setup_this();
  virtual void execute_this();
  virtual void facet(ModelEnt *curve);
  virtual void set_senses( ModelEnt *ent);
  virtual double length(  iGeom::EntityHandle vtx1, iMesh::EntityHandle vtx2);
  void set_facet_params(double faceting_tolerance = 0, double geom_resabs = 0);

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
