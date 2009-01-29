#ifndef CAMAL_GEOM_EVAL_HPP
#define CAMAL_GEOM_EVAL_HPP

#include "CMLGeomEval.hpp"
#include "iGeom.h"

//! Implement CAMAL geometry callbacks using iGeom
class CAMALGeomEval : public CMLGeomEval
{
public:
  CAMALGeomEval(iGeom_Instance geomIface, iBase_EntityHandle gent);
  virtual ~CAMALGeomEval();
  
  virtual double area();

  virtual void bounding_box(double box_min[3], double box_max[3]);

  virtual void move_to(double& x, double& y, double& z);

  virtual bool normal_at(const double& x, const double& y, const double& z,
                         double& nx, double& ny, double& nz);

  void set_mesh_size(double tmp_size);
  virtual double get_mesh_size();
  
  int get_dimension() const { return myDimension; }

private:

  iGeom_Instance geomIface;
  iBase_EntityHandle myEnt;
  int myDimension;
  double meshSize;
};

#endif
