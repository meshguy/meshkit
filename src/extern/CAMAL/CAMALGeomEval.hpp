#ifndef CAMAL_GEOM_EVAL_HPP
#define CAMAL_GEOM_EVAL_HPP

namespace MeshKit 
{

class ModelEnt;
    
#if CAMAL_VERSION < 500
#include "CMLGeomEval.hpp"

//! Implement CAMAL geometry callbacks using iGeom
class CAMALGeomEval : public CMLGeomEval
{
public:
  CAMALGeomEval(ModelEnt *ment);
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

  ModelEnt *modelEnt;
};
#else
#include "CMLSurfEval.hpp"

//! Implement CAMAL geometry callbacks using iGeom
class CAMALGeomEval : public CMLSurfEval
{
public:
  CAMALGeomEval(ModelEnt *ment);

  virtual ~CAMALGeomEval();
  
  virtual double area();

  virtual void bounding_box(double box_min[3], double box_max[3]);

  virtual void move_to_surface(double& x, double& y, double& z);

  virtual void move_to_surface(double& x, double& y, double& z,
			       double& u_guess, double& v_guess);

  virtual bool normal_at(double x, double y, double z, 
                         double& nx, double& ny, double& nz);

  virtual bool normal_at(double x, double y, double z, 
                         double& u_guess, double& v_guess,
                         double& nx, double& ny, double& nz);
  
  virtual bool is_planar();

  virtual bool is_parametric();

  virtual bool is_periodic_in_u(double& u_period);

  virtual bool is_periodic_in_v(double& v_period);

  virtual void get_param_range_u(double& u_low, double& u_high);

  virtual void get_param_range_v(double& v_low, double& v_high);
  
  virtual bool uv_from_position(double x, double y, double z, 
                                double& u, double& v);

  virtual bool uv_from_position(double x, double y, double z, 
                                double& u, double& v,
                                double& cx, double& cy, double& cz);

  virtual void position_from_uv(double u, double v, 
				double& x, double& y, double& z);

  virtual void distortion_at_uv(double u, double v, 
                                double du[3], double dv[3]);

  void set_mesh_size(double tmp_size);
  double get_mesh_size();
  
  int get_dimension() const { return myDimension; }

  // a new method to compute the intersection with a ray
  bool pierce_surface_with_ray(double & x, double & y, double & z, double dir_x,
         double dir_y, double dir_z);

private:

  ModelEnt *modelEnt;
};
#endif

} // namespace MeshKit

#endif
