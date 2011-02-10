#include <iostream>
#include "CAMALSurfEval.hpp"
#include "meshkit/ModelEnt.hpp"
#include "meshkit/MKCore.hpp"
#include "meshkit/SizingFunction.hpp"

namespace MeshKit 
{
    
CAMALSurfEval::CAMALSurfEval(ModelEnt *me)
        : modelEnt(me)
{
}

CAMALSurfEval::~CAMALSurfEval() {}

double CAMALSurfEval::area() 
{
  return (modelEnt->dimension() != 2 ? -1.0 : modelEnt->measure());
}

void CAMALSurfEval::bounding_box(double box_min[3], double box_max[3]) 
{
  iGeom::Error err = modelEnt->mk_core()->igeom_instance()->getEntBoundBox(modelEnt->geom_handle(),
                                                                           box_min[0], box_min[1], box_min[2], 
                                                                           box_max[0], box_max[1], box_max[2]);
  IBERRCHK(err, "Trouble getting entity bounding box.");
}

void CAMALSurfEval::move_to_surface(double& x, double& y, double& z) 
{
  double close[3];
  modelEnt->evaluate(x, y, z, close);
  x = close[0];
  y = close[1];
  z = close[2];
}

void CAMALSurfEval::move_to_surface(double& x, double& y, double& z,
		     double& u_guess, double& v_guess)
{
  move_to_surface(x, y, z);
}

bool CAMALSurfEval::normal_at(double x, double y, double z, 
                              double& nx, double& ny, double& nz) 
{
  double norm[3];
  modelEnt->evaluate(x, y, z, NULL, norm);
  nx = norm[0];
  ny = norm[1];
  nz = norm[2];
  return true;
}

bool CAMALSurfEval::normal_at(double x, double y, double z, 
			      double& u_guess, double& v_guess,
			      double& nx, double& ny, double& nz)
{
  double norm[3];
  modelEnt->evaluate(x, y, z, NULL, norm);
  nx = norm[0];
  ny = norm[1];
  nz = norm[2];
  return true;
}

bool CAMALSurfEval::is_planar()
{
  std::string surf_type;
  iGeom::Error err = modelEnt->mk_core()->igeom_instance()->getFaceType(modelEnt->geom_handle(), surf_type);
  IBERRCHK(err, "Trouble getting surface type");
  return (surf_type == "plane");
}

bool CAMALSurfEval::is_parametric()
{
  bool is_param = false;
  iGeom::Error err = modelEnt->mk_core()->igeom_instance()->isEntParametric(modelEnt->geom_handle(), is_param);
  IBERRCHK(err, "Trouble getting whether entity is parametric");
  return is_param;
}

bool CAMALSurfEval::is_periodic_in_u(double& u_period)
{
  bool per_u = false, per_v = false;
  iGeom::Error err = modelEnt->mk_core()->igeom_instance()->isEntPeriodic(modelEnt->geom_handle(), per_u, per_v);
  IBERRCHK(err, "Trouble getting whether entity is periodic");
  return per_u;
}

bool CAMALSurfEval::is_periodic_in_v(double& v_period)
{
  bool per_u = false, per_v = false;
  iGeom::Error err = modelEnt->mk_core()->igeom_instance()->isEntPeriodic(modelEnt->geom_handle(), per_u, per_v);
  IBERRCHK(err, "Trouble getting whether entity is periodic");
  return per_v;
}

void CAMALSurfEval::get_param_range_u(double& u_low, double& u_high)
{
  double vmin, vmax;
  iGeom::Error err = modelEnt->mk_core()->igeom_instance()->getEntUVRange(modelEnt->geom_handle(), u_low, vmin, u_high, vmax);
  IBERRCHK(err, "Trouble getting entity parameter ranges.");
}

void CAMALSurfEval::get_param_range_v(double& v_low, double& v_high)
{
  double umin, umax;
  iGeom::Error err = modelEnt->mk_core()->igeom_instance()->getEntUVRange(modelEnt->geom_handle(), umin, v_low, umax, v_high);
  IBERRCHK(err, "Trouble getting entity parameter ranges.");
}

bool CAMALSurfEval::uv_from_position(double x, double y, double z, 
				     double& u, double& v)
  
{
  iGeom::Error err = modelEnt->mk_core()->igeom_instance()->getEntXYZtoUV(modelEnt->geom_handle(), x, y, z, u, v);
  IBERRCHK(err, "Trouble getting parameters from position.");
}

bool CAMALSurfEval::uv_from_position(double x, double y, double z, 
				     double& u, double& v,
				     double& cx, double& cy, double& cz)
{
  iGeom::Error err = modelEnt->mk_core()->igeom_instance()->getEntXYZtoUV(modelEnt->geom_handle(), x, y, z, u, v);
  IBERRCHK(err, "Trouble getting parameters from position.");
  err = modelEnt->mk_core()->igeom_instance()->getEntUVtoXYZ(modelEnt->geom_handle(), u, v, cx, cy, cz);
  IBERRCHK(err, "Trouble getting position from parameters.");
}

void CAMALSurfEval::position_from_uv(double u, double v, 
				     double& x, double& y, double& z)
{
  iGeom::Error err = modelEnt->mk_core()->igeom_instance()->getEntUVtoXYZ(modelEnt->geom_handle(), u, v, x, y, z);
  IBERRCHK(err, "Trouble getting position from parameters.");
}

void CAMALSurfEval::distortion_at_uv(double u, double v, 
				     double du[3], double dv[3])
{
  iGeom::Error err = modelEnt->mk_core()->igeom_instance()->getEnt1stDrvt(modelEnt->geom_handle(), u, v, 
                                                                          du[0], du[1], du[2], 
                                                                          dv[0], dv[1], dv[2]);
  IBERRCHK(err, "Trouble getting 1st derivative from parameters.");
}

} // namespace MeshKit
