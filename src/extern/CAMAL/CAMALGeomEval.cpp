#include <iostream>
#include "CAMALGeomEval.hpp"
#include "meshkit/ModelEnt.hpp"
#include "meshkit/MKCore.hpp"
#include "meshkit/SizingFunction.hpp"

namespace MeshKit 
{
    
#if CAMAL_VERSION < 500
void CAMALGeomEval::set_mesh_size(double tmp_size) 
{
  assert(false);
  MKCHKERR(MK_BAD_INPUT, "Shouldn't be setting size on CAMALGeomEval entity, use MeshKit methods.");
}

double CAMALGeomEval::get_mesh_size() 
{
    // check for sizing function from MKCore
  int ind = modelEnt->sizing_function_index();
  if (-1 != ind && modelEnt->mk_core()->sizing_function(ind)) 
    return modelEnt->mk_core()->sizing_function(ind)->size();
  
  MKCHKERR(MK_BAD_INPUT, "Sizing function not set yet for %s %d", modelEnt->me_type(), modelEnt->id());
  return -1.0;
}

CAMALGeomEval::CAMALGeomEval(ModelEnt *ment) 
        : modelEnt(ment)
{
}

CAMALGeomEval::~CAMALGeomEval() {}

double CAMALGeomEval::area() 
{
    // find the area of this entity
  if (modelEnt->dimension() != 2)
    return -1.0;
  
  return modelEnt->measure();
}

void CAMALGeomEval::bounding_box(double box_min[3], double box_max[3]) 
{
  iGeom::Error err = modelEnt->mk_core()->igeom_instance()->
      getEntBoundBox(modelEnt->geom_handle(), 
                     box_min, box_min+1, box_min+2,
                     box_max, box_max+1, box_max+2);
  IBCHKERR(err, "Trouble getting bounding box of %s %d", modelEnt->me_type(), modelEnt->id());
}

void CAMALGeomEval::move_to(double& x, double& y, double& z) 
{
  double tmp_close[3];
  modelEnt->evaluate(x, y, z, tmp_close);
  x = tmp_close[0];
  y = tmp_close[1];
  z = tmp_close[2];
}

bool CAMALGeomEval::normal_at(const double& x, const double& y, const double& z,
                              double& nx, double& ny, double& nz) 
{
  double tmp_norm[3];
  modelEnt->evaluate(x, y, z, NULL, tmp_norm);
  nx = tmp_norm[0];
  ny = tmp_norm[1];
  nz = tmp_norm[2];

  return true;
}

#else

bool debug_surf_eval = false;

void CAMALGeomEval::set_mesh_size(double tmp_size) 
{
  SizingFunction *sf = modelEnt->mk_core()->sizing_function(tmp_size);
  modelEnt->sizing_function_index(sf->core_index());
}

double CAMALGeomEval::get_mesh_size() 
{
  return modelEnt->mesh_interval_size();
}

CAMALGeomEval::CAMALGeomEval(ModelEnt *me)
        : modelEnt(me)
{
}

CAMALGeomEval::~CAMALGeomEval() {}

double CAMALGeomEval::area() 
{
  return (modelEnt->dimension() != 2 ? -1.0 : modelEnt->measure());
}

void CAMALGeomEval::bounding_box(double box_min[3], double box_max[3]) 
{
  iGeom::Error err = modelEnt->mk_core()->igeom_instance()->getEntBoundBox(modelEnt->geom_handle(),
                                                                           box_min[0], box_min[1], box_min[2], 
                                                                           box_max[0], box_max[1], box_max[2]);
  IBERRCHK(err, "Trouble getting entity bounding box.");
}

void CAMALGeomEval::move_to_surface(double& x, double& y, double& z) 
{
  double close[3];
  modelEnt->evaluate(x, y, z, close);
  x = close[0];
  y = close[1];
  z = close[2];
}

void CAMALGeomEval::move_to_surface(double& x, double& y, double& z,
		     double& u_guess, double& v_guess)
{
  move_to_surface(x, y, z);
}

bool CAMALGeomEval::normal_at(double x, double y, double z, 
                              double& nx, double& ny, double& nz) 
{
  double norm[3];
  modelEnt->evaluate(x, y, z, NULL, norm);
  nx = norm[0];
  ny = norm[1];
  nz = norm[2];
  return true;
}

bool CAMALGeomEval::normal_at(double x, double y, double z, 
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

bool CAMALGeomEval::is_planar()
{
  std::string surf_type;
  iGeom::Error err = modelEnt->mk_core()->igeom_instance()->getFaceType(modelEnt->geom_handle(), surf_type);
  IBERRCHK(err, "Trouble getting surface type");
  return (surf_type == "plane");
}

bool CAMALGeomEval::is_parametric()
{
  bool is_param = false;
  iGeom::Error err = modelEnt->mk_core()->igeom_instance()->isEntParametric(modelEnt->geom_handle(), is_param);
  IBERRCHK(err, "Trouble getting whether entity is parametric");
  return is_param;
}

bool CAMALGeomEval::is_periodic_in_u(double& u_period)
{
  bool per_u = false, per_v = false;
  iGeom::Error err = modelEnt->mk_core()->igeom_instance()->isEntPeriodic(modelEnt->geom_handle(), per_u, per_v);
  IBERRCHK(err, "Trouble getting whether entity is periodic");
  return per_u;
}

bool CAMALGeomEval::is_periodic_in_v(double& v_period)
{
  bool per_u = false, per_v = false;
  iGeom::Error err = modelEnt->mk_core()->igeom_instance()->isEntPeriodic(modelEnt->geom_handle(), per_u, per_v);
  IBERRCHK(err, "Trouble getting whether entity is periodic");
  return per_v;
}

void CAMALGeomEval::get_param_range_u(double& u_low, double& u_high)
{
  double vmin, vmax;
  iGeom::Error err = modelEnt->mk_core()->igeom_instance()->getEntUVRange(modelEnt->geom_handle(), u_low, vmin, u_high, vmax);
  IBERRCHK(err, "Trouble getting entity parameter ranges.");
}

void CAMALGeomEval::get_param_range_v(double& v_low, double& v_high)
{
  double umin, umax;
  iGeom::Error err = modelEnt->mk_core()->igeom_instance()->getEntUVRange(modelEnt->geom_handle(), umin, v_low, umax, v_high);
  IBERRCHK(err, "Trouble getting entity parameter ranges.");
}

bool CAMALGeomEval::uv_from_position(double x, double y, double z, 
				     double& u, double& v)
  
{
  iGeom::Error err = modelEnt->mk_core()->igeom_instance()->getEntXYZtoUV(modelEnt->geom_handle(), x, y, z, u, v);
  IBERRCHK(err, "Trouble getting parameters from position.");
}

bool CAMALGeomEval::uv_from_position(double x, double y, double z, 
				     double& u, double& v,
				     double& cx, double& cy, double& cz)
{
  iGeom::Error err = modelEnt->mk_core()->igeom_instance()->getEntXYZtoUV(modelEnt->geom_handle(), x, y, z, u, v);
  IBERRCHK(err, "Trouble getting parameters from position.");
  err = modelEnt->mk_core()->igeom_instance()->getEntUVtoXYZ(modelEnt->geom_handle(), u, v, cx, cy, cz);
  IBERRCHK(err, "Trouble getting position from parameters.");
}

void CAMALGeomEval::position_from_uv(double u, double v, 
				     double& x, double& y, double& z)
{
  iGeom::Error err = modelEnt->mk_core()->igeom_instance()->getEntUVtoXYZ(modelEnt->geom_handle(), u, v, x, y, z);
  IBERRCHK(err, "Trouble getting position from parameters.");
}

void CAMALGeomEval::distortion_at_uv(double u, double v, 
				     double du[3], double dv[3])
{
  iGeom::Error err = modelEnt->mk_core()->igeom_instance()->getEnt1stDrvt(modelEnt->geom_handle(), u, v, 
                                                                          du[0], du[1], du[2], 
                                                                          dv[0], dv[1], dv[2]);
  IBERRCHK(err, "Trouble getting 1st derivative from parameters.");
}

bool CAMALGeomEval::pierce_surface_with_ray(double & x, double & y, double & z, double dir_x,
         double dir_y, double dir_z)
{
  std::vector<iGeom::EntityHandle> entities;
  std::vector<double> points, params;
  iGeom::Error err = modelEnt->mk_core()->igeom_instance()->getPntRayIntsct(x, y, z, dir_x, dir_y, dir_z,
                                                                            iBase_INTERLEAVED,
                                                                            entities, points, params);
  IBERRCHK(err, "Trouble getting ray intersection.");
  if (entities.empty() || points.size() < 3) return false;
  
    // intersection found, populate with 1st
  x = points[0];
  y = points[1];
  z = points[2];

  return true;
}

int CAMALGeomEval::get_dimension() const { return modelEnt->dimension(); }

#endif
} // namespace MeshKit
