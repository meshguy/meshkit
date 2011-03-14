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
  iGeom::Error err = modelEnt->igeom_instance()->getEntBoundBox(modelEnt->geom_handle(),
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
  iGeom::Error err = modelEnt->igeom_instance()->getFaceType(modelEnt->geom_handle(), surf_type);
  if (iBase_SUCCESS != err || (surf_type != "plane")) return false;
  else return true;
}

bool CAMALSurfEval::is_parametric()
{
  bool is_param = false;
  iGeom::Error err = modelEnt->igeom_instance()->isEntParametric(modelEnt->geom_handle(), is_param);
  IBERRCHK(err, "Trouble getting whether entity is parametric");
  return is_param;
}

bool CAMALSurfEval::is_periodic_in_u(double& u_period)
{
  bool per_u = false, per_v = false;
  iGeom::Error err = modelEnt->igeom_instance()->isEntPeriodic(modelEnt->geom_handle(), per_u, per_v);
  IBERRCHK(err, "Trouble getting whether entity is periodic");
  return per_u;
}

bool CAMALSurfEval::is_periodic_in_v(double& v_period)
{
  bool per_u = false, per_v = false;
  iGeom::Error err = modelEnt->igeom_instance()->isEntPeriodic(modelEnt->geom_handle(), per_u, per_v);
  IBERRCHK(err, "Trouble getting whether entity is periodic");
  return per_v;
}

void CAMALSurfEval::get_param_range_u(double& u_low, double& u_high)
{
  double vmin, vmax;
  iGeom::Error err = modelEnt->igeom_instance()->getEntUVRange(modelEnt->geom_handle(), u_low, vmin, u_high, vmax);
  IBERRCHK(err, "Trouble getting entity parameter ranges.");
}

void CAMALSurfEval::get_param_range_v(double& v_low, double& v_high)
{
  double umin, umax;
  iGeom::Error err = modelEnt->igeom_instance()->getEntUVRange(modelEnt->geom_handle(), umin, v_low, umax, v_high);
  IBERRCHK(err, "Trouble getting entity parameter ranges.");
}

bool CAMALSurfEval::uv_from_position(double x, double y, double z, 
				     double& u, double& v)
  
{
  iGeom::Error err = modelEnt->igeom_instance()->getEntXYZtoUV(modelEnt->geom_handle(), x, y, z, u, v);
  IBERRCHK(err, "Trouble getting parameters from position.");
}

bool CAMALSurfEval::uv_from_position(double x, double y, double z, 
				     double& u, double& v,
				     double& cx, double& cy, double& cz)
{
  iGeom::Error err = modelEnt->igeom_instance()->getEntXYZtoUV(modelEnt->geom_handle(), x, y, z, u, v);
  IBERRCHK(err, "Trouble getting parameters from position.");
  err = modelEnt->igeom_instance()->getEntUVtoXYZ(modelEnt->geom_handle(), u, v, cx, cy, cz);
  IBERRCHK(err, "Trouble getting position from parameters.");
}

void CAMALSurfEval::position_from_uv(double u, double v, 
				     double& x, double& y, double& z)
{
  iGeom::Error err = modelEnt->igeom_instance()->getEntUVtoXYZ(modelEnt->geom_handle(), u, v, x, y, z);
  IBERRCHK(err, "Trouble getting position from parameters.");
}

void CAMALSurfEval::distortion_at_uv(double u, double v, 
				     double du[3], double dv[3])
{
  iGeom::Error err = modelEnt->igeom_instance()->getEnt1stDrvt(modelEnt->geom_handle(), u, v,
                                                                          du[0], du[1], du[2], 
                                                                          dv[0], dv[1], dv[2]);
  IBERRCHK(err, "Trouble getting 1st derivative from parameters.");
}

// this could be inlined
double tArea(double * a, double *b, double *c, double * normal)
{
  double result = 0;
  // ( ( B-A ) X (C-A) ) * normal
  double AB[3] = { b[0] - a[0], b[1] - a[1], b[2] - a[2] };
  double AC[3] = { c[0] - a[0], c[1] - a[1], c[2] - a[2] };
  result += (AB[1] * AC[2] - AB[2] * AC[1]) * normal[0];
  result += (AB[2] * AC[0] - AB[0] * AC[2]) * normal[1];
  result += (AB[0] * AC[1] - AB[1] * AC[0]) * normal[2];
}

// decide if the boundary loops are positively oriented around normal
// use the normal at the first node on the boundary
// It should work even if the first node is on an interior loop, not external loop
// it should work in most cases, except the surface is highly distorted, and the projected area of an internal
// loop is higher than the projected area of an external loop
// it is specific to Camal advancing front algorithms (Paver and TriAdvance)
void CAMALSurfEval::correct_orientation(std::vector<int> & loop_sizes,
    std::vector<int> & loops, std::vector<double> & bdy_coords)
{
  // first, normal at the initial point on the boundary
  double normal[3] = { 0, 0, 0 };
  /*this->*/normal_at(bdy_coords[0], bdy_coords[1], bdy_coords[2], normal[0],
      normal[1], normal[2]);

  double oriented_area = 0.;
  unsigned int start_current_loop = 0;
  for (unsigned int k = 0; k < loop_sizes.size(); k++)
  {
    // for each loop, compute the oriented area of each triangle
    int current_loop_size = loop_sizes[k];
    unsigned int startIndex = loops[start_current_loop];
    for (unsigned int i = 1; i < current_loop_size - 1; i++)
    {
      unsigned int i1 = loops[start_current_loop + i];
      unsigned int i2 = loops[start_current_loop + (i + 1)];

      double oArea = tArea(&bdy_coords[3 * startIndex], &bdy_coords[3 * i1],
          &bdy_coords[3 * i2], normal);
      oriented_area += oArea;
    }
    start_current_loop += current_loop_size;
  }
  if (oriented_area < 0.)
  { // correct orientation
    unsigned int start_current_loop = 0;
    for (unsigned int k = 0; k < loop_sizes.size(); k++)
    {
      // for each loop, switch index 1 with last, 2, with first before last...
      int current_loop_size = loop_sizes[k];
      std::reverse(&loops[start_current_loop + 1], &loops[start_current_loop
          + current_loop_size]);
      start_current_loop += current_loop_size;
    }
  }
}

} // namespace MeshKit
