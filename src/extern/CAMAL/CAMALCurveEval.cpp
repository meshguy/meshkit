#include "CAMALCurveEval.hpp"
#include "meshkit/ModelEnt.hpp"
#include "meshkit/MKCore.hpp"
#include "meshkit/iGeom.hh"
#include "meshkit/Error.hpp"

namespace MeshKit 
{

double CAMALCurveEval::arc_length()
{
  return modelEnt->measure();
}

bool CAMALCurveEval::is_parametric()
{
}

bool CAMALCurveEval::is_periodic(double& period)
{
  bool is_peru, is_perv;
  iGeom::Error err = modelEnt->mk_core()->igeom_instance()->isEntPeriodic(modelEnt->geom_handle(), is_peru, is_perv);
  IBERRCHK(err, "Trouble calling isEntPeriodic.");
  double umin, umax;
  err = modelEnt->mk_core()->igeom_instance()->getEntURange(modelEnt->geom_handle(), umin, umax);
  IBERRCHK(err, "Trouble calling getEntURange.");
  period = umax - umin;
  return is_peru;
}

void CAMALCurveEval::get_param_range(double& u_start, double& u_end)
{
  iGeom::Error err = modelEnt->mk_core()->igeom_instance()->getEntURange(modelEnt->geom_handle(), u_start, u_end);
  IBERRCHK(err, "Trouble calling getEntURange.");
}

double CAMALCurveEval::u_from_arc_length(double u_root, double arc_length)
{
  double u_start, u_end;
  iGeom::Error err = modelEnt->mk_core()->igeom_instance()->getEntURange(modelEnt->geom_handle(), u_start, u_end);
  IBERRCHK(err, "Trouble calling getEntURange.");
  if (u_end-u_start <= 0) return arc_length/modelEnt->measure();
  else return u_root + (arc_length/modelEnt->measure())*(u_end-u_start);
}
             
bool CAMALCurveEval::position_from_u(double u, 
                             double& x, double& y, double& z )
{
  iGeom::Error err = modelEnt->mk_core()->igeom_instance()->getEntUtoXYZ(modelEnt->geom_handle(), u, x, y, z);
  IBERRCHK(err, "Trouble calling getEntUtoXYZ.");
  return true;
}

void CAMALCurveEval::move_to_curve(double& x, double& y, double& z) 
{
  double dum[3];
  modelEnt->evaluate(x, y, z, dum);
  x = dum[0];
  y = dum[1];
  z = dum[2];
}

double CAMALCurveEval::u_from_position(double x, double y, double z)
{
  double u;
  iGeom::Error err = modelEnt->mk_core()->igeom_instance()->getEntXYZtoU(modelEnt->geom_handle(), x, y, z, u);
  IBERRCHK(err, "Trouble calling getEntXYZtoU.");
  return u;
}

void CAMALCurveEval::start_coordinates(double& x, double& y, double& z)
{
  double u_start, u_end;
  iGeom::Error err = modelEnt->mk_core()->igeom_instance()->getEntURange(modelEnt->geom_handle(), u_start, u_end);
  IBERRCHK(err, "Trouble calling getEntURange.");
  position_from_u(u_start, x, y, z);
}

void CAMALCurveEval::end_coordinates(double& x, double& y, double& z)
{
  double u_start, u_end;
  iGeom::Error err = modelEnt->mk_core()->igeom_instance()->getEntURange(modelEnt->geom_handle(), u_start, u_end);
  IBERRCHK(err, "Trouble calling getEntURange.");
  position_from_u(u_end, x, y, z);
}

} // namespace MeshKit

