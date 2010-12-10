/**
 * Copyright 2006 Sandia Corporation.  Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Coroporation, the U.S. Government
 * retains certain rights in this software.
 * 
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * 
 */
#include <iostream>
#include "CAMALGeomEval.hpp"
#include "stdlib.h"  // include it just for free()

#if CAMAL_VERSION < 500
void CAMALGeomEval::set_mesh_size(double tmp_size) 
{
  meshSize = tmp_size;
}

double CAMALGeomEval::get_mesh_size() 
{
  return meshSize;
}

CAMALGeomEval::CAMALGeomEval(iGeom_Instance geom_iface, iBase_EntityHandle gent) 
  : geomIface( geom_iface ),
    myEnt( gent ),
    myDimension( -1 )
{
  int gtype, result;
  iGeom_getEntType(geomIface, gent, &gtype, &result);
  if (iBase_SUCCESS != result) {
    std::cerr << "Trouble getting dimension of gentity." << std::endl;
  }
  else 
    myDimension = gtype;
}

CAMALGeomEval::~CAMALGeomEval() {}

double CAMALGeomEval::area() 
{
    // find the area of this entity
  if (myDimension != 2) 
    return -1.0;
  
  double measure = -1.0;
  double *measure_ptr = &measure;
  int measure_size, measure_alloc = 1;
  int result;
  iGeom_measure(geomIface, &myEnt, 1, &measure_ptr, &measure_alloc, &measure_size, &result);
  if (iBase_SUCCESS != result) {
    std::cerr << "Trouble getting area of gentity." << std::endl;
    return -1.0;
  }
  
  return measure;
}

void CAMALGeomEval::bounding_box(double box_min[3], double box_max[3]) 
{
  int result;
  iGeom_getEntBoundBox(geomIface, myEnt,  
                       box_min, box_min+1, box_min+2,
                       box_max, box_max+1, box_max+2,
                       &result);
  if (iBase_SUCCESS != result) {
    std::cerr << "Trouble getting bounding box of gentity." << std::endl;
  }
}

void CAMALGeomEval::move_to(double& x, double& y, double& z) 
{
  int result;
  iGeom_getEntClosestPt(geomIface, myEnt, x, y, z, &x, &y, &z, &result);
  if (iBase_SUCCESS != result) {
    std::cerr << "Trouble getting closest point on gentity." << std::endl;
  }
}

bool CAMALGeomEval::normal_at(const double& x, const double& y, const double& z,
                              double& nx, double& ny, double& nz) 
{
  int result;
  iGeom_getEntNrmlXYZ(geomIface, myEnt, x, y, z, &nx, &ny, &nz, &result);
  if (iBase_SUCCESS != result) {
    std::cerr << "Trouble getting normal on gentity." << std::endl;
    return false;
  }

  return true;
}

#else

bool debug_surf_eval = false;

void CAMALGeomEval::set_mesh_size(double tmp_size) 
{
  if (debug_surf_eval) {
   std::cout << "set_mesh_size called." << std::endl;
  }
  meshSize = tmp_size;
}

double CAMALGeomEval::get_mesh_size() 
{
  if (debug_surf_eval) {
   std::cout << "get_mesh_size called." << std::endl;
  }
  return meshSize;
}

CAMALGeomEval::CAMALGeomEval(iGeom_Instance geom_iface, iBase_EntityHandle gent) 
  : geomIface( geom_iface ),
    myEnt( gent ),
    myDimension( -1 )
{
  int gtype, result;
  iGeom_getEntType(geomIface, gent, &gtype, &result);
  if (iBase_SUCCESS != result) {
    std::cerr << "Trouble getting dimension of gentity." << std::endl;
  }
  else 
    myDimension = gtype;

  if (debug_surf_eval) {
   std::cout << "CAMALGeomEval constructor called." << std::endl;
  }
}

CAMALGeomEval::~CAMALGeomEval() {}

double CAMALGeomEval::area() 
{
    // find the area of this entity
  if (myDimension != 2) 
    return -1.0;
  
  double measure = -1.0;
  double *measure_ptr = &measure;
  int measure_size, measure_alloc = 1;
  int result;
  iGeom_measure(geomIface, &myEnt, 1, &measure_ptr, &measure_alloc, &measure_size, &result);
  if (iBase_SUCCESS != result) {
    std::cerr << "Trouble getting area of gentity." << std::endl;
    return -1.0;
  }

  if (debug_surf_eval) {
    std::cout << "myEnt=" << myEnt << ",area=" << measure << std::endl;
  }
  
  return measure;
}

void CAMALGeomEval::bounding_box(double box_min[3], double box_max[3]) 
{
  int result;
  iGeom_getEntBoundBox(geomIface, myEnt,  
                       box_min, box_min+1, box_min+2,
                       box_max, box_max+1, box_max+2,
                       &result);
  if (iBase_SUCCESS != result) {
    std::cerr << "Trouble getting bounding box of gentity." << std::endl;
  }

  if (debug_surf_eval) {
    std::cout << "myEnt=" << myEnt << ",box_min=" << box_min[0]
	      << "," << box_min[1] << "," << box_min[2] 
	      << ",box_max=" << box_max[0]
	      << "," << box_max[1] << "," << box_max[2]<< std::endl;
  }
}

void CAMALGeomEval::move_to_surface(double& x, double& y, double& z) 
{
  int result;
  double ori[3] = {x, y, z};
  iGeom_getEntClosestPt(geomIface, myEnt, x, y, z, &x, &y, &z, &result);
  if (iBase_SUCCESS != result) {
    std::cerr << "Trouble getting closest point on gentity." << std::endl;
  }
  if (debug_surf_eval) {
    std::cout << "myEnt=" << myEnt << " is moved from "
	      << ori[0] << "," << ori[1] << "," << ori[2] << " to "
	      << x << "," << y << "," << z << std::endl;
  }
}

void CAMALGeomEval::move_to_surface(double& x, double& y, double& z,
		     double& u_guess, double& v_guess)
{
  if (debug_surf_eval) {
   std::cout << "move_to_surface called." << std::endl;
  }
}

bool CAMALGeomEval::normal_at(double x, double y, double z, 
			      double& nx, double& ny, double& nz) 
{
  int result;
  iGeom_getEntNrmlXYZ(geomIface, myEnt, x, y, z, &nx, &ny, &nz, &result);
  if (iBase_SUCCESS != result) {
    std::cerr << "Trouble getting normal on gentity." << std::endl;
    return false;
  }
  if (debug_surf_eval) {
    std::cout << "normal: at "
	      << x << "," << y << "," << z 
	      << ",to " << nx << "," << ny << ","
	      << nz << std::endl;
  }

  return true;
}

bool CAMALGeomEval::normal_at(double x, double y, double z, 
			      double& u_guess, double& v_guess,
			      double& nx, double& ny, double& nz)
{
  if (debug_surf_eval) {
   std::cout << "normal_at called." << std::endl;
  }
  return true;
}

bool CAMALGeomEval::is_planar()
{
 if (debug_surf_eval) {
   std::cout << "is_planar called." << std::endl;
 } 
 return false;
}

bool CAMALGeomEval::is_parametric()
{
  if (debug_surf_eval) {
    std::cout << "is_parametric called." << std::endl;
  } 
  return false;
}

bool CAMALGeomEval::is_periodic_in_u(double& u_period)
{
  if (debug_surf_eval) {
    std::cout << "is_periodic_in_u called." << std::endl;
  }
  return false;
}

bool CAMALGeomEval::is_periodic_in_v(double& v_period)
{
  if (debug_surf_eval) {
    std::cout << "is_periodic_in_v called." << std::endl;
  }
  return false;
}

void CAMALGeomEval::get_param_range_u(double& u_low, double& u_high)
{
  if (debug_surf_eval) {
    std::cout << "is_param_range_u called." << std::endl;
  }
}

void CAMALGeomEval::get_param_range_v(double& v_low, double& v_high)
{
}

bool CAMALGeomEval::uv_from_position(double x, double y, double z, 
				     double& u, double& v)
  
{
  if (debug_surf_eval) {
    std::cout << "uv_from_position." << std::endl;
  }
  return true;
}

bool CAMALGeomEval::uv_from_position(double x, double y, double z, 
				     double& u, double& v,
				     double& cx, double& cy, double& cz)
{
  if (debug_surf_eval) {
    std::cout << "uv_from_position." << std::endl;
  }
  return true;
}

void CAMALGeomEval::position_from_uv(double u, double v, 
				     double& x, double& y, double& z)
{
  if (debug_surf_eval) {
    std::cout << "position_from_uv." << std::endl;
  }
}

void CAMALGeomEval::distortion_at_uv(double u, double v, 
				     double du[3], double dv[3])
{
  if (debug_surf_eval) {
    std::cout << "distortion_at_uv." << std::endl;
  }
}

bool CAMALGeomEval::pierce_surface_with_ray(double & x, double & y, double & z, double dir_x,
         double dir_y, double dir_z)
{
   // will use the iGeom stuff
   /* iGeom_getPntRayIntsct( iGeom_Instance,
                                 double x,
                                 double y,
                                 double z,
                                 double dir_x,
                                 double dir_y,
                                 double dir_z,
                                 iBase_EntityHandle** intersect_entity_handles,
                                 int* intersect_entity_handles_allocated,
                                 int* intersect_entity_hangles_size,
                                 int storage_order,
                                 double** intersect_coords,
                                 int* intersect_coords_allocated,
                                 int* intersect_coords_size,
                                 double** param_coords,
                                 int* param_coords_allocated,
                                 int* param_coords_size,
                                 int* err ); */
   /*
    * // march along given direction  ; it should be at most size / 3?
  */
   iBase_EntityHandle * intersect_entity_handles = NULL;
   int intersect_entity_handles_allocated = 0, intersect_entity_handles_size = 0;
   double * intersect_coords = NULL;
   int intersect_coords_allocated =0 ,  intersect_coords_size = 0;
   double * param_coords = NULL;
   int param_coords_allocated = 0, param_coords_size =0;
   int err=0;
   iGeom_getPntRayIntsct( geomIface,
            x, y, z,
            dir_x, dir_y, dir_z,
            &intersect_entity_handles, &intersect_entity_handles_allocated,
            &intersect_entity_handles_size, iBase_INTERLEAVED,
            &intersect_coords, &intersect_coords_allocated, &intersect_coords_size,
            &param_coords, &param_coords_allocated, &param_coords_size,
            &err );
   // get the first coordinate
   if (iBase_SUCCESS != err || intersect_entity_handles_size ==0)
   {
      std::cerr << "Trouble getting intersection of a ray with geometry" << std::endl;
             return false;
   }
   // consider only the first intersection point
   x=intersect_coords[0];
   y=intersect_coords[1];
   z=intersect_coords[2];

   free(intersect_entity_handles);
   free(intersect_coords);
   free(param_coords);
   return true;
}
#endif
