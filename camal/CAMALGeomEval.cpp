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
#include "CAMALGeomEval.hpp"
#include <iostream>

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
