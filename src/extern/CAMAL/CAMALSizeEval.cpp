#if CAMAL_VERSION > 500
#include "CAMALSizeEval.hpp"
#include <iostream>


CAMALSizeEval::CAMALSizeEval(double size)
  : meshSize(size)
{
}

CAMALSizeEval::~CAMALSizeEval() {}

bool CAMALSizeEval::size_at_point(double x, double y, double z, double& size, 
			   int level)
{
  size = meshSize;
  return true;
}

bool CAMALSizeEval::stretch_vector(double loc_x, double loc_y, double loc_z,
			    double& vec_x, double& vec_y, double& vec_z)
{
  vec_x = loc_x;
  vec_y = loc_y;
  vec_z = loc_z;
  return true;
}

bool CAMALSizeEval::tensor_at_point(double x, double y, double z, 
			     double*& size)
{
  for (int i = 0; i < 9; i++) {
    size[i] = meshSize;
  }
  return true;
}

bool CAMALSizeEval::size_at_point(double u, double v, double& size, 
			   int level)
{
  size = meshSize;
  return true;
}

bool CAMALSizeEval::stretch_vector(double loc_u, double loc_v,
			    double& vec_u, double& vec_v)
{
  vec_u = loc_u;
  vec_v = loc_v;
  return true;
}

bool CAMALSizeEval::tensor_at_point(double u, double v, 
			     double*& size)
{
  for (int i = 0; i < 4; i++) {
    size[i] = meshSize;
  }
}

bool CAMALSizeEval::size_at_point(double u, double& size, 
			   int level)
{
  size = meshSize;
  return true;
}

bool CAMALSizeEval::is_anisotropic()
{
  return false;
}

#endif
