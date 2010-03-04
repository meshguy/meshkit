#ifndef CAMAL_SIZE_EVAL_HPP
#define CAMAL_SIZE_EVAL_HPP

#if CAMAL_VERSION > 500
#include "CMLSizeEval.hpp"
#include "iGeom.h"

class CAMALSizeEval : public CMLSizeEval
{
public:
  CAMALSizeEval(double size);
  
  virtual ~CAMALSizeEval();
  
  virtual bool size_at_point(double x, double y, double z, double& size, 
                             int level = -1);

  virtual bool stretch_vector(double loc_x, double loc_y, double loc_z,
                              double& vec_x, double& vec_y, double& vec_z);

  virtual bool tensor_at_point(double x, double y, double z, 
                               double*& size);

  virtual bool size_at_point(double u, double v, double& size, 
                             int level = -1);

  virtual bool stretch_vector(double loc_u, double loc_v,
			      double& vec_u, double& vec_v);

  virtual bool tensor_at_point(double u, double v, 
                               double*& size);

  virtual bool size_at_point(double u, double& size, 
                             int level = -1);

  virtual bool is_anisotropic();

private:
  double meshSize;
};

#endif
#endif
