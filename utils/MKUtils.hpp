#ifndef MESHKIT_MKUTILS_HPP
#define MESHKIT_MKUTILS_HPP

#include "iMesh.h"

class MKUtils {

public:
  MKUtils(iMesh_Instance impl) : imeshImpl(impl) {};
  
  ~MKUtils() {};

  int assign_global_ids(iBase_EntitySetHandle this_set,
                        const int dimension, 
                        const int start_id,
                        const bool largest_dim_only,
                        const bool parallel,
                        const char *tag_name = "GLOBAL_ID");

private:

  iMesh_Instance imeshImpl;
  
};

#endif
