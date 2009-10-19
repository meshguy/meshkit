#include <iostream>
#include "MKUtils.hpp"

#define ERRORR(a,b) {if (iBase_SUCCESS != err) {std::cerr << a << std::endl; return b;}}

// HJK : simple main function
// need more implementation
int main()
{
  int err;
  iMesh_Instance mesh;
  iMesh_newMesh("", &mesh, &err, 0);
  ERRORR("Couldn't create mesh.", 1);

  iBase_EntitySetHandle root_set;
  iMesh_getRootSet(mesh, &root_set, &err);
  ERRORR("Couldn't get root set.", 1);

  MKUtils mkutils(mesh);
  
  err = mkutils.assign_global_ids(root_set, 3, 1, true, false,
				  "GLOBAL_ID");
  ERRORR("Error assigning global ids.", err);

  return 0;
}

