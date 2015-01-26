#include "CopyGeom.hpp"
#include "iGeom.h"
#include <iostream>
#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <functional>

CopyGeom::CopyGeom(iGeom_Instance impl) : igeomImpl(impl)
{}

CopyGeom::~CopyGeom()
{}

// copy move entities to dx
void CopyGeom::copy(iBase_EntityHandle *entities, const int entities_ehsize,
		    const double *dx,
		    iBase_EntityHandle **new_ents,
		    int *new_ents_allocated,
		    int *new_ents_size,
		    bool do_merge)
{
  int err=0;

  for (int i =0; i < entities_ehsize; i++){

    iBase_EntityHandle temp;
    iGeom_copyEnt(igeomImpl, entities[i], &temp, &err);
    check_error(igeomImpl, err);
 
    iGeom_moveEnt(igeomImpl, temp, dx[0], dx[1], dx[2], &err);
    check_error(igeomImpl, err);
  }
}

void CopyGeom::tag_copied_sets(const char **tag_names, const char **tag_vals,
                       const int num_tags)
{
  // TODO:

}
