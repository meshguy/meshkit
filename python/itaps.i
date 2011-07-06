/* Bindings to ITAPS interfaces by way of PyTAPS */

%{
#include "python/PyTAPS/iGeom_Python.h"
#include "python/PyTAPS/iMesh_Python.h"
#include "python/PyTAPS/iRel_Python.h"

#include "meshkit/iGeom.hpp"
#include "meshkit/iMesh.hpp"
#include "meshkit/iRel.hpp"
%}

/* C ITAPS typemaps */
%typemap(in)  iGeom_Instance { $1 = ((iGeom_Object*)$input)->handle; }
%typemap(out) iGeom_Instance { $result = iGeom_FromInstance($1); }

%typemap(in)  iMesh_Instance { $1 = ((iMesh_Object*)$input)->handle; }
%typemap(out) iMesh_Instance { $result = iMesh_FromInstance($1); }

%typemap(in)  iRel_Instance  { $1 = ((iRel_Object*) $input)->handle; }
%typemap(out) iRel_Instance  { $result = iRel_FromInstance ($1); }

/* C++ ITAPS typemaps */
%typemap(out) iGeom* { $result = iGeom_FromInstance($1->instance()); }
%typemap(out) iMesh* { $result = iMesh_FromInstance($1->instance()); }
%typemap(out) iRel*  { $result = iRel_FromInstance ($1->instance()); }

%init {
  import_iBase();
  import_iGeom();
  import_iMesh();
  import_iRel();
  import_array();
}
