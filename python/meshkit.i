%module MeshKit
%include "std_vector.i"
%include "std_string.i"
%include "factory.i"

%{
#include "python/PyTAPS/iGeom_Python.h"
#include "python/PyTAPS/iMesh_Python.h"
#include "python/PyTAPS/iRel_Python.h"

#include "meshkit/iGeom.hpp"
#include "meshkit/iMesh.hpp"
#include "meshkit/iRel.hpp"

#include "meshkit/MKCore.hpp"
#include "meshkit/ModelEnt.hpp"
#include "meshkit/MeshOp.hpp"
#include "meshkit/MeshOpProxy.hpp"

#include "meshkit/Transform.hpp"
#include "meshkit/CESets.hpp"

#include "meshkit/EdgeMesher.hpp"
#include "meshkit/VertexMesher.hpp"
#include "meshkit/CopyMesh.hpp"
%}

%factory(MeshKit::MeshOp * MeshKit::MKCore::construct_meshop,
         MeshKit::EdgeMesher,
         MeshKit::VertexMesher,
         MeshKit::CopyMesh);

%template(MEntVector) std::vector<MeshKit::ModelEnt*>;

/* C ITAPS typemaps */
%typemap(in)  iGeom_Instance { $1 = ((iGeom_Object*)$input)->handle; }
%typemap(out) iGeom_Instance { $result = iGeom_FromInstance($1); }

%typemap(in)  iMesh_Instance { $1 = ((iMesh_Object*)$input)->handle; }
%typemap(out) iMesh_Instance { $result = iMesh_FromInstance($1); }

%typemap(in)  iRel_Instance  { $1 = ((iRel_Object*)$input)->handle; }
%typemap(out) iRel_Instance  { $result = iRel_FromInstance($1); }

/* C++ ITAPS typemaps */
%typemap(out) iGeom* { $result = iGeom_FromInstance($1->instance()); }
%typemap(out) iMesh* { $result = iMesh_FromInstance($1->instance()); }
%typemap(out) iRel*  { $result = iRel_FromInstance ($1->instance()); }

/* MeshKit::Vector typemaps */
%typemap(in) MeshKit::Vector<3> & (MeshKit::Vector<3> temp) {
  if (PySequence_Size($input) == 3) {
    for (size_t i=0; i<3; i++)
      temp[i] = PyFloat_AsDouble(PySequence_GetItem($input, i));
    $1 = &temp;
  }
}

%include "meshkit/Types.hpp"
%include "meshkit/MKCore.hpp"
%include "meshkit/ModelEnt.hpp"
%include "meshkit/MeshOp.hpp"
%include "meshkit/MeshOpProxy.hpp"
%include "meshkit/Matrix.hpp"

%rename(Copy_Transform) MeshKit::Copy::Transform;
%rename(Copy_Identity) MeshKit::Copy::Identity;
%rename(Copy_Translate) MeshKit::Copy::Translate;
%rename(Copy_Rotate) MeshKit::Copy::Rotate;

%rename(Extrude_Transform) MeshKit::Extrude::Transform;
%rename(Extrude_Translate) MeshKit::Extrude::Translate;
%rename(Extrude_Rotate) MeshKit::Extrude::Rotate;

%include "meshkit/TransformBase.hpp"

%template(Copy_BasicTransform_Identity)
  MeshKit::Copy::BasicTransform<MeshKit::Copy::Identity>;
%template(Copy_BasicTransform_Translate)
  MeshKit::Copy::BasicTransform<MeshKit::Copy::Translate>;
%template(Copy_BasicTransform_Rotate)
  MeshKit::Copy::BasicTransform<MeshKit::Copy::Rotate>;

%template(Extrude_BasicTransform_Translate)
  MeshKit::Extrude::BasicTransform<MeshKit::Extrude::Translate>;
%template(Extrude_BasicTransform_Rotate)
  MeshKit::Extrude::BasicTransform<MeshKit::Extrude::Rotate>;

%include "meshkit/Transform.hpp"
%include "meshkit/CESets.hpp"

%include "meshkit/EdgeMesher.hpp"
%include "meshkit/VertexMesher.hpp"
%include "meshkit/CopyMesh.hpp"

%init {
    import_iBase();
    import_iGeom();
    import_iMesh();
    import_iRel();
    import_array();
}
