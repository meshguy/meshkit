%module MeshKit
%include "std_vector.i"
%include "std_string.i"

%rename(ExtrudeTransform) MeshKit::Extrude::Transform;
%rename(ExtrudeTranslate) MeshKit::Extrude::Translate;
%rename(ExtrudeRotate) MeshKit::Extrude::Rotate;

%typemap(in) MeshKit::Vector<3> {
    $1 = MeshKit::Vector<3>();
}

/* Convert from C --> Python */
%typemap(out) MeshKit::Vector<3> {
    $result = PyInt_FromLong(0);
}

%{
#include "python/PyTAPS/iMesh_Python.h"

#include "meshkit/MKCore.hpp"
#include "meshkit/ModelEnt.hpp"
#include "meshkit/MeshOp.hpp"
#include "meshkit/MeshOpProxy.hpp"
%}

%template(MEntVector) std::vector<MeshKit::ModelEnt*>;

/* Convert from Python --> C */
%typemap(in) iMesh_Instance {
    $1 = ((iMesh_Object*)$input)->handle;
}

/* Convert from C --> Python */
%typemap(out) iMesh_Instance {
    $result = iMesh_FromInstance($1);
}

%include "meshkit/Types.hpp"
%include "meshkit/MKCore.hpp"
%include "meshkit/ModelEnt.hpp"
%include "meshkit/MeshOp.hpp"
%include "meshkit/MeshOpProxy.hpp"
