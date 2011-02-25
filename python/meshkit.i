%module MeshKit

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
#include "meshkit/MKCore.hpp"
//#include "meshkit/ModelEnt.hpp"
#include "meshkit/MeshOp.hpp"
#include "meshkit/MeshOpProxy.hpp"
%}

%include "meshkit/Types.hpp"
%include "meshkit/MKCore.hpp"
 //%include "meshkit/ModelEnt.hpp"
 //%include "meshkit/Error.hpp"
%include "meshkit/MeshOp.hpp"
%include "meshkit/MeshOpProxy.hpp"
