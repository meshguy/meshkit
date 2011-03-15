%module MeshKit
%include "std_vector.i"
%include "std_string.i"
%include "factory.i"

%rename(ExtrudeTransform) MeshKit::Extrude::Transform;
%rename(ExtrudeTranslate) MeshKit::Extrude::Translate;
%rename(ExtrudeRotate) MeshKit::Extrude::Rotate;

%{
#include "python/PyTAPS/iMesh_Python.h"

#include "meshkit/iMesh.hpp"

#include "meshkit/MKCore.hpp"
#include "meshkit/ModelEnt.hpp"
#include "meshkit/MeshOp.hpp"
#include "meshkit/MeshOpProxy.hpp"
#include "meshkit/Transform.hpp"

#include "meshkit/EdgeMesher.hpp"
#include "meshkit/VertexMesher.hpp"
%}

%factory(MeshKit::MeshOp * MeshKit::MKCore::construct_meshop,
         MeshKit::EdgeMesher,
         MeshKit::VertexMesher);

%template(MEntVector) std::vector<MeshKit::ModelEnt*>;

/* C ITAPS typemaps */
%typemap(in) iMesh_Instance {
    $1 = ((iMesh_Object*)$input)->handle;
}
%typemap(out) iMesh_Instance {
    $result = iMesh_FromInstance($1);
}

/* C++ ITAPS typemaps */
%typemap(out) iMesh* {
    $result = iMesh_FromInstance($1->instance());
}

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
%include "meshkit/Transform.hpp"

%include "meshkit/EdgeMesher.hpp"
%include "meshkit/VertexMesher.hpp"

%init {
    import_iBase();
    import_iMesh();
    import_array();
}
