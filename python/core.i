/* Declarations for the core classes in MeshKit (mostly MKCore and friends) */

%include "std_vector.i"
%include "std_string.i"

%{
#include "meshkit/MKCore.hpp"

#include "meshkit/GraphNode.hpp"
#include "meshkit/ModelEnt.hpp"
#include "meshkit/MeshOp.hpp"
#include "meshkit/MeshOpProxy.hpp"
#include "meshkit/MeshScheme.hpp"
#include "meshkit/SizingFunction.hpp"
#include "meshkit/SizingFunctionVar.hpp"
%}

/* MeshKit::Vector typemaps */
%typemap(in) MeshKit::Vector<3> & (MeshKit::Vector<3> temp) {
  if (PySequence_Size($input) == 3) {
    for (size_t i=0; i<3; i++)
      temp[i] = PyFloat_AsDouble(PySequence_GetItem($input, i));
    $1 = &temp;
  }
}

%template(MEntVector) std::vector<MeshKit::ModelEnt*>;

%exception {
  try {
    $function
  }
  catch (const MeshKit::Error &e) {
    PyErr_SetString(PyExc_RuntimeError, e.what());
    return NULL;
  }
}

%include "meshkit/Types.hpp"
%include "meshkit/MKGraph.hpp"
%include "meshkit/MKCore.hpp"

%include "meshkit/GraphNode.hpp"
%include "meshkit/ModelEnt.hpp"
%include "meshkit/MeshOp.hpp"
%include "meshkit/MeshOpProxy.hpp"
%include "meshkit/MeshScheme.hpp"
%include "meshkit/Matrix.hpp"
%include "meshkit/SizingFunction.hpp"
%include "meshkit/SizingFunctionVar.hpp"
