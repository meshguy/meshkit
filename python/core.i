/* Declarations for the core classes in MeshKit (mostly MKCore and friends) */

%include "std_vector.i"
%include "std_string.i"

%{
#include "meshkit/MKCore.hpp"
#include "meshkit/ModelEnt.hpp"
#include "meshkit/MeshOp.hpp"
#include "meshkit/MeshOpProxy.hpp"
#include "meshkit/SizingFunction.hpp"
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

%include "meshkit/Types.hpp"
%include "meshkit/MKGraph.hpp"
%include "meshkit/MKCore.hpp"
%include "meshkit/ModelEnt.hpp"
%include "meshkit/MeshOp.hpp"
%include "meshkit/MeshOpProxy.hpp"
%include "meshkit/Matrix.hpp"
%include "meshkit/SizingFunction.hpp"