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

%typemap(in) (int argc, char *argv[]) {
  int i;
  if (!PyList_Check($input)) {
    PyErr_SetString(PyExc_ValueError, "Expecting a list");
    return NULL;
  }
  $1 = PyList_Size($input);
  $2 = (char **) malloc(($1+1)*sizeof(char *));
  for (i = 0; i < $1; i++) {
    PyObject *s = PyList_GetItem($input,i);
    if (!PyString_Check(s)) {
        free($2);
        PyErr_SetString(PyExc_ValueError, "List items must be strings");
        return NULL;
    }
    $2[i] = PyString_AsString(s);
  }
  $2[i] = 0;
}

%typemap(freearg) (int argc, char *argv[]) {
   if ($2) free($2);
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
