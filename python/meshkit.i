/* Main MeshKit bindings. This just includes other SWIG files in the appropriate
   magical order to make everything work. */

%module MeshKit

%include "algs_factory.i"
%include "itaps.i"
%include "core.i"
%include "algs.i"
