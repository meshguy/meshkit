/* Binding code for algorithms.  When adding a new algorithm, add a reference
   here and in algs_factory.i! */

%rename(Copy_Transform) MeshKit::Copy::Transform;
%rename(Copy_Identity)  MeshKit::Copy::Identity;
%rename(Copy_Translate) MeshKit::Copy::Translate;
%rename(Copy_Rotate)    MeshKit::Copy::Rotate;

%rename(Extrude_Transform) MeshKit::Extrude::Transform;
%rename(Extrude_Translate) MeshKit::Extrude::Translate;
%rename(Extrude_Rotate)    MeshKit::Extrude::Rotate;

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
%include "meshkit/QslimOptions.hpp"

// TODO: AssyGen
%include "meshkit/CopyGeom.hpp"
%include "meshkit/EBMesher.hpp"
%include "meshkit/EdgeMesher.hpp"
%include "meshkit/ExtrudeMesh.hpp"
// TODO: MergeMesh
%include "meshkit/OneToOneSwept.hpp"
// TODO: ParallelMesher
// TODO: ParExchangeMesh
%include "meshkit/QslimMesher.hpp"
%include "meshkit/SCDMesh.hpp"
%include "meshkit/TFIMapping.hpp"
%include "meshkit/VertexMesher.hpp"

%include "meshkit/CAMALPaver.hpp"
%include "meshkit/CAMALTetMesher.hpp"
%include "meshkit/CAMALTriAdvance.hpp"

%include "meshkit/TriangleMesher.hpp"
%include "meshkit/MBGeomOp.hpp"
