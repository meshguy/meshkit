/* Factory code for automatic downcasting of algorithm objects. When adding a
   new algorithm, add a reference here and in algs.i! */

%include "factory.i"

%{
#include "meshkit/Transform.hpp"
#include "meshkit/CESets.hpp"
#include "meshkit/QslimOptions.hpp"

#include "meshkit/CopyGeom.hpp"
#include "meshkit/CopyMesh.hpp"
#include "meshkit/EBMesher.hpp"
#include "meshkit/EdgeMesher.hpp"
#include "meshkit/ExtrudeMesh.hpp"
#include "meshkit/OneToOneSwept.hpp"
#include "meshkit/QslimMesher.hpp"
#include "meshkit/SCDMesh.hpp"
#include "meshkit/TFIMapping.hpp"
#include "meshkit/VertexMesher.hpp"

#include "meshkit/CAMALPaver.hpp"
#include "meshkit/CAMALTetMesher.hpp"
#include "meshkit/CAMALTriAdvance.hpp"
%}

%factory(MeshKit::MeshOp * MeshKit::MKCore::construct_meshop,
         MeshKit::CopyGeom,
         MeshKit::CopyMesh,
         MeshKit::EBMesher,
         MeshKit::EdgeMesher,
         MeshKit::ExtrudeMesh,
         MeshKit::OneToOneSwept,
         MeshKit::QslimMesher,
         MeshKit::SCDMesh,
         MeshKit::TFIMapping,
         MeshKit::VertexMesher,

         MeshKit::CAMALPaver,
         MeshKit::CAMALTetMesher,
         MeshKit::CAMALTriAdvance);
