/* Factory code for automatic downcasting of algorithm objects. When adding a
   new algorithm, add a reference here and in algs.i! */

%include "factory.i"

%{
#include "meshkit/Transform.hpp"
#include "meshkit/CESets.hpp"

#include "meshkit/EdgeMesher.hpp"
#include "meshkit/VertexMesher.hpp"
#include "meshkit/CopyMesh.hpp"
#include "meshkit/OneToOneSwept.hpp"

#include "meshkit/CAMALPaver.hpp"
#include "meshkit/CAMALTetMesher.hpp"
#include "meshkit/CAMALTriAdvance.hpp"
%}

%factory(MeshKit::MeshOp * MeshKit::MKCore::construct_meshop,
         MeshKit::EdgeMesher,
         MeshKit::VertexMesher,
         MeshKit::CopyMesh,
         MeshKit::OneToOneSwept,
         MeshKit::CAMALPaver,
         MeshKit::CAMALTetMesher,
         MeshKit::CAMALTriAdvance);
