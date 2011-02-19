#include "meshkit/QuadMesher.hpp"
#include "meshkit/MKCore.hpp"
#include "meshkit/ModelEnt.hpp"
#include "meshkit/RegisterMeshOp.hpp"
#include "moab/ReadUtilIface.hpp"
#include <vector>
#include <math.h>

namespace MeshKit
{
moab::EntityType QuadMesher_tps[] = {moab::MBVERTEX, moab::MBQUAD, moab::MBMAXTYPE};

const moab::EntityType* QuadMesher::output_types()
  { return QuadMesher_tps; }

//---------------------------------------------------------------------------//
// Construction Function for Edge Mesher
QuadMesher::QuadMesher(MKCore *mk_core, const MEntVector &me_vec) 
        : MeshScheme(mk_core, me_vec), schemeType(EQUAL)
{
   imesh = mk_core()->imesh_instance();
}

// setup function: set up the number of intervals for edge meshing through the 
// sizing function
void QuadMesher::setup_this()
{
   setup_boundary();
}

Jaal::Mesh* QuadMesher :: tri_quad_conversion (Jaal::Mesh *trimesh)
{
  Jaal::Tri2Quads t2quad;

  cout << "Input: Triangle Mesh " << endl;
  cout << "# Nodes " << trimesh->getSize(0) << endl;
  cout << "# Faces " << trimesh->getSize(2) << endl;

  Jaal::Mesh *quadmesh = t2quad.getQuadMesh(trimesh, 1);

  cout << "Input: Quad Mesh " << endl;
  cout << "# Nodes " << quadmesh->getSize(0) << endl;
  cout << "# Faces " << quadmesh->getSize(2) << endl;

  return quadmesh;
}

Jaal::Mesh* tri_quad_conversion (iMesh_Instance imesh)
{
  Jaal::JaalMoabConverter meshconverter;
  Jaal::Mesh *trimesh  = meshconverter.fromMOAB(imesh);
  Jaal::Mesh* quadmesh = tri_quad_conversion ( trimesh );

  meshconverter.toMOAB(quadmesh, imesh);

  delete trimesh;
  return quadmesh;
}


void QuadMesher::execute_this()
{
   tri_quad_conversion(imesh);

/*
  std::vector<iMesh::EntityHandle> trifaces;
   int i = 0;
   for (MEntSelection::iterator mit = mentSelection.begin();
         mit != mentSelection.end(); mit++) {
      ModelEnt *me = mit->first;
      trifaces[i++] = reinterpret_cast<iBase_EntityHandle> (me->mesh_handle());
  }

    me->commit_mesh(mit->second, COMPLETE_MESH);	
  }
*/
}

