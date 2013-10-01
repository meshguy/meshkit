/*!
\example example_extrudemesh.cpp

\section example_extrudemesh_cpp_title 10 Hex ExtrudeMesh Example

\subsection example_extrudemesh_cpp_in Input
This example creates a face and extrudes it. \n
Inputs are: coords of face, sides in the face, topology of the face and translation distance dx.
\subsection example_extrudemesh_cpp_out Output
Final mesh with 10 hex elements is saved
\subsection example_extrudemesh_cpp_inf Misc. Information
\author <your-name-here>
\date 9-30-2013

\subsection example_extrudemesh_cpp_src Source Code
*/
#include "meshkit/MKCore.hpp"
#include "meshkit/MeshOp.hpp"
#include "meshkit/ExtrudeMesh.hpp"
#include "meshkit/ModelEnt.hpp"

using namespace MeshKit;


#define DEFAULT_TEST_FILE "cube.cub"

MKCore *mk;

int main(int argc, char **argv)
{
  mk = new MKCore();
  //input
  size_t nverts = 4;
  double coords[] = {
    0, 0, 0,
    1, 0, 0,
    1, 1, 0,
    0, 1, 0
  };
  iMesh_EntityTopology topo = iMesh_QUADRILATERAL;

  // working with imesh instance
  iMesh *mesh = mk->imesh_instance();
  std::vector<iMesh::EntityHandle> verts(nverts);
  iMesh::EntityHandle face;
  iMesh::EntitySetHandle set;

  mesh->createVtxArr(nverts, iBase_INTERLEAVED, coords, &verts[0]);
  mesh->createEnt(topo, &verts[0], nverts, face);
  mesh->createEntSet(true, set);
  mesh->addEntToSet(face, set);

  ModelEnt me(mk, iBase_EntitySetHandle(0), /*igeom instance*/0,
              (moab::EntityHandle)set);
  MEntVector selection;
  selection.push_back(&me);

  ExtrudeMesh *em = (ExtrudeMesh*) mk->construct_meshop("ExtrudeMesh",
                                                        selection);
  em->set_name("extrude_mesh");

  Vector<3> dx; dx[0] = 0; dx[1] = 0; dx[2] = 10;
  em->set_transform(Extrude::Translate(dx, 10));
  em->copy_faces(true);

  // put them in the graph
  mk->get_graph().addArc(mk->root_node()->get_node(), em->get_node());
  mk->get_graph().addArc(em->get_node(), mk->leaf_node()->get_node());

  em->copy_sets().add_set(set);
  em->extrude_sets().add_set(set);

  // mesh embedded boundary mesh, by calling execute
  mk->setup_and_execute();
  mk->save_mesh("out_em.exo");
  return 0;
}
