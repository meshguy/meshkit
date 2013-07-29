/*!
\example extrudemesh.cpp

\section extrudemesh_cpp_title <pretty-name-of-this-file>

\subsection extrudemesh_cpp_inf Misc. Information
\author <your-name-here>
\date 7-15-2013
\bug <placeholder>
\warning <placeholder>

\subsection extrudemesh_cpp_goal Goal

\subsection extrudemesh_cpp_cw Code Walkthrough

\subsection extrudemesh_cpp_in Input
\image html extrudemesh.in.jpg
There is no input.

\subsection extrudemesh_cpp_out Output
\image html extrudemesh.out.jpg

\subsection extrudemesh_cpp_src Source Code
*/
#include "meshkit/MKCore.hpp"
#include "meshkit/MeshOp.hpp"
#include "meshkit/ExtrudeMesh.hpp"
#include "meshkit/ModelEnt.hpp"

using namespace MeshKit;


#define DEFAULT_TEST_FILE "cube.cub"

MKCore *mk;

void test_load_and_extrude();
void test_extrude_quad();
void test_extrude_tri();

void help_test_extrude_face(iMesh_EntityTopology topo, double *coords,
                            size_t nverts);

int main(int argc, char **argv)
{
  mk = new MKCore();
  int num_fail = 0;

  test_load_and_extrude();
  test_extrude_quad();
  test_extrude_tri();

  delete mk;
  return num_fail;
}

void clear_mesh(MKCore *mk)
{
  iMesh *mesh = mk->imesh_instance();
  std::vector<iMesh::EntityHandle> ents;
  mesh->getEntities(mesh->getRootSet(), iBase_ALL_TYPES, iMesh_ALL_TOPOLOGIES,
                    ents);
  mesh->deleteEntArr(&ents[0], ents.size());
}

void test_load_and_extrude()
{
  clear_mesh(mk);
  std::string filename = TestDir + "/" + DEFAULT_TEST_FILE;
  mk->load_mesh(filename.c_str());

    // populate mesh with no attached geometry or relations
  mk->populate_model_ents(-1, 0, -1);

  // get the hexes
  MEntVector vols;
  mk->get_entities_by_dimension(2, vols);

  ExtrudeMesh *em = (ExtrudeMesh*) mk->construct_meshop("ExtrudeMesh", vols);
  em->set_name("extrude_mesh");

  // some entity tag types are always copy or expand
  em->expand_sets().add_tag("MATERIAL_SET");
  em->expand_sets().add_tag("DIRICHLET_SET");
  em->expand_sets().add_tag("NEUMANN_SET");

  Vector<3> dx; dx[0] = 10; dx[1] = 0; dx[2] = 0;
  em->set_transform(Extrude::Translate(dx, 5));

  // put them in the graph
  mk->get_graph().addArc(mk->root_node()->get_node(), em->get_node());
  mk->get_graph().addArc(em->get_node(), mk->leaf_node()->get_node());

  // mesh embedded boundary mesh, by calling execute
  mk->setup_and_execute();

  delete em;
}


void help_test_extrude_face(iMesh_EntityTopology topo, double *coords,
                            size_t nverts)
{
  clear_mesh(mk);

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

  std::vector<iMesh::EntityHandle> new_verts;
  std::vector<double> new_coords(nverts*3);
  std::vector<iMesh::EntityHandle> new_face;
  iMesh::EntitySetHandle copy_set;

  mesh->getEntSetEHData(set, em->copy_tag(), (iMesh::EntityHandle&)copy_set);

  mesh->getEntities(copy_set, iBase_ALL_TYPES, iMesh_ALL_TOPOLOGIES, new_face);
  CHECK_EQUAL(new_face.size(), size_t(1));

  mesh->getEntAdj(new_face[0], iBase_VERTEX, new_verts);
  CHECK_EQUAL(new_verts.size(), nverts);

  mesh->getVtxArrCoords(&new_verts[0], new_verts.size(), iBase_INTERLEAVED,
                        &new_coords[0]);

  for(size_t i=0; i<nverts; i++) {
    CHECK_REAL_EQUAL(coords[i*3+0],    new_coords[i*3+0], 0.00001);
    CHECK_REAL_EQUAL(coords[i*3+1],    new_coords[i*3+1], 0.00001);
    CHECK_REAL_EQUAL(coords[i*3+2]+10, new_coords[i*3+2], 0.00001);
  }

  iMesh::EntitySetHandle extrude_set;
  int n;
  mesh->getEntSetEHData(set, em->extrude_tag(),
                        (iMesh::EntityHandle&)extrude_set);
  mesh->getNumOfType(extrude_set, iBase_ALL_TYPES, n);
  CHECK_EQUAL(n, 10);


  delete em;
}

void test_extrude_quad()
{
  double coords[] = {
    0, 0, 0,
    1, 0, 0,
    1, 1, 0,
    0, 1, 0
  };

  help_test_extrude_face(iMesh_QUADRILATERAL, coords, 4);
}

void test_extrude_tri()
{
  double coords[] = {
    0, 0, 0,
    1, 0, 0,
    1, 1, 0
  };

  help_test_extrude_face(iMesh_TRIANGLE, coords, 3);
}
