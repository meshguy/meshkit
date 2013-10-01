/*!
\example example_copymesh.cpp

\section example_CopyMesh_cpp_title CopyMesh With Multiple Files

\subsection example_CopyMesh_cpp_in Input
Two input files c.exo and d.exo from data folder. Translate coordinates dx.
\subsection example_CopyMesh_cpp_out Output
Final mesh after CopyMesh operation.
\subsection example_CopyMesh_cpp_inf Misc. Information
\date 9-30-2013

\subsection example_CopyMesh_cpp_src Source Code
*/
#include "meshkit/MKCore.hpp"
#include "meshkit/MeshOp.hpp"
#include "meshkit/CopyMesh.hpp"
#include "meshkit/ModelEnt.hpp"

using namespace MeshKit;

#define DEFAULT_TEST_FILE1 "c.exo"
#define DEFAULT_TEST_FILE2 "d.exo"

MKCore *mk;

int main(int argc, char **argv)
{
  mk = new MKCore();
  std::vector <CopyMesh*> cm;
  iMesh *imesh = mk->imesh_instance();
  int err =0;
  std::vector<std::string> files;
  std::string name1 = std::string(MESH_DIR) + "/" + DEFAULT_TEST_FILE1;
  std::string name2 = std::string(MESH_DIR) + "/" + DEFAULT_TEST_FILE2;

  files.push_back(name1);
  files.push_back(name2);
  cm.resize(files.size());
  iBase_EntitySetHandle orig_set;

  for (unsigned int i = 0; i < files.size(); i++) {
      iMesh_createEntSet(imesh->instance(), 0, &orig_set, &err);
      iMesh_load(imesh->instance(), orig_set, files[i].c_str(), NULL, &err, strlen(files[i].c_str()), 0);
      ModelEnt *me;
      me = NULL;
      me = new ModelEnt(mk, iBase_EntitySetHandle(0), /*igeom instance*/0, (moab::EntityHandle)orig_set, 0);
      MEntVector assm_set;
      assm_set.push_back(me);
      cm[i] = (CopyMesh*) mk->construct_meshop("CopyMesh", assm_set);
      cm[i]->set_name("copy_move_mesh");
      cm[i]->copy_sets().add_set(orig_set);
      // some entity tag types are always copy or expand
      cm[i]->expand_sets().add_tag("MATERIAL_SET");
      cm[i]->expand_sets().add_tag("DIRICHLET_SET");
      cm[i]->expand_sets().add_tag("NEUMANN_SET");
      cm[i]->copy_sets().add_tag("GEOM_DIMENSION");
      Vector<3> dx; dx[0] = 10; dx[1] = 0; dx[2] = 0;
      cm[i]->set_transform(Copy::Translate(dx));
    }
  // copy and merge
  mk->setup_and_execute();
  for (unsigned int i = 0; i < files.size(); i++)
    delete cm[i];
  mk->save_mesh("o.h5m");
  return 0;
}
