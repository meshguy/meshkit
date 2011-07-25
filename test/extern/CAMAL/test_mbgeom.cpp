/** \file test_mbgeom.cpp
 *
 * Test the mesh based geometry in meshkit
 *
 */

#include "meshkit/MKCore.hpp"
#include "meshkit/FBiGeom.hpp"
#include "meshkit/SizingFunction.hpp"
#include "meshkit/ModelEnt.hpp"

using namespace MeshKit;

MKCore * mk;
bool save_mesh = false;
bool quadMesh = true;
double mesh_size = 0.3;
std::string file_name; //="shell.h5m";

#include "TestUtil.hpp"

void meshFB();

int main(int argc, char **argv)
{
  // start up a new iGeom engine, based on moab
  /*moab::Interface * mb = new moab::Core();
   iGeom * fbiGeom = new FBiGeom(mb); // true means smooth
   MKCore mk(fbiGeom, mb); // iMesh, iRel, will be constructed*/

  file_name = TestDir + "/" + "shell.h5m";// default test file
  if (argc < 4) {
    {
      std::cout << "Usage : filename < q, t >  <mesh_size> \n";
      std::cout << "use default options: " << file_name << " " <<
          (quadMesh ? "q " : "t " ) << mesh_size << "\n";
    }
  }
  else
  {
    file_name = argv[1];
    if (!strcmp(argv[2], "t"))
      quadMesh = false;
    mesh_size = atof(argv[3]);
    save_mesh = true;
  }

  int num_fail = 0;
  num_fail += RUN_TEST(meshFB);

  return num_fail;
}
void meshFB()
{
  mk = new MKCore;
  // just for debugging
  mk->print_graph();
  FBiGeom * fbiGeom = new FBiGeom(mk, true); // true for smooth, false for linear

  // this will do the reading of the moab db in memory

  //this also will do heavy stuff, like smoothing
  // we do not want to do it in the constructor
  // this should also populate ModelEnts in MKCore
  fbiGeom->load(file_name.c_str());

  moab::Range tris;
  moab::ErrorCode rval = mk->moab_instance()->get_entities_by_dimension(
      (moab::EntityHandle)fbiGeom->getRootSet() , 2,
      tris, /*recursive*/ true);

  int nbInitial = tris.size();
  tris.clear();
  // initial number of triangles
  std::cout << nbInitial << " initial elements in the model" << std::endl;
  //fbiGeom->Init();

  // get the surface model ents
  MEntVector surfs;
  mk->get_entities_by_dimension(2, surfs);

  // create surface mesher for them, either quad or tri mesher
  if (quadMesh)
    mk->construct_meshop("CAMALPaver", surfs);
  else
    mk->construct_meshop("CAMALTriAdvance", surfs);


  // just for debugging
  mk->print_graph();

  // size the mesh
  SizingFunction *sf = new SizingFunction(mk, -1, mesh_size);
  for (unsigned int i = 0; i < surfs.size(); i++)
    surfs[i]->sizing_function_index(sf->core_index());

  // now mesh them
  mk->setup_and_execute();

  std::cout<<" after execute:\n";
  // just for debugging
  mk->print_graph();

  // report the number of triangles generated
  // put it in a new set in moab

  rval = mk->moab_instance()->get_entities_by_dimension(0, 2, tris);
  CHECK_EQUAL(moab::MB_SUCCESS, rval);
  std::string elem =(quadMesh? " quads ": " triangles");
  std::cout << tris.size() - nbInitial << elem << " generated." << std::endl;

  delete fbiGeom;// this will trigger also deletion of FBEngine tags, etc

  if (save_mesh) {
    // output mesh for surfaces (implicitly for edges too, as edges are children of surface sets)
    std::string outfile = file_name + std::string(".h5m");
    moab::EntityHandle out_set;
    rval = mk->moab_instance()->create_meshset(moab::MESHSET_SET, out_set);
    MBERRCHK(rval, mk->moab_instance());
    for (unsigned int i = 0; i < surfs.size(); i++) {
      rval = mk->moab_instance()->add_child_meshset(out_set,
          surfs[i]->mesh_handle());
      MBERRCHK(rval, mk->moab_instance());
    }

    rval = mk->moab_instance()->write_file(outfile.c_str(), NULL, NULL,
        &out_set, 1);
    MBERRCHK(rval, mk->moab_instance());
  }
  // delete the model ents too, unload geometry

  delete mk;
}
