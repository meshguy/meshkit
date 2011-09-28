/** \file test_mbgeom.cpp
 *
 * Test the mesh based geometry in meshkit
 *
 */

#include "meshkit/MKCore.hpp"
#include "meshkit/FBiGeom.hpp"
#include "meshkit/SizingFunctionVar.hpp"
#include "meshkit/ModelEnt.hpp"

using namespace MeshKit;

MKCore * mk;
bool save_mesh = false;
bool quadMesh = true;
double mesh_size = 0.3;
std::string file_name; //="shell.h5m";
double fixedPoint[3] = {0., 0., 0.};
double a=0.1, b=0.1, c=0.1; // so size will vary from 0.3 to about 0.6, for the default model

#include "TestUtil.hpp"

void meshFBvar();

int main(int argc, char **argv)
{
  // start up a new iGeom engine, based on moab
  /*moab::Interface * mb = new moab::Core();
   iGeom * fbiGeom = new FBiGeom(mb); // true means smooth
   MKCore mk(fbiGeom, mb); // iMesh, iRel, will be constructed*/

  file_name = TestDir + "/" + "shell.h5m";// default test file
  if (argc < 4) {
    {
      std::cout << "Usage : filename < q, t >  <mesh_size>  x0  y0  z0  a  b  c \n";
      std::cout <<  " mesh size is variable with formula: \n size(x, y, z) = mesh_size+a(x-x0)+b(y-y0)+c(z-z0)\n";
      std::cout << "use default options: " << file_name << " " <<
          (quadMesh ? "q " : "t " ) << mesh_size ;
      std::cout << " 0. 0. 0. 0.1 0.1 0.1 \n";
    }
  }
  else
  {
    file_name = argv[1];
    if (!strcmp(argv[2], "t"))
      quadMesh = false;
    mesh_size = atof(argv[3]);
    save_mesh = true;
    // give some linear coefficients to the var size case
    if (argc==10)
    {
      fixedPoint[0] = atof(argv[4]);
      fixedPoint[1] = atof(argv[5]);
      fixedPoint[2] = atof(argv[6]);
      a = atof(argv[7]);
      b = atof(argv[8]);
      c = atof(argv[9]);
    }
    else
    {
      std::cout<<" wrong number of arguments: argc=" << argc << "\n";
      return 1; // fail;
    }

  }

  int num_fail = 0;
  num_fail += RUN_TEST(meshFBvar);

  return num_fail;
}

void meshFBvar()
{
  mk = new MKCore;
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



  // size the mesh
  SizingFunctionVar *sf = new SizingFunctionVar(mk, -1, mesh_size);
  for (unsigned int i = 0; i < surfs.size(); i++)
    surfs[i]->sizing_function_index(sf->core_index());

  double coeff[4] = {a, b, c, mesh_size};
  sf->set_linear_coeff(fixedPoint, coeff);

  // now mesh them
  mk->setup_and_execute();

  // report the number of triangles generated
  // put it in a new set in moab

  rval = mk->moab_instance()->get_entities_by_dimension(0, 2, tris);
  CHECK_EQUAL(moab::MB_SUCCESS, rval);
  std::string elem =(quadMesh? " quads ": " triangles");
  std::cout << tris.size() - nbInitial << elem << " generated." << std::endl;

  delete fbiGeom;// this will trigger also deletion of FBEngine tags, etc

  if (save_mesh) {
    // output mesh for surfaces (implicitly for edges too, as edges are children of surface sets)
    std::string outfile = file_name + std::string(".var.h5m");
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
