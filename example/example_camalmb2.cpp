/*!
\example example_camalmb2.cpp

\section mb2_cpp_title <pretty-name-of-this-file>

\subsection mb2_cpp_in Input
\image html mb2.in.jpg
There is no input.

\subsection mb2_cpp_out Output
\image html mb2.out.jpg

\subsection mb2_cpp_inf Misc. Information
\author <your-name-here>
\date 7-15-2013
\bug <placeholder>
\warning <placeholder>

\subsection mb2_cpp_src Source Code
*/

#include "meshkit/MKCore.hpp"
#include "meshkit/SizingFunction.hpp"
#include "meshkit/ModelEnt.hpp"

using namespace MeshKit;

MKCore * mk;
bool save_mesh = false;
bool quadMesh = true;
double mesh_size = 0.3;
std::string file_name; //="shell.h5m";
std::string output_file; //="output.h5m";



void meshFB2();

int main(int argc, char **argv)
{
  // start up a new iGeom engine, based on moab
  /*moab::Interface * mb = new moab::Core();
   iGeom * fbiGeom = new FBiGeom(mb); // true means smooth
   MKCore mk(fbiGeom, mb); // iMesh, iRel, will be constructed*/

  file_name = TestDir + "/" + "shell.h5m";// default test file
  output_file = "output.h5m";

  if (argc < 5) {
    {
      std::cout << "Usage : filename  output  < q, t >  <mesh_size> \n";
      std::cout << "use default options: " << file_name << " " << output_file <<
          (quadMesh ? " q " : " t " ) << mesh_size << "\n";
    }
  }
  else
  {
    file_name = argv[1];
    output_file = argv[2];
    if (!strcmp(argv[3], "t"))
      quadMesh = false;
    mesh_size = atof(argv[4]);
    save_mesh = true;
  }

  int num_fail = 0;
  meshFB2();

  return num_fail;
}
void meshFB2()
{
  mk = new MKCore;

  mk->load_mesh(file_name.c_str(), NULL, 0, 0, 0, false, false);
  int indx=  mk->initialize_mesh_based_geometry(0);

  moab::Range tris;
  moab::ErrorCode rval = mk->moab_instance()->get_entities_by_dimension(
      (moab::EntityHandle)0 , 2,
      tris);

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

  mk->remove_mesh_based_geometry(indx);
  if (save_mesh) {
    // output mesh for surfaces (implicitly for edges too, as edges are children of surface sets)
    mk->save_mesh_from_model_ents(output_file.c_str(), surfs);
  }
  // delete the model ents too, unload geometry

  delete mk;
}
