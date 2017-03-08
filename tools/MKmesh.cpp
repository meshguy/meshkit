/*
 * MKmesh.cpp
 *
 *  create a tool for meshop, to run different input files with various sizes
 *
 *  example of use
 *  MKmesh <input_geo_file> AF2DfltTriangleMeshOp <size> <out_mesh>  <dbg_level>
 */

// C++
#include <cstddef>
#include <iostream>
#include <string>

// MeshKit
#include "meshkit/MKCore.hpp"
#include "meshkit/iGeom.hpp"
#include "meshkit/MeshOp.hpp"
#include "meshkit/MeshOpProxy.hpp"
#include "meshkit/SizingFunction.hpp"
#include "meshkit/SizingFunctionVar.hpp"
#include "meshkit/ModelEnt.hpp"

// MeshKit testing
#include "TestUtil.hpp"


using namespace MeshKit;

int main(int argc, char **argv)
{
  // This variable is defined and used in main because a new MKCore
  // instance cannot be easily constructed after another MKCore
  // instance is deleted; there are problems with a tag left behind in
  // iGeom.
  MKCore * mk = new MeshKit::MKCore();

  std::string file_name, meshop, out_mesh;
  double size = 1.0;
  int debug=0, mesh_op_dim=3;

  if (argc < 6)
  {
    std::cout << " usage: " << argv[0] << " <geo_file> <mesh_op> <mesh_op_dim>[0, 1, 2, 3] <size> <out_file> [debug_level] \n ";
    // available mesh ops:
    std::cout << " Available mesh-ops:\n";
    const unsigned n = MKCore::num_meshops();
    for (unsigned i = 0; i < n; ++i) {
      MeshOpProxy* p = MKCore::meshop_proxy(i);
      std::cout << i << " " << p->name() << "\n";
    }
    return 1;
  }
  if (argc >= 6)
  {
    file_name = argv[1];
    meshop = argv[2];
    mesh_op_dim = atoi(argv[3]);
    size = atof(argv[4]);
    out_mesh = argv[5];
  }

  if (argc>=7)
  {
    debug = atoi(argv[6]);
  }

  //mk->load_geometry(file_name.c_str());

  iGeom* geom = mk->igeom_instance();
  iGeom::Error ierr = geom->load(file_name.c_str());
  IBERRCHK(ierr, "Trouble loading the geometry file.");
  std::vector<iBase_EntityHandle> entities_out, ini_surf, ini_curve, ini_verts;
  geom->getEntities(0, iBase_REGION, entities_out);
  geom->getEntities(0, iBase_FACE, ini_surf);
  geom->getEntities(0, iBase_EDGE, ini_curve);
  geom->getEntities(0, iBase_VERTEX, ini_verts);

  std::cout << " Reg:" << entities_out.size() << " f:" << ini_surf.size() << " c: " << ini_curve.size()
		  << " v: " << ini_verts.size() << std::endl;
  ierr = geom->imprintEnts(&entities_out[0], (int)entities_out.size());
  IBERRCHK(ierr, "Trouble imprinting ");
  double dTol = 1.0e-2;
  ierr = geom->mergeEnts(&entities_out[0],entities_out.size(), dTol);
  IBERRCHK(ierr, "Trouble merging.");
  std::vector<iBase_EntityHandle> vols, surfs, curves, verts;
  geom->getEntities(0, iBase_REGION, vols);
  geom->getEntities(0, iBase_FACE, surfs);
  geom->getEntities(0, iBase_EDGE, curves);
  geom->getEntities(0, iBase_VERTEX, verts);

  std::cout << " after imprint and merge: volumes:" << vols.size() << " faces:" << surfs.size() << " curves: " << curves.size()
        << " vertices:" << verts.size() << std::endl;

  // populate the model with model entities
  mk->populate_model_ents();
  std::cout <<" model populated successfully with entities\n";
  iGeom::TagHandle nameTag;
  ierr = geom->getTagHandle("NAME", nameTag);
  if (0 == ierr)
    for (int dim = 3; dim >= 0; dim--) {
      MEntVector ents;
      mk->get_entities_by_dimension(dim, ents);
      for (unsigned int si = 0u; si < ents.size(); ++si) {
        ModelEnt * me = ents[si];
        iGeom::EntityHandle geo_h = me->geom_handle();
        char name_ent[100];
        geom->getData(geo_h, nameTag, name_ent);
        std::cout <<" dim: " << dim << " index: " << si << "  name:" << name_ent << "\n";
      }
    }
  // get the model entities of interest, for the meshop dim selected
  // int dim_primary_dim=2; // could change for tet mesher, for example
  // this tool will work only for surf meshers (now)
  MEntVector primary_ents;
  mk->get_entities_by_dimension(mesh_op_dim, primary_ents);

  std::cout<<" mesh using " << meshop << " " << primary_ents.size() << " model entities\n";

  // Construct the AF2DfltTriangleMeshOp on the surfaces
  MeshOp * meshOp = mk->construct_meshop(meshop.c_str(), primary_ents);

  SizingFunction* sfPtr = new SizingFunction(mk, -1, size);
  for (unsigned int si = 0u; si < primary_ents.size(); ++si)
  {
    primary_ents[si]->sizing_function_index(sfPtr->core_index());
  }
  // also set some edges ?

  meshOp->set_debug_verbosity(debug);

  mk->setup_and_execute();

  delete sfPtr;

  // save output
  mk->moab_instance()->write_file(out_mesh.c_str());

  delete mk;
  return 0;
}




