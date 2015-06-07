#include "meshkit/NGTriMesher.hpp"

namespace nglib {
#include "nglib.h"
}

#include "meshkit/MeshScheme.hpp"
#include "meshkit/ModelEnt.hpp"
#include "meshkit/iGeom.hpp"
#include "moab/Interface.hpp"
#include "moab/ReadUtilIface.hpp"
#include "moab/Range.hpp"
#include "meshkit/RegisterMeshOp.hpp"
#include "RefEntity.hpp"

#ifdef PARALLEL
#ifdef HAVE_PARALLEL_MOAB
#include "moab/ParallelComm.hpp"
#endif
#endif

#include <cstddef>
#include <cstdio>
#include <vector>
#include <iostream>

namespace MeshKit
{

moab::EntityType NGTriMesher::meshTps[] = {moab::MBVERTEX, moab::MBTRI, moab::MBMAXTYPE};

NGTriMesher::NGTriMesher(MKCore *mk_core, const MEntVector &me_vec)
        : MeshScheme(mk_core, me_vec)
{
}

NGTriMesher::~NGTriMesher()
{
}

void NGTriMesher::setup_this()
{
//  std::cout << "Entering NGTriMesher.setup_this()" << std::endl;
  // just call setup_boundary, since that's the only constraint we have
  // setup_boundary();
//  std::cout << "Exiting NGTriMesher.setup_this()" << std::endl;
}

void NGTriMesher::execute_this()
{
/*
  iGeom_save(iGeom_Instance instance, const char *name, const char* options,
    int* err, int name_len, int options_len)

  // only valid option in CGMA seems to be
  // TYPE=<fileType>, where fileType can be IGES, STEP, ACIS_SAT, or OCC
  // options may need to begin and end with a space character as delimiter
  // The file type can also be determined from the extension of the file name
  // IGES as .igs, .iges
  // STEP as .stp, .step
  // ACIS_SAT as .sat
  // OCC as .occ, .brep
  // Eventually this calls CubitCompat_export_solid_model
  // . . . I don't see how the iGeom_Instance is used

*/
//  std::cout << "Entering NGTriMesher.execute_this()" << std::endl;


  for (MEntSelection::iterator sit = mentSelection.begin(); sit != mentSelection.end(); sit++) {
    // make a me, for convenience
    ModelEnt* me = (*sit).first;

    std::string fileName = "temp.brep";
    me->igeom_instance()->save(fileName.c_str(), "");

    nglib::Ng_Init();
    nglib::Ng_OCC_Geometry* ngOccGeom =
      nglib::Ng_OCC_Load_BREP(fileName.c_str());

    // delete the temporary file
    remove(fileName.c_str());

    nglib::Ng_Meshing_Parameters* meshParams =
      new nglib::Ng_Meshing_Parameters();

    double my_size = me->mesh_interval_size();
    if (0 > my_size)
    {
      // no size has been set
      meshParams->uselocalh = 1;
      meshParams->maxh = 1.0;
    }
    else
    {
      meshParams->uselocalh = 0;
      meshParams->maxh = my_size;
      // netgen optimization fails if uselocalh = 0
      meshParams->optsurfmeshenable = 0;
    }

    nglib::Ng_Result result;
    nglib::Ng_Mesh* netGenMesh = nglib::Ng_NewMesh();
    result = nglib::Ng_OCC_SetLocalMeshSize(ngOccGeom, netGenMesh, meshParams);
    if (result != nglib::NG_OK)
    {
      ECERRCHK(MK_FAILURE,
        "Failed to set local mesh size in netgen triangle mesher.");
    }
    result = nglib::Ng_OCC_GenerateEdgeMesh(ngOccGeom, netGenMesh, meshParams);
    if (result != nglib::NG_OK)
    {
      ECERRCHK(MK_FAILURE,
        "Failed to generate edge mesh in netgen triangle mesher.");
    }
    result =
      nglib::Ng_OCC_GenerateSurfaceMesh(ngOccGeom, netGenMesh, meshParams);
    if (result != nglib::NG_OK)
    {
      ECERRCHK(MK_FAILURE,
        "Failed to generate surface mesh in netgen triangle mesher.");
    }

    int numPoints = nglib::Ng_GetNP(netGenMesh);

    // temporary -- for testing
//    std::cout << "Result: " << result << std::endl;
//    std::cout << "Number of Points: " << numPoints << std::endl;
    // end testing

    // Copy the vertices resulting from the meshing operation to the
    // mesh in MKCore
    std::vector<double> coords;
    double ngPnt[3];
    for (int i = 0; i < numPoints; ++i)
    {
      nglib::Ng_GetPoint(netGenMesh, i + 1, ngPnt);
      for (int j = 0; j < 3; ++j)
      {
        coords.push_back(ngPnt[j]);
      }
    }

    moab::Range& newEntities = (*sit).second;
    moab::ErrorCode rval;
    rval = mk_core()->moab_instance()->create_vertices(&coords[0],
        numPoints, newEntities);
    MBERRCHK(rval, mk_core()->moab_instance());

    int numTriangles = nglib::Ng_GetNSE(netGenMesh);
    // temporary -- for testing
//    std::cout << "Number of Surface Elements: " << numTriangles << std::endl;
    // end testing

    // Transfer the results of the meshing operation to the
    // mesh in MKCore
    moab::ReadUtilIface *iface;
    rval = mk_core()->moab_instance()->query_interface(iface);
    MBERRCHK(rval, mk_core()->moab_instance());

    moab::EntityHandle starth, *connect;
    rval = iface->get_element_connect(numTriangles, 3, moab::MBTRI,
        1, starth, connect);
    MBERRCHK(rval, mk_core()->moab_instance());

    int ngElm[3];
    for (int i = 0; i < numTriangles; ++i)
    {
      nglib::Ng_GetSurfaceElement(netGenMesh, i + 1, ngElm);
      for (int j = 0; j < 3; ++j)
      {
        connect[3*i + j] = newEntities[ngElm[j]];
      }
    }

    // put new tris into new entity range, then commit the mesh
    newEntities.insert(starth, starth + numTriangles - 1);
    me->commit_mesh(newEntities, COMPLETE_MESH);

    nglib::Ng_DeleteMesh(netGenMesh);
    nglib::Ng_OCC_DeleteGeometry(ngOccGeom);
    delete meshParams;

    nglib::Ng_Exit();
  }

//  std::cout << "Exiting NGTriMesher.execute_this()" << std::endl;
}

} // namespace MeshKit

