/**
 * \file triMeshOp.cpp \test
 *
 * \brief Test the AF2DfltTriangleMeshOp on a variety of surfaces, both flat
 *   and curved.
 *
 * In make check the meshes are not saved, but if this test is run separately,
 * it accepts a command line argument specifying the file extension for the
 * type of meshes that it should save.
 */
// C++
#include <cstddef>
#include <iostream>
#include <string>

// MeshKit
#include "meshkit/MKCore.hpp"
#include "meshkit/AF2DfltTriangleMeshOp.hpp"
#include "meshkit/SizingFunction.hpp"
#include "meshkit/SizingFunctionVar.hpp"
#include "meshkit/ModelEnt.hpp"

// MeshKit testing
#include "TestUtil.hpp"

// define the geometry file extension depending on the geometry model
#ifdef HAVE_ACIS
std::string geomExt = ".sat";
#elif HAVE_OCC
std::string geomExt = ".stp";
#else
std::string geomExt = ".facet";
#define HAVE_FACET
#endif

using namespace MeshKit;

MEntVector constructMeshOp(std::string geomFile);
void testMesh(std::string geomFile,
    std::string meshFile, SizingFunction* sfPtr);

void testSquare();
void testHoleySurf();
void testSingleHoleSurf();
void testSingleHoleSurfImprinted();
void testSquareVarSize();
void testPieceOfTorus();
void testSphere();
void testBrick();


// This variable is at global scope because calling deleteAll on
// the MKCore geometry instance appears to cause memory inconsistencies
// with later use of the geometry instance
MKCore* mk = NULL;

// These variables are at global scope because they affect all of
// the tests.  There is no way provided to save only some of the
// meshes
bool saveMesh = false;
std::string meshExt = "";

int main(int argc, char **argv) 
{
  // This variable is defined and used in main because a new MKCore
  // instance cannot be easily constructed after another MKCore
  // instance is deleted; there are problems with a tag left behind in
  // iGeom.
  mk = new MeshKit::MKCore();

  int num_fail = 0;

  if (argc == 2) 
  {
    meshExt = argv[1];
    if (meshExt.size() > 0 && meshExt.at(0) != '.')
    {
      meshExt.insert(0, ".");
    }
    saveMesh = true;
    std::cout << "Requested meshes saved with file extension "
        << meshExt << std::endl;
  }
  num_fail += RUN_TEST(testSquare);
  num_fail += RUN_TEST(testSquareVarSize);
#ifdef HAVE_FACET
  num_fail += RUN_TEST(testBrick);
#else
  num_fail += RUN_TEST(testHoleySurf);
  num_fail += RUN_TEST(testSingleHoleSurf);
  num_fail += RUN_TEST(testSingleHoleSurfImprinted);
  num_fail += RUN_TEST(testPieceOfTorus);
  num_fail += RUN_TEST(testSphere);
#endif
  delete mk;

  return num_fail;
}

/**
 * Load a geometry file from the test directory with the specified file name
 * (after appending a  file name suffix depending on the geometry engine).
 * Then extract all of the newly loaded two-dimensional geometry entities
 * (surfaces) and create an AF2DfltTriangleMeshOp with those surfaces.
 *
 * Return an MEntVector containing the new surfaces.
 */
MEntVector constructMeshOp(std::string geomFile)
{
  // count the number of surfaces that existing before
  MEntVector temp;
  mk->get_entities_by_dimension(2, temp);
  unsigned int countExisting = temp.size();
  temp.clear();

  // load the geometry
  std::string file_name = TestDir + "/" + geomFile + geomExt;
  mk->load_geometry(file_name.c_str());

  // extract the new surfaces into an MEntVector
  MEntVector newSurfs;
  mk->get_entities_by_dimension(2, temp);
  unsigned int countAll = temp.size();
  unsigned int count = 0u;
  MEntVector::reverse_iterator surfItr = temp.rbegin();
  while (count < (countAll - countExisting))
  {
    newSurfs.push_back(*surfItr);
    ++count;
    ++surfItr;
  }

  // Construct the AF2DfltTriangleMeshOp on the surfaces
  mk->construct_meshop("AF2DfltTriangleMeshOp", newSurfs);

  return newSurfs;
}

/**
 * Load a geometry file from the test directory with the specified file name
 * (after appending a  file name suffix depending on the geometry engine).
 * Then mesh all of the newly loaded two-dimensional geometry entities
 * (surfaces) using an AF2DfltTriangleMeshOp.
 *
 * After meshing is complete, report the number of triangles, save the
 * mesh to a file with the specified mesh file name (after appending a
 * suffix depending on what the user wanted) if requested, and delete the
 * triangles from the mesh.
 */
void testMesh(std::string geomFile,
    std::string meshFile, SizingFunction* sfPtr)
{
  MEntVector surfs = constructMeshOp(geomFile);

  // Set the size on the surfaces and their children edges
  // if a size was provided
  if (sfPtr != NULL)
  {
    for (unsigned int si = 0u; si < surfs.size(); ++si)
    {
      surfs[si]->sizing_function_index(sfPtr->core_index());
    }
  }

  // setup and execute the MeshOp to mesh the surface
  mk->setup_and_execute();

  // report the number of triangles
  moab::Range tris;
  moab::ErrorCode rval =
      mk->moab_instance()->get_entities_by_dimension(0, 2, tris);
  CHECK_EQUAL(moab::MB_SUCCESS, rval);
  std::cout << tris.size() << " tris generated." << std::endl;

  // save the mesh if directed to do so
  if (saveMesh) {
    std::string outfile = meshFile + meshExt;
    moab::EntityHandle* outputSets = new moab::EntityHandle[surfs.size()];
    for (unsigned int si = 0u; si < surfs.size(); ++si)
    {
      outputSets[si] = surfs[si]->mesh_handle();
    }
    rval = mk->moab_instance()->write_file(outfile.c_str(), NULL, NULL,
        outputSets, surfs.size());
    MBERRCHK(rval, mk->moab_instance());
    delete[] outputSets;
  }

  // remove the MeshOp and the rest of the graph
  mk->clear_graph();

  // delete the triangles from the mesh
  mk->moab_instance()->delete_entities(tris);
}

void testSquare()
{
  // Provide a sizing function that will become the index 0 sizing function
  // and will be used on the edges.  Without this, there will be a failure
  // in the edge mesher.
  SizingFunction* sfPtr = new SizingFunction(mk, 5, -1);
  // Do not pass a sizing function to the surface, though.  It is not required.
  testMesh("squaresurf", "squaresurf", NULL);
  delete sfPtr;
}

void testHoleySurf() 
{
  SizingFunction esize(mk, -1, 0.25);
  testMesh("holysurf", "holysurf", &esize);
}

void testSingleHoleSurf() 
{
  SizingFunction esize(mk, -1, 0.25);
  testMesh("singleholesurf", "singleholesurf", &esize);
}

void testSingleHoleSurfImprinted()
{
  SizingFunction esize(mk, -1, 0.25);
  testMesh("singleholesurfimprinted", "singleholesurfimprinted", &esize);
}


void testSquareVarSize()
{
  // make a sizing function and set it on the surface
  SizingFunctionVar svar(mk, -1, 0.1);

  // these could be read from a file, or something
  double point0[3] = {0, 0, 0};
  double coeffs[4] = {0.05, 0.05, 0.05, 0.1};
  svar.set_linear_coeff(point0, coeffs);

  // run the test
  testMesh("squaresurf", "squaresurf_var", &svar);
}

void testPieceOfTorus()
{
  SizingFunction esize(mk, -1, 0.25);
  testMesh("pieceOfTorus01", "pieceOfTorus01", &esize);
}

void testSphere()
{
  MEntVector surfs = constructMeshOp("sphere");

  // Set the size to 1.0 on all surfaces of the sphere,
  // but don't pass the sizing down to the child edges to test
  // whether the setup method will do that
  SizingFunction esize(mk, -1, 1.0);
  for (unsigned int si = 0u; si < surfs.size(); ++si)
  {
    surfs[si]->sizing_function_index(esize.core_index(), false);
  }

  // setup and execute the MeshOp to mesh the surface
  mk->setup_and_execute();

  // report the number of triangles
  moab::Range tris;
  moab::ErrorCode rval =
      mk->moab_instance()->get_entities_by_dimension(0, 2, tris);
  CHECK_EQUAL(moab::MB_SUCCESS, rval);
  std::cout << tris.size() << " tris generated." << std::endl;

  // save the mesh if directed to do so
  if (saveMesh) {
    std::string outfile = std::string("sphere") + meshExt;
    moab::EntityHandle* outputSets = new moab::EntityHandle[surfs.size()];
    for (unsigned int si = 0u; si < surfs.size(); ++si)
    {
      outputSets[si] = surfs[si]->mesh_handle();
    }
    rval = mk->moab_instance()->write_file(outfile.c_str(), NULL, NULL,
        outputSets, surfs.size());
    MBERRCHK(rval, mk->moab_instance());
    delete[] outputSets;
  }

  // remove the MeshOp and the rest of the graph
  mk->clear_graph();

  // delete the triangles from the mesh
  mk->moab_instance()->delete_entities(tris);
}

void testBrick()
{
  MEntVector surfs = constructMeshOp("brick");

    // Set the size to 1.0 on all surfaces of the sphere,
    // but don't pass the sizing down to the child edges to test
    // whether the setup method will do that
    SizingFunction esize(mk, -1, 0.2);
    for (unsigned int si = 0u; si < surfs.size(); ++si)
    {
      surfs[si]->sizing_function_index(esize.core_index(), false);
    }

    // setup and execute the MeshOp to mesh the surface
    mk->setup_and_execute();

    // report the number of triangles
    moab::Range tris;
    moab::ErrorCode rval =
        mk->moab_instance()->get_entities_by_dimension(0, 2, tris);
    CHECK_EQUAL(moab::MB_SUCCESS, rval);
    std::cout << tris.size() << " tris generated." << std::endl;

    // save the mesh if directed to do so
    if (saveMesh) {
      std::string outfile = std::string("brick") + meshExt;
      moab::EntityHandle* outputSets = new moab::EntityHandle[surfs.size()];
      for (unsigned int si = 0u; si < surfs.size(); ++si)
      {
        outputSets[si] = surfs[si]->mesh_handle();
      }
      rval = mk->moab_instance()->write_file(outfile.c_str(), NULL, NULL,
          outputSets, surfs.size());
      MBERRCHK(rval, mk->moab_instance());
      delete[] outputSets;
    }

    // remove the MeshOp and the rest of the graph
    mk->clear_graph();

    // delete the triangles from the mesh
    mk->moab_instance()->delete_entities(tris);
}
