#include "meshkit/AF2DfltTriangleMeshOp.hpp"

// C++
#include <cstddef>
#include <list>
#include <vector>

// MOAB
#include "moab/Interface.hpp"
#include "moab/Range.hpp"
#include "moab/ReadUtilIface.hpp"
#include "moab/Types.hpp"

// MeshKit
#include "meshkit/AF2Algorithm.hpp"
#include "meshkit/AF2AlgorithmResult.hpp"
#include "meshkit/AF2DfltPlaneProjMaker.hpp"
#include "meshkit/AF2DfltTriangleRules.hpp"
#include "meshkit/AF2Point3D.hpp"
#include "meshkit/Error.hpp"
#include "meshkit/MeshScheme.hpp"
#include "meshkit/ModelEnt.hpp"
#include "meshkit/RegisterMeshOp.hpp"

namespace MeshKit
{

moab::EntityType AF2DfltTriangleMeshOp::meshTypes[] =
    {moab::MBVERTEX, moab::MBTRI, moab::MBMAXTYPE};

bool AF2DfltTriangleMeshOp::can_mesh(iBase_EntityType dimension)
{
  return dimension == iBase_FACE;
}

bool AF2DfltTriangleMeshOp::can_mesh(ModelEnt *me)
{
  return canmesh_face(me);
}

const char* AF2DfltTriangleMeshOp::name()
{
  return "AF2DfltTriangleMeshOp";
}

const moab::EntityType* AF2DfltTriangleMeshOp::output_types()
{
  return meshTypes;
}

AF2DfltTriangleMeshOp::AF2DfltTriangleMeshOp(
    MKCore *meshkitCore, const MEntVector &meshEntVec) :
    MeshScheme(meshkitCore, meshEntVec)
{
  // don't do anything beyond constructing superclass with correct arguments
}

AF2DfltTriangleMeshOp::~AF2DfltTriangleMeshOp()
{
  // don't do anything; superclass destructor will be called
}

AF2DfltTriangleMeshOp::AF2DfltTriangleMeshOp(
    const AF2DfltTriangleMeshOp& toCopy) : MeshScheme(toCopy)
{
  // don't do anything beyond constructing superclass with correct arguments
}

AF2DfltTriangleMeshOp& AF2DfltTriangleMeshOp::operator=(
    const AF2DfltTriangleMeshOp& rhs)
{
  Error notImpl(ErrorCode::MK_NOT_IMPLEMENTED);
  notImpl.set_string(
      "AF2DfltTriangleMeshOp assignment operator is not supported.");
  throw notImpl;
}

void AF2DfltTriangleMeshOp::execute_this()
{
  AF2DfltTriangleRules defaultRules;
  AF2Algorithm algorithm(defaultRules.getRules());

  int debug = get_debug_verbosity(); // will call MeshOp to get the debug verbosity level

  // if some debug level, dump out the mesh so far
  if (debug>0)
  {
    // get the moab instance for convenience
    moab::Interface* moabInstance = mk_core()->moab_instance();
    moabInstance->write_file("meshBeforeAF2.vtk");
  }
  for (MEntSelection::iterator selIt = mentSelection.begin();
      selIt != mentSelection.end(); ++selIt)
  {
    // extract the model entity from the map iterator
    ModelEnt* modelEnt = selIt->first;

    if (modelEnt->get_meshed_state() >= COMPLETE_MESH)
    {
      continue;
    }

    // get the moab instance for convenience
    moab::Interface* moabInstance = mk_core()->moab_instance();

    // query the boundary to find mesh edges
    // Note: the senses returned by this method should be either FORWARD
    //   or REVERSE, not BOTH . . .
    //   but I think there is a chance that mesh edges that are on a
    //   geometric edge that has a sense of BOTH might not agree with
    //   nearby edges and thus might not form a wire when placed end
    //   to end.  The AF2Algorithm doesn't require that the input be in
    //   wire order, and should be okay as long as the edge actually
    //   does appear with both senses.
    std::vector<moab::EntityHandle> bndryEdges;
    std::vector<int> bEdgeSenses;
    modelEnt->boundary(1, bndryEdges, &bEdgeSenses);

    // get lists of the vertex coordinates and handles
    // along with an index from the edge handles into those lists
    std::vector<int> bndryIds;
    std::vector<double> coords;
    moab::Range bndryVertsRange;
    modelEnt->get_indexed_connect_coords(
        bndryEdges, &bEdgeSenses, NULL, bndryIds, coords, &bndryVertsRange);
    std::vector<moab::EntityHandle> bndryVertsVec(
        bndryVertsRange.begin(), bndryVertsRange.end());

    // convert the edges into the unsigned int format that is needed for
    // input into the AF2Algorithm
    // Note: The tags used inside get_indexed_connect_coords are actually
    // unsigned int (as of May 2016)
    std::vector<unsigned int> inputEdges;
    for (std::vector<int>::const_iterator bidItr = bndryIds.begin();
        bidItr != bndryIds.end(); ++bidItr)
    {
      inputEdges.push_back(*bidItr);
    }

    // determine the number of boundary edges and vertices
    unsigned int numVertices = bndryVertsVec.size();
    unsigned int numEdges = inputEdges.size() / 2;

    // check that there is reasonable data to initialize the front
    if (numEdges <= 0)
    {
      Error failErr(ErrorCode::MK_FAILURE);
      failErr.set_string("There are no boundary mesh edges to use to initialize the advancing front in AF2DfltTriangleMeshOp.");
      throw failErr;
    }
    // This should never happen if there are edges, but as a sanity check
    // since the coordinates of the first vertex may be used to check sizing
    if (numVertices <= 0)
    {
      Error failErr(ErrorCode::MK_FAILURE);
      failErr.set_string(
          "No boundary vertices found in AF2DfltTriangleMeshOp.");
      throw failErr;
    }

    // configure the sizing function, if any, that will be used
    SizingFunction* sizing = NULL;
    if (modelEnt->sizing_function_index() > -1)
    {
      int sfIndex = modelEnt->sizing_function_index();
      SizingFunction* tempSizing = mk_core()->sizing_function(sfIndex);
      // check the desired size at the coordinates of the first vertex
      if (tempSizing->size(&coords[0]) > 0)
      {
        sizing = tempSizing;
      }
    }

    // construct the local transform maker and run the algorithm
    AF2DfltPlaneProjMaker localTransformMaker(
        modelEnt->igeom_instance(), modelEnt->geom_handle(), sizing);

    AF2AlgorithmResult* meshResult = algorithm.execute(&localTransformMaker,
        &coords[0], numVertices, &inputEdges[0], numEdges, &bndryVertsVec[0], debug);

    // throw failure if the algorithm did not succeed
    if (!meshResult->isSuccessful())
    {
      delete meshResult;
      Error failErr(ErrorCode::MK_FAILURE);
      failErr.set_string("AF2DfltTriangleMeshOp failed.");
      throw failErr;
    }

    // Collect the new points in a vector
    const std::list<AF2Point3D*>* pointList = meshResult->getPoints();
    std::vector<double> newPointCoords;
    typedef std::list<AF2Point3D*>::const_iterator ConstPointItr;
    int numNewPoints = 0;
    for (ConstPointItr pItr = pointList->begin();
        pItr != pointList->end(); ++pItr)
    {
      if (!(*pItr)->isCommitted())
      {
        ++numNewPoints;
        newPointCoords.push_back((*pItr)->getX());
        newPointCoords.push_back((*pItr)->getY());
        newPointCoords.push_back((*pItr)->getZ());
      }
    }
    // Commit the new points to MOAB
    moab::Range newPointsRange;
    moabInstance->create_vertices(
        &newPointCoords[0], numNewPoints, newPointsRange);
    // Set the MOAB handles for the now-committed vertices
    ConstPointItr np3dItr = pointList->begin();
    for (moab::Range::const_iterator nphItr = newPointsRange.begin();
        nphItr != newPointsRange.end(); ++nphItr)
    {
      while ((*np3dItr)->isCommitted())
      {
        ++np3dItr;
      }
      (*np3dItr)->setCommittedHandle(*nphItr);
      ++np3dItr;
    }

    // get a pointer to the list of new triangles and check how many there are
    const std::list<const AF2Polygon3D*>* faceList = meshResult->getFaces();
    int numTriangles = faceList->size();

    // pre-allocate connectivity memory to store the triangles
    moab::ReadUtilIface* readInterface;
    moab::ErrorCode moabRet;
    moabRet = moabInstance->query_interface(readInterface);
    MBERRCHK(moabRet, moabInstance);
    moab::EntityHandle firstHandle;
    moab::EntityHandle* triConnectAry;
    moabRet = readInterface->get_element_connect(
        numTriangles, 3, moab::MBTRI, 0, firstHandle, triConnectAry);
    MBERRCHK(moabRet, moabInstance);

    // fill the connectivity array
    unsigned int caIndex = 0u;
    typedef std::list<const AF2Polygon3D*>::const_iterator ConstTriangleItr;
    for (ConstTriangleItr tItr = faceList->begin();
        tItr != faceList->end(); ++tItr)
    {
      for (unsigned int vi = 0; vi < 3; ++vi)
      {
        triConnectAry[caIndex + vi] = (*tItr)->getVertex(vi)->getVertexHandle();
      }
      caIndex += 3;
    }

    delete meshResult;

    // Insert the new entities into the entity range
    selIt->second.insert(firstHandle, firstHandle + numTriangles - 1);

    // Commit the mesh
    modelEnt->commit_mesh(selIt->second, COMPLETE_MESH);
  }
}

const moab::EntityType* AF2DfltTriangleMeshOp::mesh_types_arr() const
{
  return output_types();
}

void AF2DfltTriangleMeshOp::setup_this()
{
  // Check that the dimension is correct for each model entity and pass
  // down the sizing functions that are set on the selected surfaces
  // to their boundary edges if those edges don't have sizing functions set
  for (MEntSelection::iterator selIt = mentSelection.begin();
      selIt != mentSelection.end(); ++selIt)
  {
    // extract the model entity from the map iterator
    ModelEnt* modelEnt = selIt->first;

    // check that the dimension of the selected model entity is two
    if (modelEnt->dimension() != 2)
    {
      Error dimErr(ErrorCode::MK_WRONG_DIMENSION);
      dimErr.set_string("Found a selected entity of dimension not equal to 2 in AF2DfltTriangleMeshOp.");
      throw dimErr;
    }

    // Try to check whether there is a valid sizing function set on the
    // model entity.
    int sFIndex = -1;
    if (modelEnt->sizing_function_index() > -1)
    {
      int tempSFIndex = modelEnt->sizing_function_index();
      SizingFunction* tempSizing = mk_core()->sizing_function(tempSFIndex);
      // check the desired size at the coordinates of the first vertex
      if (tempSizing->size() > 0)
      {
        sFIndex = tempSFIndex;
      }
    }

    // If the sizing function appears to be valid, pass it down to adjacent
    // children that do not have a sizing function set, because it is best
    // to have the size on the boundaries about right.
    if (sFIndex > -1)
    {
      MEntVector bndryEdges;
      modelEnt->get_adjacencies(1, bndryEdges);
      for (MEntVector::iterator beItr = bndryEdges.begin();
          beItr != bndryEdges.end(); ++beItr)
      {
        int beSFIndex = (*beItr)->sizing_function_index();
        if (beSFIndex == -1)
        {
          (*beItr)->sizing_function_index(sFIndex);
        }
      }
    }
  }

  // Set up the default MeshOp on any children that need one
  setup_boundary();

  // Make sure that this MeshOp depends on all MeshOps set on its facets
  ensure_facet_dependencies(false);
}

} // namespace MeshKit
