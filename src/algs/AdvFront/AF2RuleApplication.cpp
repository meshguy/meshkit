#include "meshkit/AF2RuleApplication.hpp"
#include "meshkit/Error.hpp"

AF2RuleApplication::AF2RuleApplication(
    std::list<const AF2Point2D*> const & newPointsList,
    std::list<const AF2Polygon2D*> const & newFacesList) :
    numNewPoints(newPointsList.size()), numNewFaces(newFacesList.size())
{
  // type definitions for two iterators that the constructor will use
  typedef std::list<const AF2Point2D*>::const_iterator ConstPntItrType;
  typedef std::list<const AF2Polygon2D*>::const_iterator ConstPolyItrType;

  // check whether there is a null-valued pointer in the list of new points
  for (ConstPntItrType itr = newPointsList.begin();
      itr != newPointsList.end(); ++itr)
  {
    if (*itr == NULL)
    {
      MeshKit::Error badArg(MeshKit::ErrorCode::MK_BAD_INPUT);
      badArg.set_string("newPointsList may not contain a null pointer in AF2RuleApplication.");
      throw badArg;
    }
  }

  // check whether there is a null-valued pointer in the list of faces
  for (ConstPolyItrType itr = newFacesList.begin();
      itr != newFacesList.end(); ++itr)
  {
    if (*itr == NULL)
    {
      MeshKit::Error badArg(MeshKit::ErrorCode::MK_BAD_INPUT);
      badArg.set_string(
          "newFacesList may not contain a null pointer in AF2RuleApplication.");
      throw badArg;
    }
  }

  // map from new points passed in to the copies of the new points
  std::map<const AF2Point2D*, const AF2Point2D*> pointPtrCopyMap;

  // now that the constructor has checked things that will cause it
  // to throw its own exceptions, proceed with allocating memory, etc.
  newPoints = new const AF2Point2D*[numNewPoints];
  int indx = 0;
  for (ConstPntItrType itr = newPointsList.begin();
      itr != newPointsList.end(); ++itr)
  {
    newPoints[indx] = new AF2Point2D(*(*itr));
    pointPtrCopyMap[*itr] = newPoints[indx];
    ++indx;
  }

  newFaces = new const AF2Polygon2D*[numNewFaces];
  indx = 0;
  for (ConstPolyItrType itr = newFacesList.begin();
      itr != newFacesList.end(); ++itr)
  {
    mapCopyPolygon(*(*itr), newFaces[indx], pointPtrCopyMap);
    ++indx;
  }
}

AF2RuleApplication::~AF2RuleApplication()
{
  for (unsigned int plygnIndx = 0; plygnIndx < numNewFaces; ++plygnIndx)
  {
    delete newFaces[plygnIndx];
  }
  delete[] newFaces;

  for (unsigned int pntIndx = 0; pntIndx < numNewPoints; ++pntIndx)
  {
    delete newPoints[pntIndx];
  }
  delete[] newPoints;
}

AF2RuleApplication::AF2RuleApplication(const AF2RuleApplication & toCopy) :
    numNewPoints(toCopy.numNewPoints), numNewFaces(toCopy.numNewFaces)
{
  std::map<const AF2Point2D*, const AF2Point2D*> pointPtrCopyMap;

  newPoints = new const AF2Point2D*[numNewPoints];
  for (unsigned int pntIndx = 0; pntIndx < numNewPoints; ++pntIndx)
  {
    newPoints[pntIndx] = new AF2Point2D(*(toCopy.newPoints[pntIndx]));
    pointPtrCopyMap[toCopy.newPoints[pntIndx]] = newPoints[pntIndx];
  }

  newFaces = new const AF2Polygon2D*[numNewFaces];
  for (unsigned int plygnIndx = 0; plygnIndx < numNewFaces; ++plygnIndx)
  {
    mapCopyPolygon(*(toCopy.newFaces[plygnIndx]),
        newFaces[plygnIndx], pointPtrCopyMap);
  }
}

AF2RuleApplication& AF2RuleApplication::operator=(
    const AF2RuleApplication & rhs)
{
  // copy constructor functionality, but copy to
  // other parts of memory rather than to this
  std::map<const AF2Point2D*, const AF2Point2D*> pointPtrCopyMap;

  const AF2Point2D** otherNewPoints = new const AF2Point2D*[rhs.numNewPoints];
  for (unsigned int pntIndx = 0; pntIndx < rhs.numNewPoints; ++pntIndx)
  {
    otherNewPoints[pntIndx] = new AF2Point2D(*(rhs.newPoints[pntIndx]));
    pointPtrCopyMap[rhs.newPoints[pntIndx]] = otherNewPoints[pntIndx];
  }

  const AF2Polygon2D** otherNewFaces = new const AF2Polygon2D*[rhs.numNewFaces];
  for (unsigned int plygnIndx = 0; plygnIndx < rhs.numNewFaces; ++plygnIndx)
  {
    mapCopyPolygon(*(rhs.newFaces[plygnIndx]),
        otherNewFaces[plygnIndx], pointPtrCopyMap);
  }

  // destructor functionality to clean up current members of this
  for (unsigned int plygnIndx = 0; plygnIndx < numNewFaces; ++plygnIndx)
  {
    delete newFaces[plygnIndx];
  }
  delete[] newFaces;

  for (unsigned int pntIndx = 0; pntIndx < numNewPoints; ++pntIndx)
  {
    delete newPoints[pntIndx];
  }
  delete[] newPoints;

  // transfer ownership from other parts of memory to this object
  numNewPoints = rhs.numNewPoints;
  newPoints = otherNewPoints;
  otherNewPoints = NULL; // not necessary, but to be explicit
  numNewFaces = rhs.numNewFaces;
  newFaces = otherNewFaces;
  otherNewFaces = NULL; // not necessary, but to be explicit

  // return reference to this
  return *this;
}

unsigned int AF2RuleApplication::getNumNewFaces() const
{
  return numNewFaces;
}

const AF2Polygon2D* AF2RuleApplication::getNewFace(
    unsigned int newFaceIndex) const
{
  if (newFaceIndex >= numNewFaces)
  {
    MeshKit::Error badArg(MeshKit::ErrorCode::MK_BAD_INPUT);
    badArg.set_string("newFaceIndex is out of range");
    throw badArg;
  }
  return newFaces[newFaceIndex];
}

unsigned int AF2RuleApplication::getNumNewPoints() const
{
  return numNewPoints;
}

const AF2Point2D* AF2RuleApplication::getNewPoint(
    unsigned int newPointIndex) const
{
  if (newPointIndex >= numNewPoints)
  {
    MeshKit::Error badArg(MeshKit::ErrorCode::MK_BAD_INPUT);
    badArg.set_string("newPointIndex is out of range");
    throw badArg;
  }
  return newPoints[newPointIndex];
}

void AF2RuleApplication::mapCopyPolygon(const AF2Polygon2D & sourcePolygon,
    const AF2Polygon2D* & targetPolygon,
    const std::map<const AF2Point2D*, const AF2Point2D*> & vertexMap)
{
  std::list<const AF2Point2D*> mapCopyVertices;
  for (unsigned int vIndx = 0; vIndx < sourcePolygon.getNumVertices(); ++vIndx)
  {
    const AF2Point2D* sourcePolyVertex = sourcePolygon.getVertex(vIndx);
    typedef std::map<const AF2Point2D*, const AF2Point2D*>::const_iterator
        ConstItrType;
    ConstItrType mappedVertex = vertexMap.find(sourcePolyVertex);
    if (mappedVertex == vertexMap.end())
    {
      // The vertex is not a new vertex from rule application.
      // Instead it is from the neighborhood and does not need to be mapped
      mapCopyVertices.push_back(sourcePolyVertex);
    }
    else
    {
      // The vertex is a vertex from the rule application and
      // this code maps it to the appropriate pointer in the object
      // that is being constructed.
      mapCopyVertices.push_back(mappedVertex->second);
    }
  }

  targetPolygon = new AF2Polygon2D(mapCopyVertices);
}
