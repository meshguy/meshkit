#include "meshkit/AF2Polygon3D.hpp"
#include "meshkit/Error.hpp"

AF2Polygon3D::AF2Polygon3D(
    std::list<const AF2Point3D*> const & polygonVertices) :
    numVertices(polygonVertices.size())
{
  typedef std::list<const AF2Point3D*>::const_iterator ItrType;

  // check for the exceptional case of a null pointer in the list of vertices
  for (ItrType itr = polygonVertices.begin();
      itr != polygonVertices.end(); ++itr)
  {
    if (*itr == NULL)
    {
      MeshKit::Error badArg(MeshKit::ErrorCode::MK_BAD_INPUT);
      badArg.set_string(
          "polygonVertices may not contain a null pointer in AF2Polygon3D.");
      throw badArg;
    }
  }

  // allocate memory and copy the vertex pointers
  vertices = new const AF2Point3D*[numVertices];
  int indx = 0;
  for (ItrType itr = polygonVertices.begin();
      itr != polygonVertices.end(); ++itr)
  {
    vertices[indx] = *itr;
    ++indx;
  }
}

AF2Polygon3D::~AF2Polygon3D()
{
  delete[] vertices;
}

AF2Polygon3D::AF2Polygon3D(const AF2Polygon3D & toCopy) :
    numVertices(toCopy.numVertices)
{
  vertices = new const AF2Point3D*[numVertices];
  for (unsigned int indx = 0; indx < numVertices; ++indx)
  {
    vertices[indx] = toCopy.vertices[indx];
  }
}

AF2Polygon3D& AF2Polygon3D::operator=(const AF2Polygon3D & rhs)
{
  // directly copy the number of vertices
  numVertices = rhs.numVertices;

  // allocate new array in temporary location
  // and copy the pointers into the temporary array
  const AF2Point3D** tempVertices = new const AF2Point3D*[numVertices];
  for (unsigned int indx = 0; indx < numVertices; ++indx)
  {
    tempVertices[indx] = rhs.vertices[indx];
  }

  // delete what used to be held by the array of vertices
  // (may be null if default constructor was used)
  delete[] vertices;

  // transfer ownership from the temporary array to a member of this object
  vertices = tempVertices;
  tempVertices = NULL; // not necessary, but to be explicit

  return *this;
}

unsigned int AF2Polygon3D::getNumVertices() const
{
  return numVertices;
}

const AF2Point3D* AF2Polygon3D::getVertex(unsigned int vtxNum) const
{
  if (vtxNum >= numVertices)
  {
    MeshKit::Error badArg(MeshKit::ErrorCode::MK_BAD_INPUT);
    badArg.set_string("The vertex index is too large.");
    throw badArg;
  }
  return vertices[vtxNum];
}
