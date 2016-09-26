#include "meshkit/AF2PntTrnsfrmLnrV.hpp"

// C++
#include <cstddef>
#include <sstream>

// MeshKit
#include "meshkit/Error.hpp"

AF2PntTrnsfrmLnrV::AF2PntTrnsfrmLnrV(
    std::list<const AF2RuleExistVertex*> vertices,
    std::list<double> xDiffCoeff, std::list<double> yDiffCoeff) :
    refVertices(vertices.begin(), vertices.end()),
    xCoeff(xDiffCoeff.begin(), xDiffCoeff.end()),
    yCoeff(yDiffCoeff.begin(), yDiffCoeff.end())
{
  if (xDiffCoeff.size() != 2 * vertices.size())
  {
    std::ostringstream errStringStream;
    errStringStream << "The list of x coefficients (size "
        << xDiffCoeff.size() << ") is not twice as long as the list of "
        << "vertices (size " << vertices.size() << ").";
    MeshKit::Error badArg(MeshKit::ErrorCode::MK_BAD_INPUT);
    badArg.set_string(errStringStream.str().c_str());
    throw badArg;
  }
  if (yDiffCoeff.size() != 2 * vertices.size())
  {
    std::ostringstream errStringStream;
    errStringStream << "The list of y coefficients (size "
        << yDiffCoeff.size() << ") is not twice as long as the list of "
        << "vertices (size " << vertices.size() << ").";
    MeshKit::Error badArg(MeshKit::ErrorCode::MK_BAD_INPUT);
    badArg.set_string(errStringStream.str().c_str());
    throw badArg;
  }
}

AF2PntTrnsfrmLnrV::AF2PntTrnsfrmLnrV(
    std::vector<const AF2RuleExistVertex*> vertices,
    std::vector<double> xDiffCoeff, std::vector<double> yDiffCoeff) :
    refVertices(vertices), xCoeff(xDiffCoeff), yCoeff(yDiffCoeff)
{
  if (xDiffCoeff.size() != 2 * vertices.size())
  {
    std::ostringstream errStringStream;
    errStringStream << "The vector of x coefficients (size "
        << xDiffCoeff.size() << ") is not twice as long as the vector of "
        << "vertices (size " << vertices.size() << ").";
    MeshKit::Error badArg(MeshKit::ErrorCode::MK_BAD_INPUT);
    badArg.set_string(errStringStream.str().c_str());
    throw badArg;
  }
  if (yDiffCoeff.size() != 2 * vertices.size())
  {
    std::ostringstream errStringStream;
    errStringStream << "The vector of y coefficients (size "
        << yDiffCoeff.size() << ") is not twice as long as the vector of "
        << "vertices (size " << vertices.size() << ").";
    MeshKit::Error badArg(MeshKit::ErrorCode::MK_BAD_INPUT);
    badArg.set_string(errStringStream.str().c_str());
    throw badArg;
  }
}

AF2PntTrnsfrmLnrV* AF2PntTrnsfrmLnrV::clone() const
{
  return new AF2PntTrnsfrmLnrV(refVertices, xCoeff, yCoeff);
}

AF2Point2D AF2PntTrnsfrmLnrV::transformPoint(AF2Point2D const & point,
    AF2Binding const & vBinding) const
{
  double xOffset = 0.0;
  double yOffset = 0.0;
  for (unsigned int i = 0; i < refVertices.size(); ++i)
  {
    const AF2RuleExistVertex* refVertex = refVertices[i];
    const AF2Point2D* boundVal = vBinding.getBoundValue(refVertex);
    if (boundVal == NULL)
    {
      MeshKit::Error badArg(MeshKit::ErrorCode::MK_BAD_INPUT);
      badArg.set_string("The binding does not have a bound value for all reference vertices.");
      throw badArg;
    }
    double xDiff = boundVal->getX() - refVertex->getX();
    double yDiff = boundVal->getY() - refVertex->getY();
    xOffset += xCoeff[2*i]*xDiff + xCoeff[2*i + 1]*yDiff;
    yOffset += yCoeff[2*i]*xDiff + yCoeff[2*i + 1]*yDiff;
  }

  AF2Point2D translated(point.getX() + xOffset, point.getY() + yOffset);
  return translated;
}
