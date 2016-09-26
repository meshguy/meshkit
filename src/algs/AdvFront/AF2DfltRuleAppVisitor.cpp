#include "meshkit/AF2DfltRuleAppVisitor.hpp"

// C++
#include <cmath>
#include <cstddef>

// MeshKit
#include "meshkit/Error.hpp"

const double AF2DfltRuleAppVisitor::eqTriAreaPerimSqRatio = sqrt(3.0) / 36.0;

AF2DfltRuleAppVisitor::AF2DfltRuleAppVisitor() :
    bestMetricVal(0), bestRuleApp(NULL)
{
}

AF2DfltRuleAppVisitor::~AF2DfltRuleAppVisitor()
{
  if (bestRuleApp != NULL)
  {
    delete bestRuleApp;
  }
}

AF2DfltRuleAppVisitor::AF2DfltRuleAppVisitor(
    const AF2DfltRuleAppVisitor & toCopy) :
    bestMetricVal(toCopy.bestMetricVal)
{
  if (toCopy.bestRuleApp != NULL)
  {
    bestRuleApp = new AF2RuleApplication(*(toCopy.bestRuleApp));
  }
}

AF2DfltRuleAppVisitor& AF2DfltRuleAppVisitor::operator=(
    const AF2DfltRuleAppVisitor & rhs)
{
  bestMetricVal = rhs.bestMetricVal;

  if (bestRuleApp == NULL)
  {
    if (rhs.bestRuleApp != NULL)
    {
      bestRuleApp = new AF2RuleApplication(*(rhs.bestRuleApp));
    }
  }
  else
  {
    AF2RuleApplication* tempBestRuleApp = NULL;
    if (rhs.bestRuleApp != NULL)
    {
      tempBestRuleApp = new AF2RuleApplication(*(rhs.bestRuleApp));
    }
    delete bestRuleApp;
    bestRuleApp = tempBestRuleApp;
  }

  return *this;
}

const AF2RuleApplication* AF2DfltRuleAppVisitor::getBestRuleApplication() const
{
  return bestRuleApp;
}

void AF2DfltRuleAppVisitor::visit(AF2RuleApplication const & ruleApp)
{
  double metricVal = 0.0;

  // loop over the faces and take the worst face metric value
  // as the rule metric value
  for (unsigned int faceIndex = 0u;
      faceIndex < ruleApp.getNumNewFaces(); ++faceIndex)
  {
    const AF2Polygon2D* facePtr = ruleApp.getNewFace(faceIndex);

    // verify that the face is a triangle
    if (facePtr->getNumVertices() != 3u)
    {
      MeshKit::Error badArg(MeshKit::ErrorCode::MK_BAD_INPUT);
      badArg.set_string(
          "AF2DfltRuleAppVisitor can visit only triangular faces.");
      throw badArg;
    }

    // access the vertices of the triangle
    const AF2Point2D* faceVtx0 = facePtr->getVertex(0u);
    const AF2Point2D* faceVtx1 = facePtr->getVertex(1u);
    const AF2Point2D* faceVtx2 = facePtr->getVertex(2u);

    // calculate the components of the vectors along triangle edges
    double faceVec01X = faceVtx1->getX() - faceVtx0->getX();
    double faceVec01Y = faceVtx1->getY() - faceVtx0->getY();
    double faceVec02X = faceVtx2->getX() - faceVtx0->getX();
    double faceVec02Y = faceVtx2->getY() - faceVtx0->getY();
    double faceVec12X = faceVtx2->getX() - faceVtx1->getX();
    double faceVec12Y = faceVtx2->getY() - faceVtx1->getY();

    // calculate the length of each edge from the vector
    double sqLen = faceVec01X * faceVec01X + faceVec01Y * faceVec01Y;
    double faceVec01Len = sqrt(sqLen);
    sqLen = faceVec02X * faceVec02X + faceVec02Y * faceVec02Y;
    double faceVec02Len = sqrt(sqLen);
    sqLen = faceVec12X * faceVec12X + faceVec12Y * faceVec12Y;
    double faceVec12Len = sqrt(sqLen);

    // caculate the perimeter and area
    double facePerim = faceVec01Len + faceVec02Len + faceVec12Len;
    double faceArea =
        0.5 * (faceVec01X * faceVec02Y - faceVec01Y * faceVec02X);
    if (faceArea <= 0.0)
    {
      // the face is inverted, so this rule application is unacceptable
      return;
    }

    // calculate the face metric value
    double faceMetricVal = 10.0 * (eqTriAreaPerimSqRatio *
        facePerim * facePerim / faceArea - 1)  + // shape error
        1.0/faceVec01Len + faceVec01Len +
        1.0/faceVec02Len + faceVec02Len +
        1.0/faceVec12Len + faceVec12Len - 6.0; // size error

    // update the rule metric value, if necessary
    if (faceMetricVal > metricVal)
    {
      metricVal = faceMetricVal;
    }
  }

  // update the best rule application if this is the only rule application
  // or is a "measurable" improvement over the previous best
  if (bestRuleApp == NULL)
  {
    bestRuleApp = new AF2RuleApplication(ruleApp);
    bestMetricVal = metricVal;
  }
  else if (metricVal < 0.99*bestMetricVal)
  {
    *bestRuleApp = ruleApp;
    bestMetricVal = metricVal;
  }
}
