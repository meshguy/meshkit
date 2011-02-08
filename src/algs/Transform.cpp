#include "meshkit/Transform.hpp"

#include <cassert>
#include <cmath>

static MeshKit::Vector<3> rodrigues(const MeshKit::Vector<3> &pt,
                                    const MeshKit::Vector<3> &z,
                                    double dtheta) {
  MeshKit::Vector<3> x = vector_product(z, pt);

  double a = cos(dtheta);
  double b = sin(dtheta);
  double c = (pt % z)*(1-a);

  return a*pt + b*x + c*z;
}

static inline double * vec2ptr(std::vector< MeshKit::Vector<3> > &v) {
  return reinterpret_cast<double *>(&v[0]);
}


namespace MeshKit {
namespace Copy {

void Transform::operator ()(iMesh &impl, iMesh::EntityHandle *src, int src_size,
                            iMesh::EntityHandle **dest, int *dest_alloc,
                            int *dest_size) const {
  std::vector< Vector<3> > coords(src_size);
  IBERRCHK(impl.getVtxArrCoords(src, src_size, iBase_INTERLEAVED,
                                vec2ptr(coords)), "FIXME");

  for (int i=0; i<coords.size(); i++)
    transform(coords[i]);

  if (dest && src == *dest) {
    IBERRCHK(impl.setVtxArrCoords(src, src_size, iBase_INTERLEAVED,
                                  vec2ptr(coords)), "FIXME");
  }
  else {
    if (dest)
      *dest = static_cast<iMesh::EntityHandle*>(
        malloc(src_size * sizeof(iMesh::EntityHandle)));

    IBERRCHK(impl.createVtxArr(src_size, iBase_INTERLEAVED, vec2ptr(coords),
                               *dest), "FIXME");
  }
}

Identity::Identity() {}

void Identity::transform(Vector<3> &coords) const {}

Translate::Translate(const Vector<3> &dv) : dv_(dv) {}

void Translate::transform(Vector<3> &coords) const {
  coords += dv_;
}

Rotate::Rotate(const Vector<3> &origin, const Vector<3> &z, double dtheta)
  : origin_(origin), z_(z/length(z)), dtheta_(dtheta) {}

void Rotate::transform(Vector<3> &coords) const {
  coords = rodrigues(coords-origin_, z_, dtheta_) + origin_;
}


} // namespace Copy

namespace Extrude {

Transform::Transform(int steps) : steps_(steps) {}

void Transform::operator ()(int step, iMesh &impl, iMesh::EntityHandle *src,
                            int src_size, iMesh::EntityHandle **dest,
                            int *dest_alloc, int *dest_size) const {
  std::vector< Vector<3> > coords(src_size);
  IBERRCHK(impl.getVtxArrCoords(src, src_size, iBase_INTERLEAVED,
                                vec2ptr(coords)), "FIXME");

  for (int i=0; i<coords.size(); i++)
    transform(step, coords[i]);

  if (dest && !dest) {
    *dest = static_cast<iMesh::EntityHandle*>(
      malloc(src_size * sizeof(iMesh::EntityHandle)));
  }
  IBERRCHK(impl.createVtxArr(src_size, iBase_INTERLEAVED, vec2ptr(coords),
                             *dest), "FIXME");
}

Translate::Translate(const Vector<3> &dv, int steps) : Transform(steps),
                                                       dv_(dv/steps) {}

void Translate::transform(int step, Vector<3> &coords) const {
  coords += dv_*step;
}

Rotate::Rotate(const Vector<3> &origin, const Vector<3> &z, double dtheta,
               int steps) : origin_(origin), z_(z/length(z)), Transform(steps),
                            dtheta_(dtheta/steps) {}

void Rotate::transform(int step, Vector<3> &coords) const {
  coords = rodrigues(coords-origin_, z_, dtheta_*step) + origin_;
}

} // namespace Extrude

} // namespace MeshKit
