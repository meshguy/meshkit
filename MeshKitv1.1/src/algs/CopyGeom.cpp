#include "meshkit/CopyGeom.hpp"
#include "meshkit/MKCore.hpp"
#include "meshkit/SizingFunction.hpp"
#include "meshkit/RegisterMeshOp.hpp"
#include "meshkit/LocalSet.hpp"
#include "meshkit/ModelEnt.hpp"

#include "iMesh_extensions.h"
#include "MBCN.h"


namespace MeshKit
{
// static registration of this  mesh scheme
moab::EntityType CopyGeom_tps[] = { moab::MBVERTEX,
                                    moab::MBEDGE,
                                    moab::MBTRI,
                                    moab::MBHEX,
                                    moab::MBMAXTYPE};
const moab::EntityType* CopyGeom::output_types()
{ return CopyGeom_tps; }

CopyGeom::CopyGeom(MKCore *mkcore, const MEntVector &me_vec)
: MeshScheme(mkcore, me_vec),
  igeomImpl(mkcore->igeom_instance())
{
  m_x[0] = 0.0;
  m_x[1] = 0.0;
  m_x[2] = 0.0;
}

CopyGeom::~CopyGeom()
{}

bool CopyGeom::add_modelent(ModelEnt *model_ent)
{
  return MeshOp::add_modelent(model_ent);
}

void CopyGeom::setup_this()
{}

void CopyGeom::execute_this()
{
  for (MEntSelection::iterator mit = mentSelection.begin();
       mit != mentSelection.end(); mit++) {
    ModelEnt *me = mit->first;
    iBase_EntityHandle ent = reinterpret_cast<iBase_EntityHandle>(
      me->geom_handle());
    iBase_EntityHandle tmp;
    IBERRCHK(igeomImpl->copyEnt(ent, tmp), *igeomImpl);
    IBERRCHK(igeomImpl->moveEnt(tmp, m_x[0], m_x[1], m_x[2]), *igeomImpl);
  }
}

// set the target location - m_x
void CopyGeom::set_location(const Vector<3> &dx)
{
  m_x = dx;
}

void CopyGeom::tag_copied_sets(const char **tag_names, const char **tag_vals,
                               const int num_tags)
{
  // TODO:
}

} // namespace MeshKit
