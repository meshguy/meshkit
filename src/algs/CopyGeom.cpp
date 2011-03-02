#include "meshkit/CopyGeom.hpp"
#include "meshkit/MKCore.hpp"
#include "meshkit/SizingFunction.hpp"
#include "meshkit/RegisterMeshOp.hpp"
#include "meshkit/LocalSet.hpp"
#include "meshkit/ModelEnt.hpp"

#include "CopyUtils.hpp"

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
    std::vector<iBase_EntityHandle> orig_ents(mentSelection.size()), new_ents(mentSelection.size());
	int i = -1;
    for (MEntSelection::iterator mit = mentSelection.begin();
         mit != mentSelection.end(); mit++) {
	  i++;
      ModelEnt *me = mit->first;
      orig_ents[i++] =  reinterpret_cast<iBase_EntityHandle>( me->geom_handle());
    }
    copy(&orig_ents[0], i, m_x);
}

// set the target location - m_x
void CopyGeom::set_location(const double x[])
{
	m_x[0] = x[0];
	m_x[1] = x[1];
	m_x[2] = x[2];
 }

// copy move entities to m_x
void CopyGeom::copy(iBase_EntityHandle *entities, const int entities_ehsize,
		    const double *dx)
{
  for (int i = 0; i < entities_ehsize; i++){
	  iBase_EntityHandle temp = NULL;
	  IBERRCHK(igeomImpl->copyEnt(entities[i], temp), *igeomImpl);

	  IBERRCHK(igeomImpl->moveEnt(temp, dx[0], dx[1], dx[2]), *igeomImpl);
  }
}

void CopyGeom::tag_copied_sets(const char **tag_names, const char **tag_vals,
                       const int num_tags)
{
  // TODO:

}

} // namespace MeshKit
