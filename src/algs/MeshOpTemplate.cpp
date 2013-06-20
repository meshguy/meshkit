#include <stdio.h>
#include "meshkit/MeshOpTemplate.hpp"
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
moab::EntityType MeshOpTemplate_tps[] = { moab::MBVERTEX,
                                    moab::MBEDGE,
                                    moab::MBTRI,
                                    moab::MBHEX,
                                    moab::MBMAXTYPE};
const moab::EntityType* MeshOpTemplate::output_types()
{ return MeshOpTemplate_tps; }

MeshOpTemplate::MeshOpTemplate(MKCore *mkcore, const MEntVector &me_vec)
: MeshScheme(mkcore, me_vec),
  igeomImpl(mkcore->igeom_instance())
{
  m_x[0] = 1.0;
  m_x[1] = 1.0;
  m_x[2] = 1.0;
}

MeshOpTemplate::~MeshOpTemplate()
{}

bool MeshOpTemplate::add_modelent(ModelEnt *model_ent)
{
  return MeshOp::add_modelent(model_ent);
}

void MeshOpTemplate::setup_this()
{
}

void MeshOpTemplate::execute_this()
{
  iBase_EntityHandle_Private * argf = NULL;
  igeomImpl->createBrick(m_x[0], m_x[1], m_x[2], argf);
  printf("Made a brick!\n");
}

// set the target location - m_x
void MeshOpTemplate::set_size(const Vector<3> &dx)
{
  m_x = dx;
}

void MeshOpTemplate::tag_copied_sets(const char **tag_names, const char **tag_vals,
                               const int num_tags)
{
  // TODO:
}

} // namespace MeshKit
