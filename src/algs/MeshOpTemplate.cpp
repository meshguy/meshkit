#include <stdio.h>
#include "meshkit/MeshOpTemplate.hpp"
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
    //igeomImpl(mkcore->igeom_instance(), mkmoab(mkcore->moab_instance())
  {
    s_x[0] = 1.0;
    s_x[1] = 1.0;
    s_x[2] = 1.0;
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
    igeomImpl->createBrick(s_x[0], s_x[1], s_x[2], argf);
    printf("Made a brick!\n");

    iBase_EntitySetHandle seth;
    iBase_TagHandle mattag;
    std::string brick = "opbrick";
    std::string name_tag_id = "NAME";
    // testing to see if this propages to mesh and how to use iRel for relating mesh and geometry mp.brep
    std::string mtag_name = "MATERIAL_TAG";
    iBase_TagHandle name_tag;
    igeomImpl->getTagHandle(name_tag_id.c_str(),name_tag);

    // create material_set tag
    igeomImpl->createTag(mtag_name.c_str(),1, iBase_INTEGER, mattag);

    igeomImpl->createEntSet(false, seth);
    igeomImpl->addEntToSet(argf, seth);

    igeomImpl->setEntSetData(seth,name_tag,brick.c_str());

    int matid = 111555;
    igeomImpl->setData(argf,mattag,(char*)(&matid));
    // This does not work?
    //igeomImpl->setEntSetData(seth,mattag, (char*)(&matid));
    // Meshing this with test_ngtetmesher does NOT produce an material set
    bool save = true;
    if (save == true)
      igeomImpl->save("mp.brep");

  }

  // set the size of the brick - s_x
  void MeshOpTemplate::set_size(const Vector<3> &dx)
  {
    s_x = dx;
  }

  void MeshOpTemplate::tag_copied_sets(const char **tag_names, const char **tag_vals,
                                       const int num_tags)
  {
    // TODO:
  }

} // namespace MeshKit
