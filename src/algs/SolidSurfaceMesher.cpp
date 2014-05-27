

#include "meshkit/MKCore.hpp"
#include "meshkit/SolidSurfaceMesher.hpp"
#include "meshkit/ModelEnt.hpp"
#include <iostream>
#include <iGeom.h>

#include "MBTagConventions.hpp"
#include "MBCore.hpp"
#include "moab/GeomTopoTool.hpp"


namespace MeshKit
{

// Output mesh types for this class
moab::EntityType SolidSurfaceMesher_types[] = {moab::MBVERTEX, moab::MBEDGE, moab::MBTRI, moab::MBMAXTYPE};
const moab::EntityType* SolidSurfaceMesher::output_types()
{ return SolidSurfaceMesher_types; }

SolidSurfaceMesher::SolidSurfaceMesher(MKCore *mk_core, const MEntVector &ments)
  : MeshScheme(mk_core, ments)
{
  mk = mk_core; 
  model_ents = ments;
}

SolidSurfaceMesher::~SolidSurfaceMesher()
{
}

void SolidSurfaceMesher::setup_this()
{

  //create a solid curve mesher 
  MeshOp *scm = (MeshOp*) mk_core()->construct_meshop("CurveMesher");
 
  MEntVector::iterator i; 
  for(i = model_ents.begin(); i != model_ents.end(); i++)
    {
      ModelEnt *me = *i;
      //do an initial check that the model entity is of the correct dimension
      if(me->dimension() != 2)
	{
	  std::cout << "Found a model entity with an incorrect dimension for this mesher" << std::endl;
          model_ents.erase(i); 
          continue; 
	}

      //make sure that the model entity's children have been meshed
      MEntVector children; 
      me->get_adjacencies(1, children); 
 
      for(MEntVector::iterator j = children.begin(); j != children.end(); j++)
	{
	  ModelEnt *child = *j; 
          if(child->is_meshops_list_empty()) scm-> add_modelent(child); 
	}

    }

  //have the children be meshed first
  mk_core()->insert_node(scm, this, mk_core()->root_node());

  //set the faceting parameters to default values if they are not already set
  
  //set the faceting tolerance
  if(!facet_tol) facet_tol = 1e-4; 
  
  //set the geom_resabs value
  if(!geom_res) geom_res = 1e-6;

}

void SolidSurfaceMesher::execute_this()
{
}

}
