
#include "meshkit/MKCore.hpp"
#include "meshkit/CurveMesher.hpp"
#include "meshkit/ModelEnt.hpp"
#include <iostream>
#include <iGeom.h>

namespace MeshKit
{
// Construction Function for CurveMesher

moab::EntityType CurveMesher_types[] = {moab::MBVERTEX, moab::MBEDGE, moab::MBMAXTYPE};
const moab::EntityType* CurveMesher::output_types()
  { return CurveMesher_types; }


CurveMesher::CurveMesher(MKCore *mk_core, const MEntVector &ments)
  : MeshScheme(mk_core, ments)
{
  mk = mk_core;
  model_ents = ments;
}

// Destructor Function for Curvemesher
CurveMesher::~CurveMesher()
{
}

void CurveMesher::setup_this()
{

  MEntVector::iterator i;
  for ( i = model_ents.begin(); i != model_ents.end(); i++)
    { 
      ModelEnt *me = *i;
      //do an initial check that the model entity is of the correct dimension
      if( me->dimension() != 1)
	{
	  std::cout << "Found an entity that is of an incorrect dimension" << std::endl;
          model_ents.erase(i);
          continue;
        } 

    }
}

void CurveMesher::execute_this()
{

}


}
