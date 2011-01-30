#include "CopyGeom.hpp"
#include "MergeMesh.hpp"
#include "TestFramework.hpp"
#include <cfloat>
#include "SimpleArray.hpp"

// This test creates 4 prisms (X) next to each other (XXXX).
// Two prism are created by create command, 
// and the two other prisms are copy moved using CopyGeom class.


int main(int argc, char **argv) 
{
  int err=0, num_geoms=2;
  iGeom_Instance geom;
  iBase_EntitySetHandle root_set, orig_set1, orig_set2, set2;
  CopyGeom **cg;
  std::vector<iBase_EntitySetHandle> assys;
  SimpleArray<iBase_EntityHandle> entities, entities1, entities2, check_ents, new_ents;
  iBase_EntityHandle assm1 = NULL, assm2 = NULL;
  double dHeight = 10.0, dSide = 1.0;

  // make a geom instance
  iGeom_newGeom("GEOM", &geom, &err, 4);
  CHECK_ERR("Couldn't create mesh.");

  iGeom_getRootSet(geom, &root_set, &err);
  CHECK_ERR("Couldn't create mesh.");

  // create cm instances for each geometry
  cg = new CopyGeom*[num_geoms];
  for (unsigned int i = 0; i < num_geoms; i++) {
    cg[i] = new CopyGeom(geom);
  }

  //  create geometry on the same igeom instance and seperate into sets.
  // First Geometry
  iGeom_createEntSet(geom, 0, &orig_set1, &err);
  CHECK_ERR("Couldn't create mesh.");

  iGeom_createPrism(geom, dHeight, 6, 
		    dSide, dSide,
		    &assm1, &err); 
  CHECK_ERR("Prism creation failed.");


  iGeom_getEntities( geom, root_set, iBase_REGION,
		     ARRAY_INOUT(entities1), &err );
  CHECK_ERR( "Problems getting entity." );

  // move this to the right, to create space for the next prism
  double dx[3];
  dx[0] = 1.732;
  dx[1] = 0.0;
  dx[2] = 0.0;

  for (int i =0; i < entities1.size(); i++){
    iGeom_moveEnt(geom, entities1[i], dx[0], dx[1], dx[2], &err);
    CHECK_ERR("Failed to move geometries.");

    // Also, store the entities in a set
    iGeom_addEntToSet(geom, entities1[i], orig_set1, &err);
    CHECK_ERR( "Problems adding entity to set." );
  }

  // store the geom set corresponding to the first geometry
  assys.push_back(orig_set1);

  // Second Geometry

  iGeom_createEntSet(geom, 0, &orig_set2, &err);
  CHECK_ERR("Couldn't create mesh.");

  iGeom_createPrism(geom, dHeight, 6, 
		    dSide, dSide,
		    &assm2, &err); 
  CHECK_ERR("Prism creation failed.");

  iGeom_getEntities( geom, root_set, iBase_REGION,
		     ARRAY_INOUT(entities2),  &err );
  CHECK_ERR( "Problems getting entity." );


  // add the entity


  for (int i =0; i < entities2.size(); i++){
    iGeom_addEntToSet(geom, entities2[i], orig_set2, &err);
    CHECK_ERR( "Problems adding entity to set." );
  }

  iGeom_createEntSet(geom, 0, &set2, &err);
  CHECK_ERR("Couldn't create ent set.");

  iGeom_subtract (geom, orig_set2, orig_set1, &set2, &err);
  CHECK_ERR( "Problem getting entities." );

  iGeom_getEntities( geom, set2, iBase_REGION,
		     ARRAY_INOUT(check_ents),  &err );
  CHECK_ERR( "Problem getting entities." );

  // store the geom set corresponding to the second geometry
  assys.push_back(set2);

  double dx_cm[3];
  dx_cm[0] = 2*1.732;
  dx_cm[1] = 0.0;
  dx_cm[2] = 0.0;

  // now copy/move the geometries using the set.
  for(int i=0; i < num_geoms; i++){
    iGeom_getEntities( geom, assys[i], iBase_ALL_TYPES,
  		       ARRAY_INOUT(entities), &err );
    CHECK_ERR("Failed to get entities from set.");

    cg[i]->copy(ARRAY_IN(entities), dx_cm, ARRAY_INOUT(new_ents));
    CHECK_ERR("Failed to copy/move.");
  }
  
#ifdef TESTSAVE
  iGeom_save(geom, "c.sat", NULL, &err, 5, 0);
  CHECK_ERR("Couldn't save geometry.");
#endif
  std::cout << "CopyGeom_test passed" << std::endl;
  return 0;
}
