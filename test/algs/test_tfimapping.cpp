/** \file test_TFIMapping.cpp
 *
 * Test the TFIMapping for a few challenging examples.
 *
 */

#include "meshkit/MKCore.hpp"
#include "meshkit/TFIMapping.hpp"
#include "meshkit/SizingFunction.hpp"
#include "meshkit/ModelEnt.hpp"
#include "meshkit/Matrix.hpp"
#include "meshkit/EdgeMesher.hpp"

using namespace MeshKit;

#include "TestUtil.hpp"

MKCore *mk = NULL;

void test_TFImapping();
void test_TFImappingcubit();

int main(int argc, char **argv) 
{
  
    // start up MK and load the geometry
  mk = new MKCore();

  int num_fail = 0;
  
  //num_fail += RUN_TEST(test_TFImapping);
  
  num_fail += RUN_TEST(test_TFImappingcubit);
  
}

void test_TFImappingcubit()
{
	std::string file_name = TestDir + "/SquareWithEdgesMeshed.cub";
  	mk->load_geometry(file_name.c_str(), 0, 0, false);

	mk->load_mesh(file_name.c_str());

	// populate mesh to relate geometry entities and mesh sets
  	mk->populate_mesh();

	//check the number of geometrical edges
	MEntVector surfs, curves, loops;
  	mk->get_entities_by_dimension(2, surfs);
	ModelEnt *this_surf = (*surfs.rbegin());

	this_surf->get_adjacencies(1, curves);

	CHECK_EQUAL(4, (int)curves.size());
	
	//check the number of mesh line segments
	moab::Range edges;
	moab::ErrorCode rval = mk->moab_instance()->get_entities_by_dimension(0, 1, edges);
	CHECK_EQUAL(moab::MB_SUCCESS, rval);
	CHECK_EQUAL(40, (int)edges.size());

	//now, do the TFIMapping
	TFIMapping *tm = (TFIMapping*)mk->construct_meshop("TFIMapping", surfs);
	mk->setup_and_execute();

	//check the number of quads
	moab::Range faces;
	rval = mk->moab_instance()->get_entities_by_dimension(0, 2, faces);
	CHECK_EQUAL(moab::MB_SUCCESS, rval);
	CHECK_EQUAL(100, (int)faces.size());

	mk->save_mesh("TFIMappingFromCubit.vtk");

	delete tm;
	mk->clear_graph();
	
}

void test_TFImapping() 
{
  std::string file_name = TestDir + "/SquareWithoutMesh.cub";
  mk->load_geometry(file_name.c_str());

    // get the surface
  MEntVector surfs, curves, loops;
  mk->get_entities_by_dimension(2, surfs);
  CHECK_EQUAL(1, (int)surfs.size());
  
    // make an edge mesher
  mk->get_entities_by_dimension(1, curves);
  //test there are 4 edges bounding the surface
  CHECK_EQUAL(4, (int)curves.size());
  EdgeMesher *em = (EdgeMesher*) mk->construct_meshop("EdgeMesher", curves);

    // make a sizing function and set it on the surface
  SizingFunction esize(mk, 10, -1);
  surfs[0]->sizing_function_index(esize.core_index());
  
  // mesh the edges, by calling execute
  mk->setup_and_execute();

    // make sure we got the right number of edges
  moab::Range edges, nodes;
  moab::ErrorCode rval = mk->moab_instance()->get_entities_by_dimension(0, 1, edges);
  CHECK_EQUAL(moab::MB_SUCCESS, rval);
  CHECK_EQUAL(40, (int)edges.size());

  rval = mk->moab_instance()->get_entities_by_dimension(0, 0, nodes);
  CHECK_EQUAL(moab::MB_SUCCESS, rval);
  CHECK_EQUAL(40, (int)nodes.size());


/*
  ///////////////////////////////////////////////////////////////////////////////////////////
  // test
  
  std::cout << "test function: test_tfimapping.cpp " << std::endl; 
  iBase_EntitySetHandle geom_root_set;
  iGeom::Error g_err;
  geom_root_set = mk->igeom_instance()->getRootSet();
  
  std::vector<iBase_EntityHandle> gEdges;
  g_err = mk->igeom_instance()->getEntities(geom_root_set, iBase_EDGE, gEdges);
  IBERRCHK(g_err, "Trouble get the geometrical edges.");
  assert(gEdges.size()==4);

  for (unsigned int i = 0; i < gEdges.size(); i++)
  {
  std::vector<iBase_EntityHandle> gVertices;
  g_err = mk->igeom_instance()->getEntAdj(gEdges[i], iBase_VERTEX, gVertices);
  IBERRCHK(g_err, "Trouble get the adjacent geometric nodes on an edge.");
    
  //test one mesh node on the geometrical edge
  iBase_EntitySetHandle TmpSet;
  std::vector<iBase_EntityHandle> tmpNodeHandle;
  iRel::Error r_err = mk->irel_pair()->getEntSetRelation(gVertices[0], 0, TmpSet);
  IBERRCHK(r_err, "Trouble get the mesh entity set from the geometrical corner 0.");
  iMesh::Error m_err = mk->imesh_instance()->getEntities(TmpSet, iBase_VERTEX, iMesh_POINT, tmpNodeHandle);	
  IBERRCHK(m_err, "Trouble get the mesh node from the geometrical corner 0.");
  assert(tmpNodeHandle.size()==1);	
  double coord[3];
  m_err = mk->imesh_instance()->getVtxCoord(tmpNodeHandle[0], coord[0], coord[1], coord[2]);
  IBERRCHK(m_err, "Trouble get the coordinates from the mesh node.");
  std::cout << "one end node is x = " << coord[0] << "\ty = " << coord[1] << "\tz = " << coord[2] << std::endl;
    
  //test the other mesh node on the geometrical edge
  r_err = mk->irel_pair()->getEntSetRelation(gVertices[1], 0, TmpSet);
  IBERRCHK(r_err, "Trouble get the mesh entity set from the geometrical corner 0.");
  tmpNodeHandle.clear();
  m_err = mk->imesh_instance()->getEntities(TmpSet, iBase_VERTEX, iMesh_POINT, tmpNodeHandle);	
  IBERRCHK(m_err, "Trouble get the mesh node from the geometrical corner 0.");
  assert(tmpNodeHandle.size()==1);	
  m_err = mk->imesh_instance()->getVtxCoord(tmpNodeHandle[0], coord[0], coord[1], coord[2]);
  IBERRCHK(m_err, "Trouble get the coordinates from the mesh node.");
  std::cout << "The other end node is x = " << coord[0] << "\ty = " << coord[1] << "\tz = " << coord[2] << std::endl;
  
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //test the edge mesh after EdgeMesher
  std::vector<iBase_EntityHandle> tmpEdgeHandle;
  r_err = mk->irel_pair()->getEntSetRelation(gEdges[i], 0, TmpSet);
  IBERRCHK(r_err, "Trouble get the mesh entity set from the geometrical edge.");
  tmpNodeHandle.clear();
  m_err = mk->imesh_instance()->getEntities(TmpSet, iBase_VERTEX, iMesh_POINT, tmpNodeHandle);	
  IBERRCHK(m_err, "Trouble get the mesh node from the geometrical edge.");
  std::cout << "node size on the geometrical edge is " << tmpNodeHandle.size() << std::endl;
  for (unsigned int i = 0; i < tmpNodeHandle.size(); i++)
  {
	m_err = mk->imesh_instance()->getVtxCoord(tmpNodeHandle[i], coord[0], coord[1], coord[2]);
	IBERRCHK(m_err, "Trouble get the mesh node from the geometrical edge.");
	std::cout << "i = " << i << "\tx=" << coord[0] << "\ty=" << coord[1] << "\tz=" << coord[2] << std::endl;

  }
  //test number of line segments
  m_err = mk->imesh_instance()->getEntities(TmpSet, iBase_EDGE, iMesh_LINE_SEGMENT, tmpEdgeHandle);
  IBERRCHK(m_err, "Trouble get the mesh line segments from the geometrical edge.");
  std::cout << "number of line segments is " << tmpEdgeHandle.size() << std::endl;
   //test how many mesh entity sets are there on the geometrical edge
   
  

   //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  }

  //////////////////////////////////////////////////////////////////////////////////////////  
*/


  //ok, we are done with edge mesher
  //now, do the TFIMapping
  TFIMapping *tm = (TFIMapping*)mk->construct_meshop("TFIMapping", surfs);
  mk->setup_and_execute();
  
  //mk->populate_mesh();

  //check whether we got the right number of quads after TFIMapping
  moab::Range faces;
  rval = mk->moab_instance()->get_entities_by_dimension(0, 2, faces);
  CHECK_EQUAL(moab::MB_SUCCESS, rval);
  CHECK_EQUAL(100, (int)faces.size());

  //output the mesh to vtk file
  mk->save_mesh("TFIMapping.vtk");

    // clean up
  delete em;
  delete tm;
  //delete mk->vertex_mesher();
  mk->clear_graph();
}
