/** \file test_solidmesher.cpp \test
 *
 * Test the SolidMesher for a basic example.
 *
 */

#include "meshkit/MKCore.hpp"
#include "meshkit/SurfaceFacetMeshReader.hpp"
#include "meshkit/CurveFacetMeshReader.hpp"
#include "meshkit/ModelEnt.hpp"
#include "MBTagConventions.hpp"
#include "moab/GeomTopoTool.hpp"
#include "MBTypes.h"

using namespace MeshKit;

#include "TestUtil.hpp"

MKCore *mk = NULL;

#ifdef HAVE_ACIS
std::string extension = ".sat";
#elif HAVE_OCC
std::string extension = ".stp";
#endif


//test functions 
void cylcube_curve_senses_test();
void cylcube_surf_senses_test();


//support functions 
int geom_id_by_handle( moab::EntityHandle ent );
void load_sat_curve_sense_data( ModelEnt* curve, std::vector<int>& curve_ids_out, std::vector<int>& senses_out );
void load_stp_curve_sense_data( ModelEnt* curve, std::vector<int>& curve_ids_out, std::vector<int>& senses_out );
void load_sat_surf_sense_data( ModelEnt *surf, std::vector<int>& vol_ids_out, std::vector<int>& senses_out );
void load_stp_surf_sense_data( ModelEnt *surf, std::vector<int>& vol_ids_out, std::vector<int>& senses_out );
void check_sense_data( std::vector<moab::EntityHandle> wrt_ents, std::vector<int> senses, 
		       std::vector<int> known_wrt_ids, std::vector<int> known_senses );

int main(int argc, char **argv) 
{
  
  // start up MK and load the geometry
  mk = new MKCore();

  std::string filename = "cylcube" ;

  filename = TestDir + "/" + filename + extension;

  mk->load_geometry(&filename[0]);

  MEntVector surfs;
  mk->get_entities_by_dimension(2,surfs);
  SurfaceFacetMeshReader *sfmr;

  sfmr = (SurfaceFacetMeshReader*) mk->construct_meshop("SurfaceFacetMeshReader", surfs);

  double facet_tol = 1e-04, geom_resabs = 1e-06;
  sfmr->set_mesh_params(facet_tol, geom_resabs);

  mk->setup();
  mk->execute();


  //RUN TESTS
  int num_fail = 0;
  num_fail += RUN_TEST(cylcube_curve_senses_test);


  return num_fail;
}

void cylcube_curve_senses_test()
{

  //check that we retrieve the correct number of curves from the model
  MEntVector curves;
  mk->get_entities_by_dimension( 1, curves); 

#ifdef HAVE_OCC
  CHECK( 18 == int(curves.size()) ); 
#else
  CHECK( 14 == int(curves.size()) ); 
#endif

    // Establish GeomTopoTool instance needed to get curve data 
  moab::GeomTopoTool gt( mk->moab_instance() ); 

    // Initialize vectors for sense checking
  std::vector<moab::EntityHandle> curve_list;
  std::vector<int> senses;  
  std::vector<int> known_curve_ids;
  std::vector<int> known_senses;

  //load_sat_curve_sense_data( 
  //for each curve get the sense data from the model
  for(MEntVector::iterator i = curves.begin(); i != curves.end(); i++)
    {

      ModelEnt *curve = *i;

      curve_list.clear(); 
      senses.clear(); 
      
      //get the senses for this curve
      moab::EntityHandle curve_sh = curve->mesh_handle(); 
      moab::ErrorCode result = gt.get_senses( curve_sh,  curve_list, senses);
      if( result != moab::MB_SUCCESS ) CHECK(false); // return failed test if this call fails      
      
      //clear reference data
      known_curve_ids.clear();
      known_senses.clear(); 

      //Load known curve-sense data
#ifdef HAVE_OCC
      load_stp_curve_sense_data( curve, known_curve_ids, known_senses);
#else
      load_sat_curve_sense_data( curve, known_curve_ids, known_senses); 
#endif
      //check that each surf and sense has a match in our reference data

    } 
}

void cylcube_surf_senses_test()
{

  MEntVector surfs;
  
  mk->get_entities_by_dimension( 2, surfs ); 

  //make sure we've retrieved the right number of surfaces for the model
#ifdef HAVE_OCC
  CHECK_EQUAL( 10, int(surfs.size()) ); 
#else
  CHECK_EQUAL( 9, int(surfs.size()) ); 
#endif

  //Create a Geomtopotool and declare containers for sense data
  moab::GeomTopoTool gt( mk->moab_instance() );
  moab::ErrorCode rval;
  std::vector<moab::EntityHandle> vols; 
  std::vector<int> senses;
  std::vector<int> known_volume_ids; 
  std::vector<int> known_senses; 

  for( unsigned int i = 0; i < surfs.size(); i++)
    {
      //clear out the old containers from the previous surface
      vols.clear();
      senses.clear();

      moab::EntityHandle sh = surfs[i]->mesh_handle(); 

      rval = gt.get_senses( sh, vols, senses);
      CHECK_ERR(rval);
      
      //clear previous reference data
      known_volume_ids.clear();
      known_senses.clear(); 
      //Load the known surface to volume sense data
#ifdef HAVE_OCC
      load_stp_surf_sense_data( surfs[i], known_volume_ids, known_senses);
#else
      load_sat_surf_sense_data( surfs[i], known_volume_ids, known_senses);
#endif

      check_sense_data( vols, senses, known_volume_ids, known_senses);
    }
}


//Loads reference surface to volume sense data into the reference vectors
void load_sat_surf_sense_data( ModelEnt *surf, std::vector<int>& vol_ids_out, std::vector<int>& senses_out ){

  int surf_id = geom_id_by_handle( surf->mesh_handle() );
  switch(surf_id)
  {
    case 1:
          vol_ids_out.push_back(1);
          senses_out.push_back(SENSE_FORWARD); 
          break;

    case 2:
          vol_ids_out.push_back(1);
          senses_out.push_back(SENSE_FORWARD); 
          break;

    case 3:
          vol_ids_out.push_back(1);
          senses_out.push_back(SENSE_FORWARD); 
          break;

    case 4:
          vol_ids_out.push_back(1);
          senses_out.push_back(SENSE_FORWARD); 
          break;

    case 5:
          vol_ids_out.push_back(1);
          senses_out.push_back(SENSE_FORWARD); 
          break;

    case 6:
          vol_ids_out.push_back(1);
          senses_out.push_back(SENSE_FORWARD); 
          break;

    case 7:
          vol_ids_out.push_back(2);
          senses_out.push_back(SENSE_FORWARD);
          break;
  
    case 8:
          vol_ids_out.push_back(2);
          senses_out.push_back(SENSE_FORWARD);
          break;

    case 9:
          vol_ids_out.push_back(2);
          senses_out.push_back(SENSE_FORWARD);
          break;
    default:
      std::cout << "Could not match the surface id to sense data." << std::endl;
      CHECK(false);

   }

}

//Loads reference surface to volume sense data into the reference vectors
void load_stp_surf_sense_data( ModelEnt *surf, std::vector<int>& vol_ids_out, std::vector<int>& senses_out ){

  int surf_id = geom_id_by_handle( surf->mesh_handle() );
  switch(surf_id)
  {
    case 1:
          vol_ids_out.push_back(1);
          senses_out.push_back(SENSE_FORWARD); 
          break;

    case 2:
          vol_ids_out.push_back(1);
          senses_out.push_back(SENSE_FORWARD); 
          break;

    case 3:
          vol_ids_out.push_back(1);
          senses_out.push_back(SENSE_FORWARD); 
          break;

    case 4:
          vol_ids_out.push_back(1);
          senses_out.push_back(SENSE_FORWARD); 
          break;

    case 5:
          vol_ids_out.push_back(1);
          senses_out.push_back(SENSE_FORWARD); 
          break;

    case 6:
          vol_ids_out.push_back(1);
          senses_out.push_back(SENSE_FORWARD); 
          break;

    case 7:
          vol_ids_out.push_back(2);
          senses_out.push_back(SENSE_FORWARD);
          break;
  
    case 8:
          vol_ids_out.push_back(2);
          senses_out.push_back(SENSE_FORWARD);
          break;

    case 9:
          vol_ids_out.push_back(2);
          senses_out.push_back(SENSE_FORWARD);
          break;

    case 10:
          vol_ids_out.push_back(2);
          senses_out.push_back(SENSE_FORWARD);
          break;
    default:
      std::cout << "Could not match surface id to sense data." << std::endl;
      CHECK(false);
   }
}

//Loads two vectors with reference curve and curve_sense data
void load_sat_curve_sense_data( ModelEnt* curve, std::vector<int>& surf_ids_out, std::vector<int>& senses_out )
{

  int curve_id = geom_id_by_handle( curve->mesh_handle() );
  switch(curve_id)
  {
    case 1:
          surf_ids_out.push_back(1); surf_ids_out.push_back(6);
          senses_out.push_back(SENSE_FORWARD); senses_out.push_back(SENSE_REVERSE);
          break;

    case 2:
          surf_ids_out.push_back(1); surf_ids_out.push_back(5);
          senses_out.push_back(SENSE_FORWARD); senses_out.push_back(SENSE_REVERSE);
          break;

    case 3:
          surf_ids_out.push_back(1); surf_ids_out.push_back(4);
          senses_out.push_back(SENSE_FORWARD); senses_out.push_back(SENSE_REVERSE);
          break;

    case 4:
          surf_ids_out.push_back(1); surf_ids_out.push_back(3);
          senses_out.push_back(SENSE_FORWARD); senses_out.push_back(SENSE_REVERSE);
          break;

    case 5:
          surf_ids_out.push_back(2); surf_ids_out.push_back(6);
          senses_out.push_back(SENSE_FORWARD); senses_out.push_back(SENSE_REVERSE);
          break;

    case 6:
          surf_ids_out.push_back(2); surf_ids_out.push_back(3);
          senses_out.push_back(SENSE_FORWARD); senses_out.push_back(SENSE_REVERSE);
          break;

    case 7:
          surf_ids_out.push_back(2); surf_ids_out.push_back(4);
          senses_out.push_back(SENSE_FORWARD); senses_out.push_back(SENSE_REVERSE);
          break;

    case 8:
          surf_ids_out.push_back(2); surf_ids_out.push_back(5);
          senses_out.push_back(SENSE_FORWARD); senses_out.push_back(SENSE_REVERSE);
          break;

    case 9:
          surf_ids_out.push_back(3); surf_ids_out.push_back(4);
          senses_out.push_back(SENSE_FORWARD); senses_out.push_back(SENSE_REVERSE);
          break;

    case 10:
          surf_ids_out.push_back(3); surf_ids_out.push_back(6);
          senses_out.push_back(SENSE_REVERSE); senses_out.push_back(SENSE_FORWARD);
          break;

    case 11:
          surf_ids_out.push_back(4); surf_ids_out.push_back(5);
          senses_out.push_back(SENSE_FORWARD); senses_out.push_back(SENSE_REVERSE);
          break;

    case 12:
          surf_ids_out.push_back(5); surf_ids_out.push_back(6);
          senses_out.push_back(SENSE_FORWARD); senses_out.push_back(SENSE_REVERSE);
          break;

    case 13:
      surf_ids_out.push_back(7); surf_ids_out.push_back(8);
      senses_out.push_back(SENSE_REVERSE); senses_out.push_back(SENSE_FORWARD);
      break;

    case 14:
      surf_ids_out.push_back(7); surf_ids_out.push_back(9);
      senses_out.push_back(SENSE_REVERSE); senses_out.push_back(SENSE_FORWARD);
      break;
    default:
      std::cout << "Could not match curve id to sense data" << std::endl;
      CHECK(false);

  } 

}


//Loads two vectors with reference curve and curve_sense data
void load_stp_curve_sense_data( ModelEnt* curve, std::vector<int>& surf_ids_out, std::vector<int>& senses_out )
{

  int curve_id = geom_id_by_handle( curve->mesh_handle() );
  switch(curve_id)
  {
    case 1:
          surf_ids_out.push_back(1); surf_ids_out.push_back(6);
          senses_out.push_back(SENSE_FORWARD); senses_out.push_back(SENSE_REVERSE);
          break;

    case 2:
          surf_ids_out.push_back(1); surf_ids_out.push_back(5);
          senses_out.push_back(SENSE_FORWARD); senses_out.push_back(SENSE_REVERSE);
          break;

    case 3:
          surf_ids_out.push_back(1); surf_ids_out.push_back(4);
          senses_out.push_back(SENSE_FORWARD); senses_out.push_back(SENSE_REVERSE);
          break;

    case 4:
          surf_ids_out.push_back(1); surf_ids_out.push_back(3);
          senses_out.push_back(SENSE_FORWARD); senses_out.push_back(SENSE_REVERSE);
          break;

    case 5:
          surf_ids_out.push_back(2); surf_ids_out.push_back(6);
          senses_out.push_back(SENSE_FORWARD); senses_out.push_back(SENSE_REVERSE);
          break;

    case 6:
          surf_ids_out.push_back(2); surf_ids_out.push_back(3);
          senses_out.push_back(SENSE_FORWARD); senses_out.push_back(SENSE_REVERSE);
          break;

    case 7:
          surf_ids_out.push_back(2); surf_ids_out.push_back(4);
          senses_out.push_back(SENSE_FORWARD); senses_out.push_back(SENSE_REVERSE);
          break;

    case 8:
          surf_ids_out.push_back(2); surf_ids_out.push_back(5);
          senses_out.push_back(SENSE_FORWARD); senses_out.push_back(SENSE_REVERSE);
          break;

    case 9:
          surf_ids_out.push_back(3); surf_ids_out.push_back(4);
          senses_out.push_back(SENSE_FORWARD); senses_out.push_back(SENSE_REVERSE);
          break;

    case 10:
          surf_ids_out.push_back(3); surf_ids_out.push_back(6);
          senses_out.push_back(SENSE_REVERSE); senses_out.push_back(SENSE_FORWARD);
          break;

    case 11:
          surf_ids_out.push_back(4); surf_ids_out.push_back(5);
          senses_out.push_back(SENSE_FORWARD); senses_out.push_back(SENSE_REVERSE);
          break;

    case 12:
          surf_ids_out.push_back(5); surf_ids_out.push_back(6);
          senses_out.push_back(SENSE_FORWARD); senses_out.push_back(SENSE_REVERSE);
          break;

    case 13:
      surf_ids_out.push_back(7); surf_ids_out.push_back(8);
      senses_out.push_back(SENSE_REVERSE); senses_out.push_back(SENSE_FORWARD);
      break;

    case 14:
      surf_ids_out.push_back(7); surf_ids_out.push_back(9);
      senses_out.push_back(SENSE_REVERSE); senses_out.push_back(SENSE_FORWARD);
      break;
    case 15:
      surf_ids_out.push_back(7); surf_ids_out.push_back(8);
      senses_out.push_back(SENSE_REVERSE); senses_out.push_back(SENSE_FORWARD);
      break;
    case 16:
      surf_ids_out.push_back(7); surf_ids_out.push_back(10);
      senses_out.push_back(SENSE_REVERSE); senses_out.push_back(SENSE_FORWARD);
      break;
    case 17:
      surf_ids_out.push_back(8); surf_ids_out.push_back(10);
      senses_out.push_back(SENSE_REVERSE); senses_out.push_back(SENSE_FORWARD);
      break;
    case 18:
      surf_ids_out.push_back(8); surf_ids_out.push_back(9);
      senses_out.push_back(SENSE_REVERSE); senses_out.push_back(SENSE_FORWARD);
      break;
    default:
      std::cout << "Could not match curve id to sense data" << std::endl;
      CHECK(false);

  } 

}


void check_sense_data( std::vector<moab::EntityHandle> wrt_ents, std::vector<int> senses, 
		       std::vector<int> known_wrt_ids, std::vector<int> known_senses )
{
  
  //Get ID's of the wrt entities
  std::vector<int> wrt_ent_ids;

  for(unsigned int i=0 ; i<wrt_ents.size() ; i++)
  {
      wrt_ent_ids.push_back( geom_id_by_handle( wrt_ents[i] ) );
  }

  for(unsigned int i=0; i< wrt_ent_ids.size() ; i++)
  {
     for(unsigned int j=0; j< known_wrt_ids.size(); j++)
     {
       if( wrt_ent_ids[i] == known_wrt_ids [j] )
         {
          // Make sure the senses of the matching wrt entities
          // are correct
          CHECK_EQUAL( senses[i], known_senses[j] );
          //Once a wrt entity is matched with a known entity,
          // remove it from the list
          wrt_ent_ids.erase( wrt_ent_ids.begin()+i );
          senses.erase( senses.begin()+i );
         }
     }
  }

  // After both loops are complete, known_wrt_ents should be empty 
  int leftovers = wrt_ent_ids.size();
  CHECK_EQUAL( leftovers, 0 );

}


int geom_id_by_handle( moab::EntityHandle ent )
{
  
  moab::ErrorCode rval;
  //Get the id_tag handle
  moab::Tag id_tag;
  rval = mk->moab_instance()->tag_get_handle( GLOBAL_ID_TAG_NAME, 1, MB_TYPE_INTEGER, id_tag, moab::MB_TAG_DENSE );
  CHECK_ERR(rval);

  //Load the ID for the EntHandle given to the function                  
  int id;
  rval = mk->moab_instance()->tag_get_data( id_tag, &ent, 1, &id );                  
  CHECK_ERR(rval);                        

  return id;

 }
  
