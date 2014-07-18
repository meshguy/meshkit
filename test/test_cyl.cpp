///
/// Patrick Shriwise
/// September 2013
/// This program is designed to run a set of tests on the make_watertight algorithm.
/// This will be a stand-alone piece of code that uses MOAB to open, modify
/// (break), and re-seal geometries/ 
/// input: cyl.h5m file (found in ../make_watertight/test/)
/// output: pass/fail for each of the tests


#include <iostream>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <vector>
#include <set>
#include <algorithm>
#include "MBCore.hpp"
#include "MBTagConventions.hpp"
#include "MBRange.hpp"
#include "MBSkinner.hpp"

#include "mw_func.hpp"
#include "cw_func.hpp"
#include "gen.hpp"
#include "arc.hpp"
#include "zip.hpp"

MBErrorCode move_vert( MBEntityHandle vertex, double dx, double dy, double dz, bool verbose = false );
MBErrorCode rand_vert_move( MBEntityHandle vertex, double tol, bool verbose = false);

MBInterface *MBI();

// struct to hold coordinates of skin edge, it's surface id, and a matched flag
struct coords_and_id {
  double x1;
  double y1;
  double z1;
  double x2;
  double y2;
  double z2;
  int  surf_id;
  bool matched;
  MBEntityHandle vert1;
  MBEntityHandle vert2;
};


// used to clear all mesh data and reload the file as original
MBErrorCode reload_mesh(const char* filename,  MBEntityHandle &meshset, bool debug = false) {

MBErrorCode result;

  // delete meshset
  result= MBI() -> delete_mesh();
  if(gen::error(MB_SUCCESS!=result, "could not delete the mesh set")) return result;



  // re-initialize meshset
  result= MBI() -> create_meshset(MESHSET_SET, meshset);
  if(gen::error(MB_SUCCESS!=result, "could not create the new meshset")) return result; 
 
  //reload the file
  result = MBI() -> load_file(filename, &meshset);
  if(gen::error(MB_SUCCESS!=result, "could not re-load the file into the mesh set")) return result; 

  //check that something was loaded into the meshset
  int num_meshsets;   

  result= MBI() -> num_contained_meshsets ( meshset, &num_meshsets);
  if(gen::error(MB_SUCCESS!=result, "could not get the number of meshsets")) return result; 
  
  if(debug) std::cout << "number of mesh sets after reload = " << num_meshsets << std::endl;

  if(num_meshsets==0) return MB_FAILURE;

  return result;
}


// bumps the last vertex in the model by the x,y,z values given to the problem 
MBErrorCode single_vert_bump( MBRange verts, double bump_dist_x, double bump_dist_y, double bump_dist_z,  std::string root_name, bool verbose = false ) {
 
  MBErrorCode result; 
  MBEntityHandle vertex=verts.back();
  //move the desired vertex by the allotted distance
  result = move_vert( vertex, bump_dist_x, bump_dist_y, bump_dist_z);
  if(gen::error(MB_SUCCESS!=result, "could not move single vert")) return result;


  // write mesh to a new file (might not be necesarry any longer as we move to doing tests on a moab-instance basis)
  // alter output filename
  std::string output_filename = root_name + "_mod.h5m";

  //write file
  result=MBI()->write_mesh(output_filename.c_str());
  if(gen::error(MB_SUCCESS!=result,"could not write the mesh to output_filename")) return result; 

  return MB_SUCCESS;
}

// moves the last vertex in the model in a random direction by a distance less than the faceting tolerance
MBErrorCode rand_vert_bump( MBRange verts, double facet_tol, std::string root_name, bool verbose = false ) {
 
  MBErrorCode result; 
  MBEntityHandle vertex=verts.back();
  //move the desired vertex by the allotted distance
  result = rand_vert_move( vertex, facet_tol, verbose);
  if(gen::error(MB_SUCCESS!=result, "could not move single vert")) return result;


  // write mesh to a new file (might not be necesarry any longer as we move to doing tests on a moab-instance basis)
  // alter output filename
  std::string output_filename = root_name + "_mod.h5m";

  //write file
  result=MBI()->write_mesh(output_filename.c_str());
  if(gen::error(MB_SUCCESS!=result,"could not write the mesh to output_filename")) return result; 

  return MB_SUCCESS;
}

// FOR CYLINDER TESTING ONLY 
// moves the last vertex in the model along the curve of the cylinder 
MBErrorCode theta_vert_bump( MBRange verts, double bump_dist_theta, double tolerance, std::string root_name, bool verbose = false ) {
 
  MBErrorCode result; 
  
  //get vertex coordinates
  double coords[3];
  MBEntityHandle vertex=verts.back();
  result= MBI()-> get_coords( &vertex , 1 , coords );
 
  // determine radius
  double radius=sqrt(pow(coords[0],2)+pow(coords[1],2));
  //std::cout << "radius = " << radius << std::endl;

  // get the current theta value
  double theta=asin(coords[1]/radius);
  //std::cout << "theta = " << theta << std::endl;
  
  // set the vertex bump distance
  double dtheta = tolerance/(radius);
  //std::cout << "dtheta = " << dtheta << std::endl;

 
 if(verbose)
 {
  //original coordinates
  std::cout << std::endl << "Original Coordinates" << std::endl;
  std::cout << "x = " << coords[0] << std::endl;
  std::cout << "y = " << coords[1] << std::endl;
  std::cout << "z = " << coords[2] << std::endl;

 // create new x and y values
  coords[0]=radius*cos(theta+dtheta); 
  std::cout << "new x value = " << coords[0] << std::endl;
  coords[1]=radius*cos(theta+dtheta);
  std::cout << "new y value = " << coords[1] << std::endl;
  
  //altered coordinates
  std::cout << std::endl << "Modified Coordinates" << std::endl;
  std::cout << "x = " << coords[0] << std::endl;
  std::cout << "y = " << coords[1] << std::endl;
  std::cout << "z = " << coords[2] << std::endl;
 }
  //write new coordinates to the mesh
  // might not be necesarry any longer as we move to doing tests on a moab-instance basis
  result=MBI()-> set_coords( &vertex, 1, coords);
  if (gen::error(MB_SUCCESS!=result, "could not set the vertex coordinates")) return result;
  // alter output filename
  std::string output_filename = root_name + "_mod.h5m";
  //write file
  result=MBI()->write_mesh(output_filename.c_str());
  if(gen::error(MB_SUCCESS!=result,"could not write the mesh to output_filename")) return result; 

  return MB_SUCCESS;
}

// takes an entity handle vertex, gets the original coordinates, changes and resets the vertex coordinates in the mesh
MBErrorCode move_vert( MBEntityHandle vertex, double dx, double dy, double dz, bool verbose) {

 MBErrorCode result;
 
 //get coordinates from the mesh
 double coords[3];
 result= MBI()-> get_coords( &vertex , 1 , coords );
 if(gen::error(MB_SUCCESS!=result, "could not get the vertex coordinates")) return result; 

 if(verbose)
 {
  //original coordinates
  std::cout << std::endl << "Original Coordinates" << std::endl;
  std::cout << "x = " << coords[0] << std::endl;
  std::cout << "y = " << coords[1] << std::endl;
  std::cout << "z = " << coords[2] << std::endl;

  coords[0]+=dx;
  coords[1]+=dy;
  coords[2]+=dz;
  
  //altered coordinates
  std::cout << std::endl << "Modified Coordinates" << std::endl;
  std::cout << "x = " << coords[0] << std::endl;
  std::cout << "y = " << coords[1] << std::endl;
  std::cout << "z = " << coords[2] << std::endl;
 }
 
  //write new coordinates to the mesh
  result=MBI()-> set_coords( &vertex, 1, coords);
  if (gen::error(MB_SUCCESS!=result, "could not set the vertex coordinates")) return result;

return MB_SUCCESS;
}

// takes an entity handle vertex, gets the original coordinates, changes and resets the vertex coordinates in the mesh
MBErrorCode rand_vert_move( MBEntityHandle vertex, double tol, bool verbose) {

 MBErrorCode result;
 
 //get coordinates from the mesh
 double coords[3];
 result= MBI()-> get_coords( &vertex , 1 , coords );
 if(gen::error(MB_SUCCESS!=result, "could not get the vertex coordinates")) return result; 

 // get random values for the changes in x,y,z
 double dx,dy,dz;
 dx=rand(); dy=rand(); dz=rand();

 double mag= sqrt(pow(dx,2)+pow(dy,2)+pow(dz,2));

 // set the change in the vertex to be a unit vector
 dx/=mag; dy/=mag; dz/=mag;

 // set the change in the vertex to be something slightly less than the facet tolerance
 dx*=tol*0.9; dy*=tol*0.9; dz*=tol*0.9;


 if(verbose)
 {
  //original coordinates
  std::cout << std::endl << "Original Coordinates" << std::endl;
  std::cout << "x = " << coords[0] << std::endl;
  std::cout << "y = " << coords[1] << std::endl;
  std::cout << "z = " << coords[2] << std::endl;

  coords[0]+=dx;
  coords[1]+=dy;
  coords[2]+=dz;
  
  //altered coordinates
  std::cout << std::endl << "Modified Coordinates" << std::endl;
  std::cout << "x = " << coords[0] << std::endl;
  std::cout << "y = " << coords[1] << std::endl;
  std::cout << "z = " << coords[2] << std::endl;
 }
 
  //write new coordinates to the mesh
  result=MBI()-> set_coords( &vertex, 1, coords);
  if (gen::error(MB_SUCCESS!=result, "could not set the vertex coordinates")) return result;

return MB_SUCCESS;
}


int main(int argc, char **argv)
{

//open moab instance
MBInterface *MBI();

//for unit testing purposes, we don't care about the output. Just PASS or FAIL. 
bool verbose=false;

// ******************************************************************
  // Load the h5m file and create tags.
  // ******************************************************************

  clock_t start_time;
  start_time = clock();
  //save the mesh to a new filename
  std::string input_name=argv[1];
  std::string root_name=argv[1];
  int len = root_name.length();
  root_name.erase(len-4);
  // check input args
  if( argc < 2 || argc > 5) 
    {
    std::cout << "To check using topology of facet points:              " << std::endl;
    std::cout << "./test_cyl <filename> <verbose(true or false)>" << std::endl;
    std::cout << "To check using geometry tolerance of facet points:    " << std::endl;
    std::cout << "./test_cyl <filename> <verbose(true or false)> <tolerance>" << std::endl;
    return 1;
    }

  // load file and get tolerance from input argument
  MBErrorCode result;
  std::string filename = argv[1]; //set filename
  MBEntityHandle input_set;
  result = MBI()->create_meshset( MESHSET_SET, input_set ); //create handle to meshset
  if(MB_SUCCESS != result) 
    {
      return result;
    }

result = MBI()->load_file( filename.c_str(), &input_set ); //load the file into the meshset
  if(MB_SUCCESS != result) 
    {
      // failed to load the file
      std::cout << "could not load file" << std::endl;
      return result;
    }

  /// get faceting tolerance ////
  double tol;

  MBTag faceting_tol_tag;
  //get faceting tolerance handle from file
  result = MBI()->tag_get_handle( "FACETING_TOL", 1, MB_TYPE_DOUBLE,
        faceting_tol_tag , moab::MB_TAG_SPARSE|moab::MB_TAG_CREAT );
  if(gen::error(MB_SUCCESS!=result, "could not get the faceting tag handle")) return result;
  
  //get the faceting tolerance of any entity
  MBRange file_set;
  result = MBI()->get_entities_by_type_and_tag( 0, MBENTITYSET, 
                        &faceting_tol_tag, NULL, 1, file_set );

  //get facetint tolerance value
  result = MBI()->tag_get_data( faceting_tol_tag, &file_set.front(), 1, &tol );
  if(gen::error(MB_SUCCESS!=result, "could not get the faceting tolerance")) return result; 

  // create tags on geometry
  MBTag geom_tag, id_tag;
  result = MBI()->tag_get_handle( GEOM_DIMENSION_TAG_NAME, 1, 
                            MB_TYPE_INTEGER, geom_tag, moab::MB_TAG_DENSE|moab::MB_TAG_CREAT );
  if(gen::error(MB_SUCCESS != result, "could not get GEOM_DIMENSION_TAG_NAME handle")) return result; 

  result = MBI()->tag_get_handle( GLOBAL_ID_TAG_NAME, 1, 
                            MB_TYPE_INTEGER, id_tag, moab::MB_TAG_DENSE|moab::MB_TAG_CREAT );
  if(gen::error(MB_SUCCESS != result, "could not get GLOBAL_ID_TAG_NAME handle")) return result;
  
  // get surface and volume sets
  MBRange surf_sets, vol_sets; // MBRange of set of surfaces and volumes
  // surface sets
  int dim = 2;
  void* input_dim[] = {&dim};
  result = MBI()->get_entities_by_type_and_tag( input_set, MBENTITYSET, &geom_tag, 
                                                input_dim, 1, surf_sets);
  if(MB_SUCCESS != result) 
    {
      return result;
    }

  // volume sets
  dim = 3;
  result = MBI()->get_entities_by_type_and_tag( input_set, MBENTITYSET, &geom_tag, 
                                                input_dim, 1, vol_sets);
  if(MB_SUCCESS != result)
    {
      return result;
    }
    
  //vertex sets
  dim= 0;
  MBRange verts;
  result = MBI()->get_entities_by_dimension(input_set, dim, verts, false);
  if(gen::error(MB_SUCCESS!=result, " could not get vertices from the mesh")) return result;
  
  if(gen::error(MB_SUCCESS!=result, "could not get vertex coordinates")) return result; 

if(verbose)
{
  std::cout<< "number of verticies= " << verts.size() << std::endl;  
  std::cout<< "number of surfaces= " << surf_sets.size() << std::endl;
  std::cout<< "number of volumes= "  << vol_sets.size() << std::endl;
}

///////////////BEGIN 1st TEST////////////////////////

 std::cout << "Test 1: vertex bump in z direction:";

  //perform single vertex bump and test
  result=single_vert_bump(verts, 0 , 0 , 10*tol, root_name, verbose ); 
  if(gen::error(MB_SUCCESS!=result, "could not create single vertex bump test")) return result;
 
  
  bool check_topology, test, sealed;
  check_topology=false;
  test=true;
  
  /*
  //check that geometry is unsealed
  result=cw_func::check_mesh_for_watertightness( input_set, tol, sealed, test);
  if(gen::error(MB_SUCCESS!=result, "could not check model for watertightness")) return result;
  */

  if(sealed)
  {
   std::cout << "Warning: geometry was not broken by single bump vert" << std::endl;
  }
  
  //seal the mesh
  double facet_tol;
  result=mw_func::make_mesh_watertight (input_set, facet_tol, verbose);
  if(gen::error(MB_SUCCESS!=result, "could not make the mesh watertight")) return result;
  
  // Lastly Check to see if make_watertight fixed the model
  result=cw_func::check_mesh_for_watertightness( input_set, tol, sealed, test);
  if(gen::error(MB_SUCCESS!=result, "could not check model for watertightness")) return result;
  
  if(sealed)
  {
   std::cout << "PASS" << std::endl;
  }
  else
  {
   std::cout << "FAIL" << std::endl;
   exit(0);
  }
 

///////////////BEGIN 2ND TEST////////////////////////

  std::cout << "Test 2: vertex bump in y direction:";

  // Clear the mesh and reload original geometry for the next test
  result=reload_mesh( filename.c_str(), input_set);
  if (gen::error(MB_SUCCESS!=result, "could not reload the mesh" )) return result; 

  // retrieve the verticies again so the model can be broken
  result = MBI()->get_entities_by_dimension(input_set, dim, verts, false);
  if(gen::error(MB_SUCCESS!=result, " could not get vertices from the mesh")) return result;

  // "Break" the geometry
  result=single_vert_bump(verts, 0 , 10*tol , 0, root_name ); 
  if(gen::error(MB_SUCCESS!=result, "could not create single vertex bump test")) return result;
  
  // Seal the mesh
  result=mw_func::make_mesh_watertight (input_set, facet_tol, verbose);
  if(gen::error(MB_SUCCESS!=result, "could not make the mesh watertight")) return result;

// Lastly Check to see if make_watertight fixed the model
  result=cw_func::check_mesh_for_watertightness( input_set, tol, sealed, test);
  if(gen::error(MB_SUCCESS!=result, "could not check model for watertightness")) return result;
  

    if(sealed)
  {
   std::cout << "PASS" << std::endl;
  }
  else
  {
   std::cout << "FAIL" << std::endl;
   exit(0);
  }

///////////////BEGIN 3RD TEST////////////////////////

  std::cout << "Test 3: vertex bump in x direction:" ;

  // Clear the mesh and reload original geometry for the next test
  result=reload_mesh( filename.c_str(), input_set);
  if (gen::error(MB_SUCCESS!=result, "could not reload the mesh" )) return result; 

  // retrieve the verticies again so the model can be broken
  result = MBI()->get_entities_by_dimension(input_set, dim, verts, false);
  if(gen::error(MB_SUCCESS!=result, " could not get vertices from the mesh")) return result;

  // "Break" the geometry
  result = single_vert_bump(verts, 0.1 , 0 , 0, root_name ); 
  if(gen::error(MB_SUCCESS!=result, "could not create single vertex bump test")) return result;
  
  // Seal the mesh
  result = mw_func::make_mesh_watertight (input_set, facet_tol, verbose);
  if(gen::error(MB_SUCCESS!=result, "could not make the mesh watertight")) return result;

// Lastly Check to see if make_watertight fixed the model
  result=cw_func::check_mesh_for_watertightness( input_set, tol, sealed, test);
  if(gen::error(MB_SUCCESS!=result, "could not check model for watertightness")) return result;
  

    if(sealed)
  {
   std::cout << "PASS" << std::endl;
  }
  else
  {
   std::cout << "FAIL" << std::endl;
   exit(0);
  }

///////////////BEGIN 4TH TEST////////////////////////

  std::cout << "Test 4: vertex bump in theta direction:" ;

  // Clear the mesh and reload original geometry for the next test
  result=reload_mesh( filename.c_str(), input_set);
  if (gen::error(MB_SUCCESS!=result, "could not reload the mesh" )) return result; 

  // retrieve the verticies again so the model can be broken
  result = MBI()->get_entities_by_dimension(input_set, dim, verts, false);
  if(gen::error(MB_SUCCESS!=result, " could not get vertices from the mesh")) return result;

  
  // "Break" the geometry
  result=theta_vert_bump(verts, 0, tol, root_name ); 
  if(gen::error(MB_SUCCESS!=result, "could not create single vertex bump test")) return result;
  
  // Seal the mesh
  result=mw_func::make_mesh_watertight (input_set, facet_tol, verbose);
  if(gen::error(MB_SUCCESS!=result, "could not make the mesh watertight")) return result;

// Lastly Check to see if make_watertight fixed the model
  result=cw_func::check_mesh_for_watertightness( input_set, tol, sealed, test);
  if(gen::error(MB_SUCCESS!=result, "could not check model for watertightness")) return result;
  

    if(sealed)
  {
   std::cout << "PASS" << std::endl;
  }
  else
  {
   std::cout << "FAIL" << std::endl;
   exit(0);
  }

///////////////BEGIN 5TH TEST////////////////////////

  std::cout << "Test 5: vertex bump in rand direction:" ;

  // Clear the mesh and reload original geometry for the next test
  result=reload_mesh( filename.c_str(), input_set);
  if (gen::error(MB_SUCCESS!=result, "could not reload the mesh" )) return result; 

  // retrieve the verticies again so the model can be broken
  result = MBI()->get_entities_by_dimension(input_set, dim, verts, false);
  if(gen::error(MB_SUCCESS!=result, " could not get vertices from the mesh")) return result;

  
  // "Break" the geometry
  result=rand_vert_bump(verts, tol, root_name, true ); 
  if(gen::error(MB_SUCCESS!=result, "could not create single vertex bump test")) return result;
  
  // Seal the mesh
  result=mw_func::make_mesh_watertight (input_set, facet_tol, verbose);
  if(gen::error(MB_SUCCESS!=result, "could not make the mesh watertight")) return result;

// Lastly Check to see if make_watertight fixed the model
  result=cw_func::check_mesh_for_watertightness( input_set, tol, sealed, test);
  if(gen::error(MB_SUCCESS!=result, "could not check model for watertightness")) return result;
  

    if(sealed)
  {
   std::cout << "PASS" << std::endl;
  }
  else
  {
   std::cout << "FAIL" << std::endl;
   exit(0);
  }

}

MBInterface* MBI() {
 static MBCore instance;
 return &instance;
}

