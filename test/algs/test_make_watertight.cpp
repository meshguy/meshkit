//
// Patrick Shriwise
// September 2013
// This program is designed to run a set of tests on the make_watertight algorithm.
// This will be a stand-alone piece of code that uses MOAB to open, modify
// (break), and re-seal geometries/ 
// input: cyl.h5m file (found in ../make_watertight/test/)
// output: pass/fail for each of the tests


#include <iostream>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <vector>
#include <set>
#include <algorithm>
#include "MBCore.hpp"
#include "moab::TagConventions.hpp"
#include "moab/Range.hpp"
#include "moab/Skinner.hpp"

#include "meshkit/mw_func.hpp"
#include "meshkit/cw_func.hpp"
#include "meshkit/gen.hpp"
#include "meshkit/arc.hpp"
#include "meshkit/zip.hpp"

#include "TestUtil.hpp"

moab::ErrorCode move_vert( moab::EntityHandle vertex, double dx, double dy, double dz, bool verbose = false );
moab::ErrorCode rand_vert_move( moab::EntityHandle vertex, double tol, bool verbose = false);
moab::ErrorCode move_vert_theta( moab::EntityHandle vertex, double tolerance, bool verbose = false);
moab::ErrorCode move_vert_R ( moab::EntityHandle vertex, double tol, bool verbose = false ) ;
moab::ErrorCode rand_adj_pair( moab::Range verts, moab::EntityHandle &vert1, moab::EntityHandle &vert2);


/// struct to hold coordinates of skin edge, it's surface id, and a matched flag
struct coords_and_id {
  double x1;
  double y1;
  double z1;
  double x2;
  double y2;
  double z2;
  int  surf_id;
  bool matched;
  moab::EntityHandle vert1;
  moab::EntityHandle vert2;
};

/// takes an entity handle vertex, gets the original coordinates, changes and resets the vertex coordinates in the mesh
moab::ErrorCode move_vert( moab::EntityHandle vertex, double dx, double dy, double dz, bool verbose) {

 moab::ErrorCode result;
 
 //get coordinates from the mesh
 double coords[3];
 result= MBI()-> get_coords( &vertex , 1 , coords );
 if(gen::error(moab::MB_SUCCESS!=result, "could not get the vertex coordinates")) return result; 

 if(verbose)
 {
  //original coordinates
  std::cout << std::endl << "Original Coordinates" << std::endl;
  std::cout << "x = " << coords[0] << std::endl;
  std::cout << "y = " << coords[1] << std::endl;
  std::cout << "z = " << coords[2] << std::endl;
 }
  coords[0]+=dx;
  coords[1]+=dy;
  coords[2]+=dz;
  
 if(verbose)
 {
  //altered coordinates
  std::cout << std::endl << "Modified Coordinates" << std::endl;
  std::cout << "x = " << coords[0] << std::endl;
  std::cout << "y = " << coords[1] << std::endl;
  std::cout << "z = " << coords[2] << std::endl;
 }
 
  //write new coordinates to the mesh
  result=MBI()-> set_coords( &vertex, 1, coords);
  if (gen::error(moab::MB_SUCCESS!=result, "could not set the vertex coordinates")) return result;

return moab::MB_SUCCESS;
}

/// takes an entity handle vertex, gets the original coordinates, changes and resets the vertex coordinates in the mesh
/// setup to move the vert no further than the faceting tolerance
moab::ErrorCode rand_vert_move( moab::EntityHandle vertex, double tol, bool verbose) {

 moab::ErrorCode result;
 
 //get coordinates from the mesh
 double coords[3];
 result= MBI()-> get_coords( &vertex , 1 , coords );
 if(gen::error(moab::MB_SUCCESS!=result, "could not get the vertex coordinates")) return result; 

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
 }
  coords[0]+=dx;
  coords[1]+=dy;
  coords[2]+=dz;

 if(verbose)
 {
  //altered coordinates
  std::cout << std::endl << "Modified Coordinates" << std::endl;
  std::cout << "x = " << coords[0] << std::endl;
  std::cout << "y = " << coords[1] << std::endl;
  std::cout << "z = " << coords[2] << std::endl;
 }
 
  //write new coordinates to the mesh
  result=MBI()-> set_coords( &vertex, 1, coords);
  if (gen::error(moab::MB_SUCCESS!=result, "could not set the vertex coordinates")) return result;

return moab::MB_SUCCESS;
}

/// moves a vertex along the rim of the cylinder in the theta direction a distance equal to the faceting_tolerance
moab::ErrorCode move_vert_theta( moab::EntityHandle vertex, double tolerance, bool verbose) {

 moab::ErrorCode result; 
  
  //get vertex coordinates
  double coords[3];
  result= MBI()-> get_coords( &vertex , 1 , coords );
 
  // determine radius
  double radius=sqrt(pow(coords[0],2)+pow(coords[1],2));
  //std::cout << "radius = " << radius << std::endl;

  // get the current theta value
  // need both because of the oddness/evenness of the inverse functions
  double theta_x=acos(coords[0]/radius);
  double theta_y=asin(coords[1]/radius);

  //std::cout << "theta = " << theta << std::endl;
  
  // set the vertex bump distance
  double dtheta = tolerance/(radius);
  //std::cout << "dtheta = " << dtheta << std::endl;

  if(verbose){
  //original coordinates
  std::cout << std::endl << "Original Coordinates" << std::endl;
  std::cout << "x = " << coords[0] << std::endl;
  std::cout << "y = " << coords[1] << std::endl;
  std::cout << "z = " << coords[2] << std::endl;
 }
 // create new x and y values
  coords[0]=radius*cos(theta_x+dtheta); 
  //std::cout << "new x value = " << coords[0] << std::endl;
  coords[1]=radius*sin(theta_y+dtheta);
  //std::cout << "new y value = " << coords[1] << std::endl;
  if(verbose){
  //altered coordinates
  std::cout << std::endl << "Modified Coordinates" << std::endl;
  std::cout << "x = " << coords[0] << std::endl;
  std::cout << "y = " << coords[1] << std::endl;
  std::cout << "z = " << coords[2] << std::endl;
 }
  //set new vertex coordinates  
  result=MBI()-> set_coords( &vertex, 1, coords);
  if (gen::error(moab::MB_SUCCESS!=result, "could not set the vertex coordinates")) return result;

  return moab::MB_SUCCESS;
}

/// moves the vertex in R some distance less than tol
moab::ErrorCode move_vert_R( moab::EntityHandle vertex, double tol, bool verbose ) {

  moab::ErrorCode result; 
  
  //get vertex coordinates
  double coords[3];
  result= MBI()-> get_coords( &vertex , 1 , coords );
  if(gen::error(moab::MB_SUCCESS!=result, "could not get the vertex coordinates"));

  if(verbose){
  //original coordinates
  std::cout << "Vertex ID: " << gen::geom_id_by_handle(vertex) << std::endl;
  std::cout << std::endl << "Original Coordinates" << std::endl;
  std::cout << "x = " << coords[0] << std::endl;
  std::cout << "y = " << coords[1] << std::endl;
  std::cout << "z = " << coords[2] << std::endl;
  }
  double radius = sqrt(pow(coords[0],2)+pow(coords[1],2));
  //get unit vector in x-y plane
  coords[0]/=radius;
  coords[1]/=radius;
  
  //alter radius to new value of radius+tol
  radius-=tol;
  coords[0]*=radius;
  coords[1]*=radius;
 
  if(verbose){
  //altered coordinates
  std::cout << std::endl << "Modified Coordinates" << std::endl;
  std::cout << "x = " << coords[0] << std::endl;
  std::cout << "y = " << coords[1] << std::endl;
  std::cout << "z = " << coords[2] << std::endl;
  }
  
  //set new vertex coordinates  
  result=MBI()-> set_coords( &vertex, 1, coords);
  if (gen::error(moab::MB_SUCCESS!=result, "could not set the vertex coordinates")) return result;

return moab::MB_SUCCESS;
}


/// appends "_mod" to the original file name and writes to a new .h5m 
moab::ErrorCode write_mod_file( std::string filename ) {

 moab::ErrorCode result;
 std::string output_filename = filename + "_mod.h5m";
 //write file
 result=MBI()->write_mesh(output_filename.c_str());
 if(gen::error(moab::MB_SUCCESS!=result,"could not write the mesh to output_filename")) return result; 

  return moab::MB_SUCCESS;
}


// used to clear all mesh data and reload the file as original
moab::ErrorCode reload_mesh(const char* filename,  moab::EntityHandle &meshset, bool debug = false) {

  moab::ErrorCode result;

  // delete meshset
  result= MBI() -> delete_mesh();
  if(gen::error(moab::MB_SUCCESS!=result, "could not delete the mesh set")) return result;



  // re-initialize meshset
  result= MBI() -> create_meshset(MESHSET_SET, meshset);
  if(gen::error(moab::MB_SUCCESS!=result, "could not create the new meshset")) return result; 
 
  //reload the file
  result = MBI() -> load_file(filename, &meshset);
  if(gen::error(moab::MB_SUCCESS!=result, "could not re-load the file into the mesh set")) return result; 

  //check that something was loaded into the meshset
  int num_meshsets;   

  result= MBI() -> num_contained_meshsets ( meshset, &num_meshsets);
  if(gen::error(moab::MB_SUCCESS!=result, "could not get the number of meshsets")) return result; 
  
  if(debug) std::cout << "number of mesh sets after reload = " << num_meshsets << std::endl;

  if(num_meshsets==0) return MB_FAILURE;

  return result;
}


/// bumps the last vertex in the model by the x,y,z values given to the problem 
moab::ErrorCode single_vert_bump( moab::Range verts, double bump_dist_x, double bump_dist_y, double bump_dist_z, bool verbose = false ) {
 
  moab::ErrorCode result; 
  moab::EntityHandle vertex=verts.back();
  //move the desired vertex by the allotted distance
  result = move_vert( vertex, bump_dist_x, bump_dist_y, bump_dist_z );
  if(gen::error(moab::MB_SUCCESS!=result, "could not move single vert")) return result;

  return moab::MB_SUCCESS;
}

/// bumps the last vertex in the cylinder model in the R direction
moab::ErrorCode single_vert_bump_R( moab::Range verts, double facet_tol, bool verbose = false ) {
 
  moab::ErrorCode result; 
  moab::EntityHandle vertex=verts.back();
  //move the desired vertex by the allotted distance
  result = move_vert_R( vertex, facet_tol, verbose );
  if(gen::error(moab::MB_SUCCESS!=result, "could not move single vert")) return result;

  return moab::MB_SUCCESS;
}

/// moves the last two verticies in the model the same distance in x, y, and z 
moab::ErrorCode locked_pair_bump( moab::Range verts, double bump_dist_x, double bump_dist_y, double bump_dist_z,  std::string root_name, bool verbose = false ) {
 
  moab::ErrorCode result; 

  moab::EntityHandle vertex1=verts.back();
  moab::EntityHandle vertex2=(verts.back()-1);
 
  //move the desired verticies by the allotted distance(s)
  result = move_vert( vertex1, bump_dist_x, bump_dist_y, bump_dist_z, verbose);
  if(gen::error(moab::MB_SUCCESS!=result, "could not move vertex1")) return result;
  result = move_vert( vertex2, bump_dist_x, bump_dist_y, bump_dist_z, verbose);
  if(gen::error(moab::MB_SUCCESS!=result, "could not move vertex2")) return result;


  return moab::MB_SUCCESS;
}
/// selects a random pair of verticies from verts and moves them in random directions some distance less than the faceting tolerance
moab::ErrorCode rand_locked_pair_bump_rand( moab::Range verts, double facet_tol , std::string root_name, bool verbose = false ) {
 
  moab::ErrorCode result; 

  //select random verticies from verts
  moab::EntityHandle vertex1;
  moab::EntityHandle vertex2;

  result = rand_adj_pair( verts, vertex1, vertex2);
  
  //move the desired verticies by the allotted distance(s)
  result = rand_vert_move( vertex1, facet_tol , verbose);
  if(gen::error(moab::MB_SUCCESS!=result, "could not move vertex1")) return result;
  result = rand_vert_move( vertex2, facet_tol, verbose);
  if(gen::error(moab::MB_SUCCESS!=result, "could not move vertex2")) return result;

  return moab::MB_SUCCESS;
}

/// selects a random pair of verticies and moves them along theta a distance less than the faceting tolerance
moab::ErrorCode rand_locked_pair_bump_theta( moab::Range verts, double facet_tol , std::string root_name, bool verbose = false ) {
 
  moab::ErrorCode result; 

  //select random verticies from verts
  moab::EntityHandle vertex1;
  moab::EntityHandle vertex2;

  result = rand_adj_pair(verts, vertex1, vertex2);
  
  //move the desired verticies by the allotted distance(s)
  result = move_vert_theta( vertex1, facet_tol , verbose);
  if(gen::error(moab::MB_SUCCESS!=result, "could not move vertex1")) return result;
  result = move_vert_theta( vertex2, facet_tol, verbose);
  if(gen::error(moab::MB_SUCCESS!=result, "could not move vertex2")) return result;

  return moab::MB_SUCCESS;
}
/// selects a random pair of adjacent verticies and bumps them in x, y, and z 
moab::ErrorCode rand_locked_pair_bump( moab::Range verts, double bump_dist_x, double bump_dist_y, double bump_dist_z,  std::string root_name, bool verbose = false ) {
 
  moab::ErrorCode result; 

  //select random verticies from verts
  moab::EntityHandle vertex1;
  moab::EntityHandle vertex2;

  result = rand_adj_pair( verts, vertex1, vertex2);
  
  //move the desired verticies by the allotted distance(s)
  result = move_vert( vertex1, bump_dist_x, bump_dist_y, bump_dist_z, verbose);
  if(gen::error(moab::MB_SUCCESS!=result, "could not move vertex1")) return result;
  result = move_vert( vertex2, bump_dist_x, bump_dist_y, bump_dist_z, verbose);
  if(gen::error(moab::MB_SUCCESS!=result, "could not move vertex2")) return result;

  return moab::MB_SUCCESS;
}
/// moves the last two verticies in the model in the same direction some distance less than the faceting tolerance
moab::ErrorCode locked_pair_bump_rand( moab::Range verts, double facet_tol,  std::string root_name, bool verbose = false ) {
 
  moab::ErrorCode result; 

  moab::EntityHandle vertex1=verts.back();
  moab::EntityHandle vertex2=(verts.back()-1);
  
  //move the desired verticies by the allotted distance(s)
  result = rand_vert_move( vertex1, facet_tol, verbose);
  if(gen::error(moab::MB_SUCCESS!=result, "could not move vertex1")) return result;
  result = rand_vert_move( vertex2, facet_tol, verbose);
  if(gen::error(moab::MB_SUCCESS!=result, "could not move vertex2")) return result;

  return moab::MB_SUCCESS;
}


/// moves the last vertex in the model in a random direction by a distance less than the faceting tolerance
moab::ErrorCode rand_vert_bump( moab::Range verts, double facet_tol, std::string root_name, bool verbose = false ) {
 
  moab::ErrorCode result; 
  moab::EntityHandle vertex = verts.back();
  //move the desired vertex by the allotted distance
  result = rand_vert_move( vertex, facet_tol, verbose);
  if(gen::error(moab::MB_SUCCESS!=result, "could not move single vert")) return result;


  return moab::MB_SUCCESS;
}

// FOR CYLINDER TESTING ONLY 
/// moves the last vertex in the model along the curve of the cylinder some distance bump distance theta
moab::ErrorCode theta_vert_bump( moab::Range verts, double bump_dist_theta, double tolerance, std::string root_name, bool verbose = false ) {
 
  moab::ErrorCode result; 
  
  //get vertex coordinates
  double coords[3];
  moab::EntityHandle vertex=verts.back();
  result= MBI()-> get_coords( &vertex , 1 , coords );
 
  // determine radius
  double radius=sqrt(pow(coords[0],2)+pow(coords[1],2));
  //std::cout << "radius = " << radius << std::endl;

  // get the current theta value
  double theta=asin(coords[1]/radius);
  //std::cout << "theta = " << theta << std::endl;
  
  // set the vertex bump distance
  double dtheta = 0.5*tolerance/(radius);
  //std::cout << "dtheta = " << dtheta << std::endl;

 
 if(verbose)
 {
  //original coordinates
  std::cout << std::endl << "Original Coordinates" << std::endl;
  std::cout << "x = " << coords[0] << std::endl;
  std::cout << "y = " << coords[1] << std::endl;
  std::cout << "z = " << coords[2] << std::endl;
  }
 // create new x and y values
  coords[0]=radius*cos(theta+dtheta); 
  //std::cout << "new x value = " << coords[0] << std::endl;
  coords[1]=radius*sin(theta+dtheta);
  //std::cout << "new y value = " << coords[1] << std::endl;
  if(verbose)
  {
  //altered coordinates
  std::cout << std::endl << "Modified Coordinates" << std::endl;
  std::cout << "x = " << coords[0] << std::endl;
  std::cout << "y = " << coords[1] << std::endl;
  std::cout << "z = " << coords[2] << std::endl;
  }
  //write new coordinates to the mesh
  // might not be necesarry any longer as we move to doing tests on a moab-instance basis
  result=MBI()-> set_coords( &vertex, 1, coords);
  if (gen::error(moab::MB_SUCCESS!=result, "could not set the vertex coordinates")) return result;
  // alter output filename
 
  return moab::MB_SUCCESS;
}

/// moves two adjacent vertices along theta a distance equal to the faceting tolerance
moab::ErrorCode locked_pair_move_theta( moab::Range verts, double tolerance, std::string root_name,  bool verbose = false) {

  moab::ErrorCode result;

  //get vertex coordinates
  moab::EntityHandle vertex1=verts.back();
  moab::EntityHandle vertex2=verts.back()-1;
  
  result= move_vert_theta( vertex1, tolerance, verbose);
  if(gen::error(moab::MB_SUCCESS!=result,"could not move vertex1 along theta")) return result; 
  result= move_vert_theta( vertex2, tolerance, verbose);
  if(gen::error(moab::MB_SUCCESS!=result,"could not move vertex1 along theta")) return result; 

  return moab::MB_SUCCESS;
}

/// moves the third to last and the last vertices in the model the same distance in x, y, and z
moab::ErrorCode adjplone_locked_pair_bump( moab::Range verts, double bump_dist_x, double bump_dist_y, double bump_dist_z,  std::string root_name, bool verbose = false ) {
 
  moab::ErrorCode result; 

  moab::EntityHandle vertex1=verts.back();
  moab::EntityHandle vertex2=(verts.back()-2);
 
  //move the desired verticies by the allotted distance(s)
  result = move_vert( vertex1, bump_dist_x, bump_dist_y, bump_dist_z, verbose);
  if(gen::error(moab::MB_SUCCESS!=result, "could not move vertex1")) return result;
  result = move_vert( vertex2, bump_dist_x, bump_dist_y, bump_dist_z, verbose);
  if(gen::error(moab::MB_SUCCESS!=result, "could not move vertex2")) return result;


  return moab::MB_SUCCESS;
}

/// moves the third to last and the last verticies in the model in theta the same distance along theta equal to the faceting tolerance
moab::ErrorCode adjplone_locked_pair_bump_theta( moab::Range verts, double facet_tol , std::string root_name, bool verbose = false ) {
 
  moab::ErrorCode result; 

  moab::EntityHandle vertex1=verts.back();
  moab::EntityHandle vertex2=(verts.back()-2);
 
  //move the desired verticies by the allotted distance(s)
  result = move_vert_theta( vertex1, facet_tol, verbose);
  if(gen::error(moab::MB_SUCCESS!=result, "could not move vertex1")) return result;
  result = move_vert_theta( vertex2, facet_tol, verbose);
  if(gen::error(moab::MB_SUCCESS!=result, "could not move vertex2")) return result;


  return moab::MB_SUCCESS;
}

/// moves the third to last and the last verticies in the model in rand directions some distance less than the facet_tolerance
moab::ErrorCode adjplone_locked_pair_bump_rand( moab::Range verts, double facet_tol , std::string root_name, bool verbose = false ) {
 
  moab::ErrorCode result; 

  moab::EntityHandle vertex1=verts.back();
  moab::EntityHandle vertex2=(verts.back()-2);
 
  //move the desired verticies by the allotted distance(s)
  result = rand_vert_move( vertex1, facet_tol, verbose);
  if(gen::error(moab::MB_SUCCESS!=result, "could not move vertex1")) return result;
  result = rand_vert_move( vertex2, facet_tol, verbose);
  if(gen::error(moab::MB_SUCCESS!=result, "could not move vertex2")) return result;


  return moab::MB_SUCCESS;
}


/// selects a random pair of adjacent verticies and bumps them the same distance in x, y, and z
moab::ErrorCode nonadj_locked_pair_bump( moab::Range verts, double bump_dist_x, double bump_dist_y, double bump_dist_z,  std::string root_name, bool verbose = false ) {
 
  moab::ErrorCode result; 

  //select random verticies from verts
  int number_of_verts = verts.size();
  double num = rand()/static_cast<double>(RAND_MAX);
  //std::cout << "max index = " << (number_of_verts-2)<< std::endl;
  int index = static_cast<int>(num*((number_of_verts-2)));
  //std::cout << "index = " << index << std::endl;

  moab::EntityHandle vertex1=verts.back();
  moab::EntityHandle vertex2=(verts.back()-index);
  
  //move the desired verticies by the allotted distance(s)
  result = move_vert( vertex1, bump_dist_x, bump_dist_y, bump_dist_z, verbose);
  if(gen::error(moab::MB_SUCCESS!=result, "could not move vertex1")) return result;
  result = move_vert( vertex2, bump_dist_x, bump_dist_y, bump_dist_z, verbose);
  if(gen::error(moab::MB_SUCCESS!=result, "could not move vertex2")) return result;

  return moab::MB_SUCCESS;
}

/// selects a random pair of adjacent verticies and bumps them along the theta direction a distance equal to the faceting tolerance
moab::ErrorCode nonadj_locked_pair_bump_theta( moab::Range verts, double facet_tol , std::string root_name, bool verbose = false ) {
 
  moab::ErrorCode result; 

  //select random verticies from verts
  int number_of_verts = verts.size();
  double num = rand()/static_cast<double>(RAND_MAX);
  //std::cout << "max index = " << (number_of_verts-2)<< std::endl;
  int index = static_cast<int>(num*((number_of_verts-2)));
  //std::cout << "index = " << index << std::endl;

  moab::EntityHandle vertex1=verts.back();
  moab::EntityHandle vertex2=(verts.back()-index);
  
  //move the desired verticies by the allotted distance(s)
  result = move_vert_theta( vertex1, facet_tol, verbose);
  if(gen::error(moab::MB_SUCCESS!=result, "could not move vertex1")) return result;
  result = move_vert_theta( vertex2, facet_tol , verbose);
  if(gen::error(moab::MB_SUCCESS!=result, "could not move vertex2")) return result;

  return moab::MB_SUCCESS;
}

/// selects a random pair of adjacent verticies and bumps them along the theta direction a distance equal to the faceting tolerance
moab::ErrorCode nonadj_locked_pair_bump_rand( moab::Range verts, double facet_tol , std::string root_name, bool verbose = false ) {
 
  moab::ErrorCode result; 

  //select random verticies from verts
  int number_of_verts = verts.size();
  double num = rand()/static_cast<double>(RAND_MAX);
  //std::cout << "max index = " << (number_of_verts-2)<< std::endl;
  int index = static_cast<int>(num*((number_of_verts-2)));
  //std::cout << "index = " << index << std::endl;

  moab::EntityHandle vertex1=verts.back();
  moab::EntityHandle vertex2=(verts.back()-index);
  
  //move the desired verticies by the allotted distance(s)
  result = rand_vert_move( vertex1, facet_tol, verbose);
  if(gen::error(moab::MB_SUCCESS!=result, "could not move vertex1")) return result;
  result = rand_vert_move( vertex2, facet_tol , verbose);
  if(gen::error(moab::MB_SUCCESS!=result, "could not move vertex2")) return result;

  return moab::MB_SUCCESS;
}

moab::ErrorCode locked_pair_bump_R( moab::Range verts, double facet_tol,  std::string root_name, bool verbose = false ) {
 
  moab::ErrorCode result; 

  moab::EntityHandle vertex1=verts.back();
  moab::EntityHandle vertex2=(verts.back()-1);
 
  //move the desired verticies by the allotted distance(s)
  result = move_vert_R( vertex1, facet_tol, verbose);
  if(gen::error(moab::MB_SUCCESS!=result, "could not move vertex1")) return result;
  result = move_vert_R( vertex2, facet_tol, verbose);
  if(gen::error(moab::MB_SUCCESS!=result, "could not move vertex2")) return result;


  return moab::MB_SUCCESS;
}

/// selects random verticies from verts and moves them in R a distance equal to the faceting tolerance
moab::ErrorCode rand_locked_pair_bump_R( moab::Range verts, double facet_tol , std::string root_name, bool verbose = false ) {
 
  moab::ErrorCode result; 

  //select random verticies from verts
  moab::EntityHandle vertex1;
  moab::EntityHandle vertex2;

  result = rand_adj_pair( verts, vertex1, vertex2);
  
  //move the desired verticies by the allotted distance(s)
  result = move_vert_R( vertex1, facet_tol , verbose);
  if(gen::error(moab::MB_SUCCESS!=result, "could not move vertex1")) return result;
  result = move_vert_R( vertex2, facet_tol, verbose);
  if(gen::error(moab::MB_SUCCESS!=result, "could not move vertex2")) return result;

  return moab::MB_SUCCESS;
}

/// selects a the last vertex and third to last vertex in the model and moves them in R a distance equal to the faceting tolerance
moab::ErrorCode adjplone_locked_pair_bump_R( moab::Range verts, double facet_tol , std::string root_name, bool verbose = false ) {
 
  moab::ErrorCode result; 

  moab::EntityHandle vertex1=verts.back();
  moab::EntityHandle vertex2=(verts.back()-2);
 
  //move the desired verticies by the allotted distance(s)
  result = move_vert_R( vertex1, facet_tol, verbose);
  if(gen::error(moab::MB_SUCCESS!=result, "could not move vertex1")) return result;
  result = move_vert_R( vertex2, facet_tol, verbose);
  if(gen::error(moab::MB_SUCCESS!=result, "could not move vertex2")) return result;


  return moab::MB_SUCCESS;
}

/// moves the last vertex in the model and a randomly selected, non-adjacent vertex and moves them both in R a distance equal to the faceting tolerance
moab::ErrorCode nonadj_locked_pair_bump_R( moab::Range verts, double facet_tol , std::string root_name, bool verbose = false ) {
 
  moab::ErrorCode result; 

  //select random verticies from verts
  int number_of_verts = verts.size();
  double num = rand()/static_cast<double>(RAND_MAX);
  //std::cout << "max index = " << (number_of_verts-2)<< std::endl;
  int index = static_cast<int>(num*((number_of_verts-2)));
  //std::cout << "index = " << index << std::endl;

  moab::EntityHandle vertex1=verts.back();
  moab::EntityHandle vertex2=(verts.back()-index);
  
  //move the desired verticies by the allotted distance(s)
  result = move_vert_R( vertex1, facet_tol, verbose);
  if(gen::error(moab::MB_SUCCESS!=result, "could not move vertex1")) return result;
  result = move_vert_R( vertex2, facet_tol , verbose);
  if(gen::error(moab::MB_SUCCESS!=result, "could not move vertex2")) return result;

  return moab::MB_SUCCESS;
}
 
moab::ErrorCode rand_adj_pair( moab::Range verts, moab::EntityHandle &vert1, moab::EntityHandle &vert2)
{


  //select random verticies from verts
  int number_of_verts = verts.size();
  double num = rand()/static_cast<double>(RAND_MAX);
 
  int index = static_cast<int>(num*(number_of_verts-2));


  vert1=(verts.front()+index);
  vert2=(verts.front()+(index+1));


  return moab::MB_SUCCESS;
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

  //clock_t start_time;
  //start_time = clock();
  //save the mesh to a new filename
  std::string input_name="cyl.h5m";
  std::string root_name=input_name;
  int len = root_name.length();
  root_name.erase(len-4);

  // load file and get tolerance from input argument
  moab::ErrorCode result;
  std::string filename = TestDir + "/" + input_name; //set filename
  moab::EntityHandle input_set;
  result = MBI()->create_meshset( MESHSET_SET, input_set ); //create handle to meshset
  if(moab::MB_SUCCESS != result) 
    {
      return result;
    }

result = MBI()->load_file( filename.c_str(), &input_set ); //load the file into the meshset
  if(moab::MB_SUCCESS != result) 
    {
      // failed to load the file
      std::cout << "could not load file" << std::endl;
      return result;
    }

  //// get faceting tolerance ////
  double facet_tolerance;

  moab::Tag faceting_tol_tag;
  //get faceting tolerance handle from file
  result = MBI()->tag_get_handle( "FACETING_TOL", 1, MB_TYPE_DOUBLE,
        faceting_tol_tag , moab::MB_TAG_SPARSE|moab::MB_TAG_CREAT );
  if(gen::error(moab::MB_SUCCESS!=result, "could not get the faceting tag handle")) return result;
  
  //get the faceting tolerance of any entity
  moab::Range file_set;
  result = MBI()->get_entities_by_type_and_tag( 0, MBENTITYSET, 
                        &faceting_tol_tag, NULL, 1, file_set );

  //get facetint tolerance value
  result = MBI()->tag_get_data( faceting_tol_tag, &file_set.front(), 1, &facet_tolerance );
  if(gen::error(moab::MB_SUCCESS!=result, "could not get the faceting tolerance")) return result; 

  // create tags on geometry
  moab::Tag geom_tag, id_tag;
  result = MBI()->tag_get_handle( GEOM_DIMENSION_TAG_NAME, 1, 
                            MB_TYPE_INTEGER, geom_tag, moab::MB_TAG_DENSE|moab::MB_TAG_CREAT );
  if(gen::error(moab::MB_SUCCESS != result, "could not get GEOM_DIMENSION_TAG_NAME handle")) return result; 

  result = MBI()->tag_get_handle( GLOBAL_ID_TAG_NAME, 1, 
                            MB_TYPE_INTEGER, id_tag, moab::MB_TAG_DENSE|moab::MB_TAG_CREAT );
  if(gen::error(moab::MB_SUCCESS != result, "could not get GLOBAL_ID_TAG_NAME handle")) return result;
  
  // get surface and volume sets
  moab::Range surf_sets, vol_sets; // moab::Range of set of surfaces and volumes
  // surface sets
  int dim = 2;
  void* input_dim[] = {&dim};
  result = MBI()->get_entities_by_type_and_tag( input_set, MBENTITYSET, &geom_tag, 
                                                input_dim, 1, surf_sets);
  if(moab::MB_SUCCESS != result) 
    {
      return result;
    }

  // volume sets
  dim = 3;
  result = MBI()->get_entities_by_type_and_tag( input_set, MBENTITYSET, &geom_tag, 
                                                input_dim, 1, vol_sets);
  if(moab::MB_SUCCESS != result)
    {
      return result;
    }
    
  //vertex sets
  dim= 0;
  moab::Range verts;
  result = MBI()->get_entities_by_dimension(input_set, dim, verts, false);
  if(gen::error(moab::MB_SUCCESS!=result, " could not get vertices from the mesh")) return result;
  
  if(gen::error(moab::MB_SUCCESS!=result, "could not get vertex coordinates")) return result; 
  
if(verbose)
{
  std::cout<< "number of verticies= " << verts.size() << std::endl;  
  std::cout<< "number of surfaces= " << surf_sets.size() << std::endl;
  std::cout<< "number of volumes= "  << vol_sets.size() << std::endl;
}

//initialize booleans to pass to make_mesh_watertight

  //bool check_topology;
  bool test, sealed;
       //check_topology=false;
       test=true;
  
// initialize boolean for each set of tests
  bool test_set_result=true;

///////////Single Verticie Movement Tests////////////////////
  std::cout << "SINGLE VERTEXT MOVEMENT TESTS" << std::endl;

///////////////BEGIN 1st TEST////////////////////////

 std::cout << "Test 1: vertex bump in x direction:";

  //perform single vertex bump and test
  result=single_vert_bump(verts, 0.9*facet_tolerance , 0 , 0 ); 
  if(gen::error(moab::MB_SUCCESS!=result, "could not create single vertex bump test")) return result;
  
  // seal the model using make_watertight 
  result=mw_func::make_mesh_watertight (input_set, facet_tolerance, verbose);
  if(gen::error(moab::MB_SUCCESS!=result, "could not make the mesh watertight")) return result;
  
  // Lastly Check to see if make_watertight fixed the model
  result=cw_func::check_mesh_for_watertightness( input_set, facet_tolerance, sealed, test);
  if(gen::error(moab::MB_SUCCESS!=result, "could not check model for watertightness")) return result;
  
  if(sealed)
  {
   std::cout << "PASS" << std::endl;
  }
  else
  {
   std::cout << "FAIL" << std::endl;
   test_set_result=false;
  }

///////////////BEGIN 2ND TEST////////////////////////

  std::cout << "Test 2: vertex bump in y direction:";

  // Clear the mesh and reload original geometry for the next test
  result=reload_mesh( filename.c_str(), input_set);
  if (gen::error(moab::MB_SUCCESS!=result, "could not reload the mesh" )) return result; 

  // retrieve the verticies again so the model can be broken
  result = MBI()->get_entities_by_dimension(input_set, dim, verts, false);
  if(gen::error(moab::MB_SUCCESS!=result, " could not get vertices from the mesh")) return result;

  // "Break" the geometry
  result=single_vert_bump(verts, 0 , 0.9*facet_tolerance , 0); 
  if(gen::error(moab::MB_SUCCESS!=result, "could not create single vertex bump test")) return result;
  
  // Seal the mesh
  result=mw_func::make_mesh_watertight (input_set, facet_tolerance, verbose);
  if(gen::error(moab::MB_SUCCESS!=result, "could not make the mesh watertight")) return result;

// Lastly Check to see if make_watertight fixed the model
  result=cw_func::check_mesh_for_watertightness( input_set, facet_tolerance, sealed, test);
  if(gen::error(moab::MB_SUCCESS!=result, "could not check model for watertightness")) return result;
  

    if(sealed)
  {
   std::cout << "PASS" << std::endl;
  }
  else
  {
   std::cout << "FAIL" << std::endl;
   test_set_result=false;
  }


///////////////BEGIN 3RD TEST////////////////////////

  std::cout << "Test 3: vertex bump in z direction:" ;

  // Clear the mesh and reload original geometry for the next test
  result=reload_mesh( filename.c_str(), input_set);
  if (gen::error(moab::MB_SUCCESS!=result, "could not reload the mesh" )) return result; 

  // retrieve the verticies again so the model can be broken
  result = MBI()->get_entities_by_dimension(input_set, dim, verts, false);
  if(gen::error(moab::MB_SUCCESS!=result, " could not get vertices from the mesh")) return result;

  // "Break" the geometry
  result = single_vert_bump(verts, 0 , 0 , 0.9*facet_tolerance ); 
  if(gen::error(moab::MB_SUCCESS!=result, "could not create single vertex bump test")) return result;
  
  // Seal the mesh
  result = mw_func::make_mesh_watertight (input_set, facet_tolerance, verbose);
  if(gen::error(moab::MB_SUCCESS!=result, "could not make the mesh watertight")) return result;

// Lastly Check to see if make_watertight fixed the model
  result=cw_func::check_mesh_for_watertightness( input_set, facet_tolerance, sealed, test);
  if(gen::error(moab::MB_SUCCESS!=result, "could not check model for watertightness")) return result;
  

    if(sealed)
  {
   std::cout << "PASS" << std::endl;
  }
  else
  {
   std::cout << "FAIL" << std::endl;
   test_set_result=false;
  }


///////////////BEGIN 4th TEST ////////////////////////

  std::cout << "Test 4: vertex bump in R direction:" ;

  // Clear the mesh and reload original geometry for the next test
  result=reload_mesh( filename.c_str(), input_set);
  if (gen::error(moab::MB_SUCCESS!=result, "could not reload the mesh" )) return result; 

  // retrieve the verticies again so the model can be broken
  result = MBI()->get_entities_by_dimension(input_set, dim, verts, false);
  if(gen::error(moab::MB_SUCCESS!=result, " could not get vertices from the mesh")) return result;

    // "Break" the geometry
  result=single_vert_bump_R(verts, facet_tolerance ); 
  if(gen::error(moab::MB_SUCCESS!=result, "could not create single vertex bump test")) return result;
  

  // Seal the mesh
  result=mw_func::make_mesh_watertight (input_set, facet_tolerance, verbose);
  if(gen::error(moab::MB_SUCCESS!=result, "could not make the mesh watertight")) return result;

  // Lastly Check to see if make_watertight fixed the model
  result=cw_func::check_mesh_for_watertightness( input_set, facet_tolerance, sealed, test);
  if(gen::error(moab::MB_SUCCESS!=result, "could not check model for watertightness")) return result;
  

    if(sealed)
  {
   std::cout << "PASS" << std::endl;
  }
  else
  {
   std::cout << "FAIL" << std::endl;
   test_set_result=false;
  }
///////////////BEGIN 5TH TEST////////////////////////

  std::cout << "Test 5: vertex bump in theta direction:" ;

  // Clear the mesh and reload original geometry for the next test
  result=reload_mesh( filename.c_str(), input_set);
  if (gen::error(moab::MB_SUCCESS!=result, "could not reload the mesh" )) return result; 

  // retrieve the verticies again so the model can be broken
  result = MBI()->get_entities_by_dimension(input_set, dim, verts, false);
  if(gen::error(moab::MB_SUCCESS!=result, " could not get vertices from the mesh")) return result;

  
  // "Break" the geometry
  result=theta_vert_bump(verts, 0, facet_tolerance, root_name ); 
  if(gen::error(moab::MB_SUCCESS!=result, "could not create single vertex bump test")) return result;
  
  // Seal the mesh
  result=mw_func::make_mesh_watertight (input_set, facet_tolerance, verbose);
  if(gen::error(moab::MB_SUCCESS!=result, "could not make the mesh watertight")) return result;

  // Lastly Check to see if make_watertight fixed the model
  result=cw_func::check_mesh_for_watertightness( input_set, facet_tolerance, sealed, test);
  if(gen::error(moab::MB_SUCCESS!=result, "could not check model for watertightness")) return result;
  

    if(sealed)
  {
   std::cout << "PASS" << std::endl;
  }
  else
  {
   std::cout << "FAIL" << std::endl;
   test_set_result=false;
  }

///////////////BEGIN 6TH TEST////////////////////////

  std::cout << "Test 6: vertex bump in rand direction:" ;

  // Clear the mesh and reload original geometry for the next test
  result=reload_mesh( filename.c_str(), input_set);
  if (gen::error(moab::MB_SUCCESS!=result, "could not reload the mesh" )) return result; 

  // retrieve the verticies again so the model can be broken
  result = MBI()->get_entities_by_dimension(input_set, dim, verts, false);
  if(gen::error(moab::MB_SUCCESS!=result, " could not get vertices from the mesh")) return result;

  
  // "Break" the geometry
  result=rand_vert_bump(verts, facet_tolerance, root_name); 
  if(gen::error(moab::MB_SUCCESS!=result, "could not create single vertex bump test")) return result;
  
  // Seal the mesh
  result=mw_func::make_mesh_watertight (input_set, facet_tolerance, verbose);
  if(gen::error(moab::MB_SUCCESS!=result, "could not make the mesh watertight")) return result;

// Lastly Check to see if make_watertight fixed the model
  result=cw_func::check_mesh_for_watertightness( input_set, facet_tolerance, sealed, test);
  if(gen::error(moab::MB_SUCCESS!=result, "could not check model for watertightness")) return result;
  

  if(sealed)
  {
   std::cout << "PASS" << std::endl << std::endl ;
  }
  else
  {
   std::cout << "FAIL" << std::endl << std::endl ;
   test_set_result=false;
  }

//////////END SINGLE VERTEX MOVEMENT TESTS///////////

  if(test_set_result)
  {
   std::cout << "SINGLE VERTEX MOVEMENT TESTS PASSED" << std::endl << std::endl;
  }
  else
  {
   std::cout << "SINGLE VERTEX MOVEMENT TESTS FAILED" << std::endl << std::endl;
   exit(0);
  }

//////////LOCKED VERTEX PAIR MOVEMENT TESTS/////////////////

std::cout << "LOCKED VERTEX PAIR MOVEMENT TESTS" << std::endl;

///////////////BEGIN 1ST TEST////////////////////////

  std::cout << "Test 1: locked pair of verticies move in x:" ;

  // Clear the mesh and reload original geometry for the next test
  result=reload_mesh( filename.c_str(), input_set);
  if (gen::error(moab::MB_SUCCESS!=result, "could not reload the mesh" )) return result; 

  // retrieve the verticies again so the model can be broken
  result = MBI()->get_entities_by_dimension(input_set, dim, verts, false);
  if(gen::error(moab::MB_SUCCESS!=result, " could not get vertices from the mesh")) return result;

  
  // "Break" the geometry
  result=locked_pair_bump(verts, 0.9*facet_tolerance , 0 , 0 , root_name ); 
  if(gen::error(moab::MB_SUCCESS!=result, "could not create locked pair vertex bump test")) return result;


   // Seal the mesh
  result=mw_func::make_mesh_watertight (input_set, facet_tolerance, verbose);
  if(gen::error(moab::MB_SUCCESS!=result, "could not make the mesh watertight")) return result;

  // Lastly Check to see if make_watertight fixed the model
  result=cw_func::check_mesh_for_watertightness( input_set, facet_tolerance, sealed, test);
  if(gen::error(moab::MB_SUCCESS!=result, "could not check model for watertightness")) return result;
  
    if(sealed)
  {
   std::cout << "PASS" << std::endl;
  }
  else
  {
   std::cout << "FAIL" << std::endl;
   test_set_result=false;
  }

///////////////BEGIN 2ND TEST////////////////////////

  std::cout << "Test 2: locked pair of verticies move in y:" ;

  // Clear the mesh and reload original geometry for the next test
  result=reload_mesh( filename.c_str(), input_set);
  if (gen::error(moab::MB_SUCCESS!=result, "could not reload the mesh" )) return result; 

  // retrieve the verticies again so the model can be broken
  result = MBI()->get_entities_by_dimension(input_set, dim, verts, false);
  if(gen::error(moab::MB_SUCCESS!=result, " could not get vertices from the mesh")) return result;

  
  // "Break" the geometry
  result=locked_pair_bump(verts, 0 , 0.9*facet_tolerance , 0 , root_name ); 
  if(gen::error(moab::MB_SUCCESS!=result, "could not create locked pair vertex bump test")) return result;


   // Seal the mesh
  result=mw_func::make_mesh_watertight (input_set, facet_tolerance, verbose);
  if(gen::error(moab::MB_SUCCESS!=result, "could not make the mesh watertight")) return result;

  // Lastly Check to see if make_watertight fixed the model
  result=cw_func::check_mesh_for_watertightness( input_set, facet_tolerance, sealed, test);
  if(gen::error(moab::MB_SUCCESS!=result, "could not check model for watertightness")) return result;
  
    if(sealed)
  {
   std::cout << "PASS" << std::endl;
  }
  else
  {
   std::cout << "FAIL" << std::endl;
   test_set_result=false;
  }


///////////////BEGIN 3RD TEST////////////////////////

  std::cout << "Test 3: locked pair of verticies move in z:" ;

  // Clear the mesh and reload original geometry for the next test
  result=reload_mesh( filename.c_str(), input_set);
  if (gen::error(moab::MB_SUCCESS!=result, "could not reload the mesh" )) return result; 

  // retrieve the verticies again so the model can be broken
  result = MBI()->get_entities_by_dimension(input_set, dim, verts, false);
  if(gen::error(moab::MB_SUCCESS!=result, " could not get vertices from the mesh")) return result;

  
  // "Break" the geometry
  result=locked_pair_bump(verts, 0 , 0 , 0.9*facet_tolerance , root_name ); 
  if(gen::error(moab::MB_SUCCESS!=result, "could not create locked pair vertex bump test")) return result;


   // Seal the mesh
  result=mw_func::make_mesh_watertight (input_set, facet_tolerance, verbose);
  if(gen::error(moab::MB_SUCCESS!=result, "could not make the mesh watertight")) return result;

  // Lastly Check to see if make_watertight fixed the model
  result=cw_func::check_mesh_for_watertightness( input_set, facet_tolerance, sealed, test);
  if(gen::error(moab::MB_SUCCESS!=result, "could not check model for watertightness")) return result;
  
    if(sealed)
  {
   std::cout << "PASS" << std::endl;
  }
  else
  {
   std::cout << "FAIL" << std::endl;
   test_set_result=false;
  }

///////////////BEGIN 4TH TEST////////////////////////

  std::cout << "Test 4: locked pair of verticies move in R:" ;

  // Clear the mesh and reload original geometry for the next test
  result=reload_mesh( filename.c_str(), input_set);
  if (gen::error(moab::MB_SUCCESS!=result, "could not reload the mesh" )) return result; 

  // retrieve the verticies again so the model can be broken
  result = MBI()->get_entities_by_dimension(input_set, dim, verts, false);
  if(gen::error(moab::MB_SUCCESS!=result, " could not get vertices from the mesh")) return result;

  
  // "Break" the geometry
  result=locked_pair_bump_R(verts, facet_tolerance , root_name); 
  if(gen::error(moab::MB_SUCCESS!=result, "could not create locked pair vertex bump test")) return result;

  
  // Seal the mesh
  result=mw_func::make_mesh_watertight (input_set, facet_tolerance, verbose);
  if(gen::error(moab::MB_SUCCESS!=result, "could not make the mesh watertight")) return result;

  // Lastly Check to see if make_watertight fixed the model
  result=cw_func::check_mesh_for_watertightness( input_set, facet_tolerance, sealed, test);
  if(gen::error(moab::MB_SUCCESS!=result, "could not check model for watertightness")) return result;
  
    if(sealed)
  {
   std::cout << "PASS" << std::endl;
  }
  else
  {
   std::cout << "FAIL" << std::endl;
   test_set_result=false;
  }

///////////////BEGIN 5TH TEST////////////////////////

  std::cout << "Test 5: locked pair of verticies move in theta:" ;

  // Clear the mesh and reload original geometry for the next test
  result=reload_mesh( filename.c_str(), input_set);
  if (gen::error(moab::MB_SUCCESS!=result, "could not reload the mesh" )) return result; 

  // retrieve the verticies again so the model can be broken
  result = MBI()->get_entities_by_dimension(input_set, dim, verts, false);
  if(gen::error(moab::MB_SUCCESS!=result, " could not get vertices from the mesh")) return result;

  
  // "Break" the geometry
  result=locked_pair_move_theta(verts, facet_tolerance , root_name); 
  if(gen::error(moab::MB_SUCCESS!=result, "could not create locked pair vertex bump test")) return result;


   // Seal the mesh
  result=mw_func::make_mesh_watertight (input_set, facet_tolerance, verbose);
  if(gen::error(moab::MB_SUCCESS!=result, "could not make the mesh watertight")) return result;

  // Lastly Check to see if make_watertight fixed the model
  result=cw_func::check_mesh_for_watertightness( input_set, facet_tolerance, sealed, test);
  if(gen::error(moab::MB_SUCCESS!=result, "could not check model for watertightness")) return result;
  
    if(sealed)
  {
   std::cout << "PASS" << std::endl;
  }
  else
  {
   std::cout << "FAIL" << std::endl;
   test_set_result=false;
  }

///////////////BEGIN 6TH TEST////////////////////////

  std::cout << "Test 6: adj. pair of verticies random move:" ;

  // Clear the mesh and reload original geometry for the next test
  result=reload_mesh( filename.c_str(), input_set);
  if (gen::error(moab::MB_SUCCESS!=result, "could not reload the mesh" )) return result; 

  // retrieve the verticies again so the model can be broken
  result = MBI()->get_entities_by_dimension(input_set, dim, verts, false);
  if(gen::error(moab::MB_SUCCESS!=result, " could not get vertices from the mesh")) return result;

  
  // "Break" the geometry
  result=locked_pair_bump_rand(verts, facet_tolerance , root_name); 
  if(gen::error(moab::MB_SUCCESS!=result, "could not create locked pair vertex bump test")) return result;

  

   // Seal the mesh
  result=mw_func::make_mesh_watertight (input_set, facet_tolerance, verbose);
  if(gen::error(moab::MB_SUCCESS!=result, "could not make the mesh watertight")) return result;


  // Lastly Check to see if make_watertight fixed the model
  result=cw_func::check_mesh_for_watertightness( input_set, facet_tolerance, sealed, test);
  if(gen::error(moab::MB_SUCCESS!=result, "could not check model for watertightness")) return result;
  

  if(sealed)
  {
   std::cout << "PASS" << std::endl << std::endl ;
  }
  else
  {
   std::cout << "FAIL" << std::endl << std::endl ;
   test_set_result=false;
  }


//////////END LOCKED VERTEX PAIR MOVEMENT TESTS///////////

  if(test_set_result)
  {
   std::cout << "LOCKED VERTEX PAIR MOVEMENT TESTS PASSED" << std::endl << std::endl;
  }
  else
  {
   std::cout << "LOCKED VERTEX PAIR MOVEMENT TESTS FAILED" << std::endl << std::endl;
   exit(0);
  }


////////////RAND PAIR MOVEMENT TESTS//////////////////

std::cout << "RAND PAIR MOVEMENT TESTS" << std::endl;

///////////////BEGIN 1ST TEST////////////////////////

  std::cout << "Test 1: random locked pair of verticies move in x:" ;

  // Clear the mesh and reload original geometry for the next test
  result=reload_mesh( filename.c_str(), input_set);
  if (gen::error(moab::MB_SUCCESS!=result, "could not reload the mesh" )) return result; 

  // retrieve the verticies again so the model can be broken
  result = MBI()->get_entities_by_dimension(input_set, dim, verts, false);
  if(gen::error(moab::MB_SUCCESS!=result, " could not get vertices from the mesh")) return result;
  
  // "Break" the geometry
  result=rand_locked_pair_bump(verts, 0.9*facet_tolerance , 0 , 0 , root_name ); 
  if(gen::error(moab::MB_SUCCESS!=result, "could not create locked pair vertex bump test")) return result;

  // Seal the mesh
  result=mw_func::make_mesh_watertight (input_set, facet_tolerance, verbose);
  if(gen::error(moab::MB_SUCCESS!=result, "could not make the mesh watertight")) return result;

  // Lastly Check to see if make_watertight fixed the model
  result=cw_func::check_mesh_for_watertightness( input_set, facet_tolerance, sealed, test);
  if(gen::error(moab::MB_SUCCESS!=result, "could not check model for watertightness")) return result;
  
    if(sealed)
  {
   std::cout << "PASS" << std::endl;
  }
  else
  {
   std::cout << "FAIL" << std::endl;
   test_set_result=false;
  }

///////////////BEGIN 2ND TEST////////////////////////

  std::cout << "Test 2: random locked pair of verticies move in y:" ;

  // Clear the mesh and reload original geometry for the next test
  result=reload_mesh( filename.c_str(), input_set);
  if (gen::error(moab::MB_SUCCESS!=result, "could not reload the mesh" )) return result; 

  // retrieve the verticies again so the model can be broken
  result = MBI()->get_entities_by_dimension(input_set, dim, verts, false);
  if(gen::error(moab::MB_SUCCESS!=result, " could not get vertices from the mesh")) return result;

  
  // "Break" the geometry
  result=rand_locked_pair_bump(verts, 0 , 0.9*facet_tolerance , 0 , root_name ); 
  if(gen::error(moab::MB_SUCCESS!=result, "could not create locked pair vertex bump test")) return result;

  // Seal the mesh
  result=mw_func::make_mesh_watertight (input_set, facet_tolerance, verbose);
  if(gen::error(moab::MB_SUCCESS!=result, "could not make the mesh watertight")) return result;

  // Lastly Check to see if make_watertight fixed the model
  result=cw_func::check_mesh_for_watertightness( input_set, facet_tolerance, sealed, test);
  if(gen::error(moab::MB_SUCCESS!=result, "could not check model for watertightness")) return result;
  
    if(sealed)
  {
   std::cout << "PASS" << std::endl;
  }
  else
  {
   std::cout << "FAIL" << std::endl;
   test_set_result=false;
  }

///////////////BEGIN 3RD TEST////////////////////////

  std::cout << "Test 3: random locked pair of verticies move in z:" ;

  // Clear the mesh and reload original geometry for the next test
  result=reload_mesh( filename.c_str(), input_set);
  if (gen::error(moab::MB_SUCCESS!=result, "could not reload the mesh" )) return result; 

  // retrieve the verticies again so the model can be broken
  result = MBI()->get_entities_by_dimension(input_set, dim, verts, false);
  if(gen::error(moab::MB_SUCCESS!=result, " could not get vertices from the mesh")) return result;
  
  // "Break" the geometry
  result=rand_locked_pair_bump(verts, 0 , 0 , 0.9*facet_tolerance , root_name ); 
  if(gen::error(moab::MB_SUCCESS!=result, "could not create locked pair vertex bump test")) return result;

  // Seal the mesh
  result=mw_func::make_mesh_watertight (input_set, facet_tolerance, verbose);
  if(gen::error(moab::MB_SUCCESS!=result, "could not make the mesh watertight")) return result;

  // Lastly Check to see if make_watertight fixed the model
  result=cw_func::check_mesh_for_watertightness( input_set, facet_tolerance, sealed, test);
  if(gen::error(moab::MB_SUCCESS!=result, "could not check model for watertightness")) return result;
  
    if(sealed)
  {
   std::cout << "PASS" << std::endl;
  }
  else
  {
   std::cout << "FAIL" << std::endl;
   test_set_result=false;
  }

///////////////BEGIN 4TH TEST////////////////////////

  std::cout << "Test 4: random locked pair of verticies move in R:" ;

  // Clear the mesh and reload original geometry for the next test
  result=reload_mesh( filename.c_str(), input_set);
  if (gen::error(moab::MB_SUCCESS!=result, "could not reload the mesh" )) return result; 

  // retrieve the verticies again so the model can be broken
  result = MBI()->get_entities_by_dimension(input_set, dim, verts, false);
  if(gen::error(moab::MB_SUCCESS!=result, " could not get vertices from the mesh")) return result;

  // "Break" the geometry
  result=rand_locked_pair_bump_R(verts, facet_tolerance , root_name); 
  if(gen::error(moab::MB_SUCCESS!=result, "could not create locked pair vertex bump test")) return result;

  // Seal the mesh
  result=mw_func::make_mesh_watertight (input_set, facet_tolerance, verbose);
  if(gen::error(moab::MB_SUCCESS!=result, "could not make the mesh watertight")) return result;

  // Lastly Check to see if make_watertight fixed the model
  result=cw_func::check_mesh_for_watertightness( input_set, facet_tolerance, sealed, test);
  if(gen::error(moab::MB_SUCCESS!=result, "could not check model for watertightness")) return result;
  
    if(sealed)
  {
   std::cout << "PASS" << std::endl;
  }
  else
  {
   std::cout << "FAIL" << std::endl;
   test_set_result=false;
  }

///////////////BEGIN 5TH TEST////////////////////////

  std::cout << "Test 5: random locked pair of verticies move in theta:" ;

  // Clear the mesh and reload original geometry for the next test
  result=reload_mesh( filename.c_str(), input_set);
  if (gen::error(moab::MB_SUCCESS!=result, "could not reload the mesh" )) return result; 

  // retrieve the verticies again so the model can be broken
  result = MBI()->get_entities_by_dimension(input_set, dim, verts, false);
  if(gen::error(moab::MB_SUCCESS!=result, " could not get vertices from the mesh")) return result;

  // "Break" the geometry
  result=rand_locked_pair_bump_theta(verts, facet_tolerance , root_name); 
  if(gen::error(moab::MB_SUCCESS!=result, "could not create locked pair vertex bump test")) return result;

  // Seal the mesh
  result=mw_func::make_mesh_watertight (input_set, facet_tolerance, verbose);
  if(gen::error(moab::MB_SUCCESS!=result, "could not make the mesh watertight")) return result;

  // Lastly Check to see if make_watertight fixed the model
  result=cw_func::check_mesh_for_watertightness( input_set, facet_tolerance, sealed, test);
  if(gen::error(moab::MB_SUCCESS!=result, "could not check model for watertightness")) return result;
  
    if(sealed)
  {
   std::cout << "PASS" << std::endl;
  }
  else
  {
   std::cout << "FAIL" << std::endl;
   test_set_result=false;
  }

///////////////BEGIN 6TH TEST////////////////////////

  std::cout << "Test 6: random pair of verticies move in random dir:" ;

  // Clear the mesh and reload original geometry for the next test
  result=reload_mesh( filename.c_str(), input_set);
  if (gen::error(moab::MB_SUCCESS!=result, "could not reload the mesh" )) return result; 

  // retrieve the verticies again so the model can be broken
  result = MBI()->get_entities_by_dimension(input_set, dim, verts, false);
  if(gen::error(moab::MB_SUCCESS!=result, " could not get vertices from the mesh")) return result;

  // "Break" the geometry
  result=rand_locked_pair_bump_rand(verts, facet_tolerance , root_name); 
  if(gen::error(moab::MB_SUCCESS!=result, "could not create locked pair vertex bump test")) return result;

  // Seal the mesh
  result=mw_func::make_mesh_watertight (input_set, facet_tolerance, verbose);
  if(gen::error(moab::MB_SUCCESS!=result, "could not make the mesh watertight")) return result;

  // Lastly Check to see if make_watertight fixed the model
  result=cw_func::check_mesh_for_watertightness( input_set, facet_tolerance, sealed, test);
  if(gen::error(moab::MB_SUCCESS!=result, "could not check model for watertightness")) return result;
  
    if(sealed)
  {
   std::cout << "PASS" << std::endl << std::endl ;
  }
  else
  {
   std::cout << "FAIL" << std::endl << std::endl ;
   test_set_result=false;
  }


//////////END RAND LOCKED VERTEX PAIR MOVEMENT TESTS///////////

  if(test_set_result)
  {
   std::cout << "RANDOM LOCKED VERTEX PAIR MOVEMENT TESTS PASSED" << std::endl << std::endl;
  }
  else
  {
   std::cout << "RANDOM LOCKED VERTEX PAIR MOVEMENT TESTS FAILED" << std::endl << std::endl;
   exit(0);
  }


/////////////ADJACENT PLUS ONE TESTS//////////////////

std::cout << "ADJACENT PLUS ONE TESTS" << std::endl;

///////////////BEGIN 1ST TEST////////////////////////

  std::cout << "Test 1: locked pair of verticies (adj+1) move in x:" ;

  // Clear the mesh and reload original geometry for the next test
  result=reload_mesh( filename.c_str(), input_set);
  if (gen::error(moab::MB_SUCCESS!=result, "could not reload the mesh" )) return result; 

  // retrieve the verticies again so the model can be broken
  result = MBI()->get_entities_by_dimension(input_set, dim, verts, false);
  if(gen::error(moab::MB_SUCCESS!=result, " could not get vertices from the mesh")) return result;

  //"Break" the geometry
  result=adjplone_locked_pair_bump(verts, 0.9*facet_tolerance , 0 , 0 , root_name ); 
  if(gen::error(moab::MB_SUCCESS!=result, "could not create locked pair vertex bump test")) return result;

  // Seal the mesh
  result=mw_func::make_mesh_watertight (input_set, facet_tolerance, verbose);
  if(gen::error(moab::MB_SUCCESS!=result, "could not make the mesh watertight")) return result;

  // Lastly Check to see if make_watertight fixed the model
  result=cw_func::check_mesh_for_watertightness( input_set, facet_tolerance, sealed, test);
  if(gen::error(moab::MB_SUCCESS!=result, "could not check model for watertightness")) return result;
  
    if(sealed)
  {
   std::cout << "PASS" << std::endl;
  }
  else
  {
   std::cout << "FAIL" << std::endl;
   test_set_result=false;
  }

///////////////BEGIN 2ND TEST////////////////////////

  std::cout << "Test 2: locked pair of verticies (adj+1) move in y:" ;

  // Clear the mesh and reload original geometry for the next test
  result=reload_mesh( filename.c_str(), input_set);
  if (gen::error(moab::MB_SUCCESS!=result, "could not reload the mesh" )) return result; 

  // retrieve the verticies again so the model can be broken
  result = MBI()->get_entities_by_dimension(input_set, dim, verts, false);
  if(gen::error(moab::MB_SUCCESS!=result, " could not get vertices from the mesh")) return result;

  //"Break" the geometry
  result=adjplone_locked_pair_bump(verts, 0 , 0.9*facet_tolerance , 0 , root_name ); 
  if(gen::error(moab::MB_SUCCESS!=result, "could not create locked pair vertex bump test")) return result;

  // Seal the mesh
  result=mw_func::make_mesh_watertight (input_set, facet_tolerance, verbose);
  if(gen::error(moab::MB_SUCCESS!=result, "could not make the mesh watertight")) return result;

  // Lastly Check to see if make_watertight fixed the model
  result=cw_func::check_mesh_for_watertightness( input_set, facet_tolerance, sealed, test);
  if(gen::error(moab::MB_SUCCESS!=result, "could not check model for watertightness")) return result;
  
    if(sealed)
  {
   std::cout << "PASS" << std::endl;
  }
  else
  {
   std::cout << "FAIL" << std::endl;
   test_set_result=false;
  }

///////////////BEGIN 3RD TEST////////////////////////

  std::cout << "Test 3: locked pair of verticies (adj+1) move in z:" ;

  // Clear the mesh and reload original geometry for the next test
  result=reload_mesh( filename.c_str(), input_set);
  if (gen::error(moab::MB_SUCCESS!=result, "could not reload the mesh" )) return result; 

  // retrieve the verticies again so the model can be broken
  result = MBI()->get_entities_by_dimension(input_set, dim, verts, false);
  if(gen::error(moab::MB_SUCCESS!=result, " could not get vertices from the mesh")) return result;

  //"Break" the geometry
  result=adjplone_locked_pair_bump(verts, 0 , 0 , 0.9*facet_tolerance , root_name ); 
  if(gen::error(moab::MB_SUCCESS!=result, "could not create locked pair vertex bump test")) return result;

  // Seal the mesh
  result=mw_func::make_mesh_watertight (input_set, facet_tolerance, verbose);
  if(gen::error(moab::MB_SUCCESS!=result, "could not make the mesh watertight")) return result;

  // Lastly Check to see if make_watertight fixed the model
  result=cw_func::check_mesh_for_watertightness( input_set, facet_tolerance, sealed, test);
  if(gen::error(moab::MB_SUCCESS!=result, "could not check model for watertightness")) return result;
  
    if(sealed)
  {
   std::cout << "PASS" << std::endl;
  }
  else
  {
   std::cout << "FAIL" << std::endl;
   test_set_result=false;
  }

///////////////BEGIN 4TH TEST////////////////////////

  std::cout << "Test 4: locked pair of verticies (adj+1) move in R:" ;

  // Clear the mesh and reload original geometry for the next test
  result=reload_mesh( filename.c_str(), input_set);
  if (gen::error(moab::MB_SUCCESS!=result, "could not reload the mesh" )) return result; 

  // retrieve the verticies again so the model can be broken
  result = MBI()->get_entities_by_dimension(input_set, dim, verts, false);
  if(gen::error(moab::MB_SUCCESS!=result, " could not get vertices from the mesh")) return result;

  //"Break" the geometry
  result=adjplone_locked_pair_bump_R(verts, facet_tolerance , root_name ); 
  if(gen::error(moab::MB_SUCCESS!=result, "could not create locked pair vertex bump test")) return result;
 
  // Seal the mesh
  result=mw_func::make_mesh_watertight (input_set, facet_tolerance, verbose);
  if(gen::error(moab::MB_SUCCESS!=result, "could not make the mesh watertight")) return result;

  // Lastly Check to see if make_watertight fixed the model
  result=cw_func::check_mesh_for_watertightness( input_set, facet_tolerance, sealed, test);
  if(gen::error(moab::MB_SUCCESS!=result, "could not check model for watertightness")) return result;
  
    if(sealed)
  {
   std::cout << "PASS" << std::endl;
  }
  else
  {
   std::cout << "FAIL" << std::endl;
   test_set_result=false;
  }
///////////////BEGIN 5TH TEST////////////////////////

  std::cout << "Test 5: locked pair of verticies (adj+1) move in theta:" ;

  // Clear the mesh and reload original geometry for the next test
  result=reload_mesh( filename.c_str(), input_set);
  if (gen::error(moab::MB_SUCCESS!=result, "could not reload the mesh" )) return result; 

  // retrieve the verticies again so the model can be broken
  result = MBI()->get_entities_by_dimension(input_set, dim, verts, false);
  if(gen::error(moab::MB_SUCCESS!=result, " could not get vertices from the mesh")) return result;

  //"Break" the geometry
  result=adjplone_locked_pair_bump_theta(verts, facet_tolerance , root_name ); 
  if(gen::error(moab::MB_SUCCESS!=result, "could not create locked pair vertex bump test")) return result;

  // Seal the mesh
  result=mw_func::make_mesh_watertight (input_set, facet_tolerance, verbose);
  if(gen::error(moab::MB_SUCCESS!=result, "could not make the mesh watertight")) return result;

  // Lastly Check to see if make_watertight fixed the model
  result=cw_func::check_mesh_for_watertightness( input_set, facet_tolerance, sealed, test);
  if(gen::error(moab::MB_SUCCESS!=result, "could not check model for watertightness")) return result;
  
    if(sealed)
  {
   std::cout << "PASS" << std::endl;
  }
  else
  {
   std::cout << "FAIL" << std::endl;
   test_set_result=false;
  }

///////////////BEGIN 6TH TEST////////////////////////

  std::cout << "Test 6: pair of verticies (adj+1) move in random dir:" ;

  // Clear the mesh and reload original geometry for the next test
  result=reload_mesh( filename.c_str(), input_set);
  if (gen::error(moab::MB_SUCCESS!=result, "could not reload the mesh" )) return result; 

  // retrieve the verticies again so the model can be broken
  result = MBI()->get_entities_by_dimension(input_set, dim, verts, false);
  if(gen::error(moab::MB_SUCCESS!=result, " could not get vertices from the mesh")) return result;

  //"Break" the geometry
  result=adjplone_locked_pair_bump_rand(verts, facet_tolerance , root_name ); 
  if(gen::error(moab::MB_SUCCESS!=result, "could not create locked pair vertex bump test")) return result;

  // Seal the mesh
  result=mw_func::make_mesh_watertight (input_set, facet_tolerance, verbose);
  if(gen::error(moab::MB_SUCCESS!=result, "could not make the mesh watertight")) return result;

  // Lastly Check to see if make_watertight fixed the model
  result=cw_func::check_mesh_for_watertightness( input_set, facet_tolerance, sealed, test);
  if(gen::error(moab::MB_SUCCESS!=result, "could not check model for watertightness")) return result;
  
    if(sealed)
  {
   std::cout << "PASS" << std::endl << std::endl ;
  }
  else
   
  {
   std::cout << "FAIL" << std::endl << std::endl ;
   test_set_result=false;
  }

//////////END ADJACENT + 1 VERTEX MOVEMENT TESTS///////////

  if(test_set_result)
  {
   std::cout << "ADJACENT PLUS ONE VERTEX PAIR MOVEMENT TESTS PASSED" << std::endl << std::endl;
  }
  else
  {
   std::cout << "ADJACENT PLUS ONE VERTEX PAIR MOVEMENT TESTS FAILED" << std::endl << std::endl;
   exit(0);
  }


//////////NON-ADJACENT LOCKED PAIR TESTS//////////////

std::cout << "NON-ADJACENT LOCKED PAIR TESTS" << std::endl;

///////////////BEGIN 1ST TEST////////////////////////

  std::cout << "Test 1: non-adjacent locked pair of verticies move in x:" ;

  // Clear the mesh and reload original geometry for the next test
  result=reload_mesh( filename.c_str(), input_set);
  if (gen::error(moab::MB_SUCCESS!=result, "could not reload the mesh" )) return result; 

  // retrieve the verticies again so the model can be broken
  result = MBI()->get_entities_by_dimension(input_set, dim, verts, false);
  if(gen::error(moab::MB_SUCCESS!=result, " could not get vertices from the mesh")) return result;

  //"Break" the geometry
  result= nonadj_locked_pair_bump(verts,  0.9*facet_tolerance , 0 , 0 , root_name ); 
  if(gen::error(moab::MB_SUCCESS!=result, "could not create locked pair vertex bump test")) return result;

  // Seal the mesh
  result=mw_func::make_mesh_watertight (input_set, facet_tolerance, verbose);
  if(gen::error(moab::MB_SUCCESS!=result, "could not make the mesh watertight")) return result;

  // Lastly Check to see if make_watertight fixed the model
  result=cw_func::check_mesh_for_watertightness( input_set, facet_tolerance, sealed, test);
  if(gen::error(moab::MB_SUCCESS!=result, "could not check model for watertightness")) return result;
  
    if(sealed)
  {
   std::cout << "PASS" << std::endl;
  }
  else
  {
   std::cout << "FAIL" << std::endl;
   test_set_result=false;
  }

///////////////BEGIN 2ND TEST////////////////////////

  std::cout << "Test 2: non-adjacent locked pair of verticies move in y:" ;

  // Clear the mesh and reload original geometry for the next test
  result=reload_mesh( filename.c_str(), input_set);
  if (gen::error(moab::MB_SUCCESS!=result, "could not reload the mesh" )) return result; 

  // retrieve the verticies again so the model can be broken
  result = MBI()->get_entities_by_dimension(input_set, dim, verts, false);
  if(gen::error(moab::MB_SUCCESS!=result, " could not get vertices from the mesh")) return result;

  //"Break" the geometry
  result= nonadj_locked_pair_bump(verts, 0 ,  0.9*facet_tolerance , 0 , root_name ); 
  if(gen::error(moab::MB_SUCCESS!=result, "could not create locked pair vertex bump test")) return result;

  // Seal the mesh
  result=mw_func::make_mesh_watertight (input_set, facet_tolerance, verbose);
  if(gen::error(moab::MB_SUCCESS!=result, "could not make the mesh watertight")) return result;

  // Lastly Check to see if make_watertight fixed the model
  result=cw_func::check_mesh_for_watertightness( input_set, facet_tolerance, sealed, test);
  if(gen::error(moab::MB_SUCCESS!=result, "could not check model for watertightness")) return result;
  
    if(sealed)
  {
   std::cout << "PASS" << std::endl;
  }
  else
  {
   std::cout << "FAIL" << std::endl;
   test_set_result=false;
  }


///////////////BEGIN 3RD TEST////////////////////////

  std::cout << "Test 3: non-adjacent locked pair of verticies move in z:" ;

  // Clear the mesh and reload original geometry for the next test
  result=reload_mesh( filename.c_str(), input_set);
  if (gen::error(moab::MB_SUCCESS!=result, "could not reload the mesh" )) return result; 

  // retrieve the verticies again so the model can be broken
  result = MBI()->get_entities_by_dimension(input_set, dim, verts, false);
  if(gen::error(moab::MB_SUCCESS!=result, " could not get vertices from the mesh")) return result;

  //"Break" the geometry
  result= nonadj_locked_pair_bump(verts, 0 , 0 ,  0.9*facet_tolerance  , root_name ); 
  if(gen::error(moab::MB_SUCCESS!=result, "could not create locked pair vertex bump test")) return result;

  // Seal the mesh
  result=mw_func::make_mesh_watertight (input_set, facet_tolerance, verbose);
  if(gen::error(moab::MB_SUCCESS!=result, "could not make the mesh watertight")) return result;

  // Lastly Check to see if make_watertight fixed the model
  result=cw_func::check_mesh_for_watertightness( input_set, facet_tolerance, sealed, test);
  if(gen::error(moab::MB_SUCCESS!=result, "could not check model for watertightness")) return result;
  
    if(sealed)
  {
   std::cout << "PASS" << std::endl;
  }
  else
  {
   std::cout << "FAIL" << std::endl;
   test_set_result=false;
  }

///////////////BEGIN 4TH TEST////////////////////////

  std::cout << "Test 4: non-adjacent locked pair of verticies move in R:" ;

  // Clear the mesh and reload original geometry for the next test
  result=reload_mesh( filename.c_str(), input_set);
  if (gen::error(moab::MB_SUCCESS!=result, "could not reload the mesh" )) return result; 

  // retrieve the verticies again so the model can be broken
  result = MBI()->get_entities_by_dimension(input_set, dim, verts, false);
  if(gen::error(moab::MB_SUCCESS!=result, " could not get vertices from the mesh")) return result;

  //"Break" the geometry
  result= nonadj_locked_pair_bump_R(verts, facet_tolerance , root_name ); 
  if(gen::error(moab::MB_SUCCESS!=result, "could not create locked pair vertex bump test")) return result;

  // write a new file for checking broken geometry
  result = write_mod_file( root_name);
  if(gen::error(moab::MB_SUCCESS!=result, "could not write new file")) return result;

  // Seal the mesh
  result=mw_func::make_mesh_watertight (input_set, facet_tolerance, verbose);
  if(gen::error(moab::MB_SUCCESS!=result, "could not make the mesh watertight")) return result;

  // Lastly Check to see if make_watertight fixed the model
  result=cw_func::check_mesh_for_watertightness( input_set, facet_tolerance, sealed, test);
  if(gen::error(moab::MB_SUCCESS!=result, "could not check model for watertightness")) return result;
  
    if(sealed)
  {
   std::cout << "PASS" << std::endl;
  }
  else
  {
   std::cout << "FAIL" << std::endl;
   test_set_result=false;
  }

///////////////BEGIN 5TH TEST////////////////////////

  std::cout << "Test 5: non-adjacent locked pair of verticies move in theta:" ;

  // Clear the mesh and reload original geometry for the next test
  result=reload_mesh( filename.c_str(), input_set);
  if (gen::error(moab::MB_SUCCESS!=result, "could not reload the mesh" )) return result; 

  // retrieve the verticies again so the model can be broken
  result = MBI()->get_entities_by_dimension(input_set, dim, verts, false);
  if(gen::error(moab::MB_SUCCESS!=result, " could not get vertices from the mesh")) return result;

  //"Break" the geometry
  result= nonadj_locked_pair_bump_theta(verts, facet_tolerance , root_name ); 
  if(gen::error(moab::MB_SUCCESS!=result, "could not create locked pair vertex bump test")) return result;

  // Seal the mesh
  result=mw_func::make_mesh_watertight (input_set, facet_tolerance, verbose);
  if(gen::error(moab::MB_SUCCESS!=result, "could not make the mesh watertight")) return result;

  // Lastly Check to see if make_watertight fixed the model
  result=cw_func::check_mesh_for_watertightness( input_set, facet_tolerance, sealed, test);
  if(gen::error(moab::MB_SUCCESS!=result, "could not check model for watertightness")) return result;
  
    if(sealed)
  {
   std::cout << "PASS" << std::endl;
  }
  else
  {
   std::cout << "FAIL" << std::endl;
   test_set_result=false;
  }


///////////////BEGIN 6TH TEST////////////////////////

  std::cout << "Test 6: non-adjacent pair of verticies move in random dir:" ;

  // Clear the mesh and reload original geometry for the next test
  result=reload_mesh( filename.c_str(), input_set);
  if (gen::error(moab::MB_SUCCESS!=result, "could not reload the mesh" )) return result; 

  // retrieve the verticies again so the model can be broken
  result = MBI()->get_entities_by_dimension(input_set, dim, verts, false);
  if(gen::error(moab::MB_SUCCESS!=result, " could not get vertices from the mesh")) return result;

  //"Break" the geometry
  result= nonadj_locked_pair_bump_rand(verts, facet_tolerance , root_name ); 
  if(gen::error(moab::MB_SUCCESS!=result, "could not create locked pair vertex bump test")) return result;

  // Seal the mesh
  result=mw_func::make_mesh_watertight (input_set, facet_tolerance, verbose);
  if(gen::error(moab::MB_SUCCESS!=result, "could not make the mesh watertight")) return result;

  // Lastly Check to see if make_watertight fixed the model
  result=cw_func::check_mesh_for_watertightness( input_set, facet_tolerance, sealed, test);
  if(gen::error(moab::MB_SUCCESS!=result, "could not check model for watertightness")) return result;
  
    if(sealed)
  {
   std::cout << "PASS" << std::endl << std::endl ;
  }
  else
  {
   std::cout << "FAIL" << std::endl << std::endl ;
   test_set_result=false;
  }


//////////END NON-ADJACNET LOCKED VERTEX PAIR MOVEMENT TESTS///////////

  if(test_set_result)
  {
   std::cout << "NON-ADJACENT LOCKED VERTEX PAIR MOVEMENT TESTS PASSED" << std::endl << std::endl;
  }
  else
  {
   std::cout << "NON-ADJACENT LOCKED VERTEX PAIR MOVEMENT TESTS FAILED" << std::endl << std::endl;
   exit(0);
  }


}


