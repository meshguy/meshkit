/*
 * coregen program
 *
 * Program for generating sixth-core hexagonal lattice of hexagonal assemblies
 * in a fast reactor core.  The core is arranged around a central hexagonal
 * assembly (oriented with points of the hexagon pointing up and down, and flat
 * sides on the left and right), with subsequent rings of assemblies around 
 * that.  The inner-most assembly is ring #1.  The program allows an arbitrary
 * number of assembly types, each read from a separate mesh file.
 *
 * This program takes as input:
 * NRINGS: number of rings in the core
 * PACK_TYPE: 0: unit cell packing in hexagonal NRINGS-ring lattice, with
 *               an outer duct wall (specified in second mesh file)
 *            1: assembly-type lattice
 * PITCH: horizontal pitch between unit cells or assemblies
 * SYMM: 1/symmetry, usually 6 here
 * assembly types in lattice: NRINGS*(NRINGS+1) assembly type integers,
 *   specifying the integer assembly types in the 1st row, 2nd row, etc.
 *   There are NRINGS assemblies in the 1st row, (NRINGS-1) in the 2nd row, 
 *   etc.  Each row above the 1st starts with x and y displacements of
 *   PITCH*cos(pi/3) and PITCH*sin(pi/3), resp, from the start of the last row.
 *   An assembly type of -1 indicates no assembly at this position.  Assembly
 *   types are indexed from zero.
 * assembly mesh files: the mesh file containing the mesh for each assembly
 *   type; should be n+1 file names, one on each line, where n is the maximum
 *   assembly index in the assembly types lattice above.
 * background assembly: y=yes, n=no
 * if background assembly = y, background mesh file
 * output mesh file: file name to which the final mesh is written
 *
 * Sample input file (after you remove the '*' on each line):
 * 
 * 9
 * 1
 * 14.685
 * 6
 * 6   4   4   4   4   2   2   3   -1
 *   5   1   0   0   2   2   3   -1
 *     5   0   1   0   2   2   3
 *       5   0   0   2   2   3
 *         5   2   2   2   3
 *           2   2   2   3
 *             2   3   3
 *               3   -1
 *                 -1
 * fuel_assy.h5m
 * ctrl_assy.h5m
 * fuel_assy.h5m
 * fuel_assy.h5m
 * fuel_assy_half.h5m
 * fuel_assy_half_angle.h5m
 * ctrl_assy_sixth.h5m
 * y
 * sodium_all.h5m
 * sixth.h5m
 * 
 * NOTE: Assembly types are arranged in a wedge, but the actual arrangement
 * of assemblies in the core will look more like an inverted wedge (a "delta"),
 * since assemblies are displaced in +x and +y direction in increasing rows.
 * Also, in the above mesh, a sixth-assembly is used at the center, half-
 * assemblies along the bottom row, and half_angle assemblies as the first
 * assembly in each row; this makes the mesh resolve the exact symmetry planes
 * of a 1/6 core.  In other words, you should be able to copy/rotate this 1/6
 * mesh 5 times, merge, and get a complete 360 degree core (though I haven't
 * verified this).
 * 
 */
#include <cfloat>
#include <math.h>
#include "MKUtils.hpp"
#include "MergeMesh.hpp"
#include "CopyMesh.hpp"
#include "parser.hpp"
#include "fileio.hpp"
#include "clock.hpp"
#include <sstream>
#include <string>
;
int get_copy_expand_sets(CopyMesh *cm,
                         iBase_EntitySetHandle orig_set, 
                         const char **tag_names, const char **tag_vals, 
                         int num_tags, int copy_or_expand);

int extend_expand_sets(CopyMesh *cm);

int copy_move_assys(CopyMesh *cm,
                    const int nrings, const int pack_type,
                    const double pitch,
                    const int symm,
                    std::vector<std::string> &core_alias,
                    std::vector<iBase_EntitySetHandle> &assys);

int copy_move_sq_assys(CopyMesh *cm,
		       const int nrings, const int pack_type,
		       const double pitch,
		       const int symm,
		        std::vector<std::string> &core_alias,
		       std::vector<iBase_EntitySetHandle> &assys);

int copy_move_hex_assys(CopyMesh *cm,
		       const int nrings, const int pack_type,
		       const double pitch,
		       const int symm,
		        std::vector<std::string> &core_alias,
		       std::vector<iBase_EntitySetHandle> &assys);

int read_inputs();

int banner();

int find_assm(const int i, int &assm_index);

int del_orig_mesh(std::vector<iBase_EntitySetHandle> &assys,
                  const bool back_mesh);

int prepareIO (int argc, char *argv[]);

int write_makefile();

iMesh_Instance impl;
iBase_EntitySetHandle root_set;
#define DEFAULT_TEST_FILE "twoassm"

const int UNITCELL_DUCT = 0, ASSY_TYPES = 1;

// declare variables read in the inputs
int nrings, nringsx, nringsy, pack_type =1, symm = 6;
double pitch, pitchx, pitchy;
bool global_ids = true, back_mesh;
std::vector<std::string> files;
std::string outfile;
int nassys; // the number of mesh files
int tot_assys; // total no. of assms forming core
int set_DIM = 3;

// file related
std::ifstream file_input;	// File Input
std::ofstream make_file;        // File Output
std::string iname, ifile, mfile, geometry, back_meshfile; 
int linenumber;

// variables in file
std::string card,geom_type, meshfile, mf_alias, temp_alias;
std::vector<std::string> assm_alias;
std::vector<std::string> core_alias;
  
  // parsing related
  std::string input_string;
  std::string comment = "!";
  int MAXCHARS = 300;




int main(int argc, char **argv) 
{
  int err;

  // print banner
  err = banner();

  // start the timer 
  CClock Timer;
  std::string szDateTime;
  Timer.GetDateTime (szDateTime);
  std::cout << "\nStarting out at : " << szDateTime << "\n";

  // read inputs and open makefile for writing
  err = prepareIO (argc, argv);
  err = read_inputs();
  err = write_makefile();
  // make a mesh instance
  iMesh_newMesh("MOAB", &impl, &err, 4);
  ERRORR("Failed to create instance.", 1);
  
  iMesh_getRootSet(impl, &root_set, &err);

  CopyMesh *cm = new CopyMesh(impl);
  MergeMesh *mm = new MergeMesh(impl);

  std::vector<iBase_EntitySetHandle> assys;
  iBase_EntitySetHandle orig_set;

  // loop reading all mesh files
  for (unsigned int i = 0; i < files.size(); i++) {
    iMesh_createEntSet(impl, 0, &orig_set, &err);
    ERRORR("Couldn't create file set.", err);
  
    iMesh_load(impl, orig_set, files[i].c_str(), NULL, &err, 
               strlen(files[i].c_str()), 0);
    ERRORR("Couldn't read mesh file.", err);

    assys.push_back(orig_set);
  }
  std::cout << "Loaded mesh files." << std::endl;

  // move the assys based on the geometry type
  if(!strcmp(geom_type.c_str(), "hexagonal")){
    err = copy_move_hex_assys(cm, nrings, pack_type, pitch, symm, 
			  core_alias, assys);
    ERRORR("Failed in copy/move step.", err);
  }
  else if(!strcmp(geom_type.c_str(),"cartesian")){
    err = copy_move_sq_assys(cm, nrings, pack_type, pitch, symm, 
			  core_alias, assys);
    ERRORR("Failed in copy/move step.", err);
  } 
  else if(!strcmp(geom_type.c_str(),"default")){
    err = copy_move_assys(cm, nrings, pack_type, pitch, symm, 
			  core_alias, assys);
    ERRORR("Failed in copy/move step.", err);
  }
  
  // delete all original meshes
  err = del_orig_mesh(assys, back_mesh);
  ERRORR("Failed in delete step.", err);

  //get entities for merging
  iBase_EntityHandle *ents = NULL; 
  int ents_alloc = 0, ents_size;

  if(set_DIM ==2){ // if surface geometry specified
    iMesh_getEntities(impl, root_set,
                    iBase_FACE, iMesh_ALL_TOPOLOGIES,
                    &ents, &ents_alloc, &ents_size, &err);
  }   
  else{
  iMesh_getEntities(impl, root_set,
                    iBase_REGION, iMesh_ALL_TOPOLOGIES,
                    &ents, &ents_alloc, &ents_size, &err);
  }
  ERRORR("Trouble getting entities for merge.", err);
  
  const double merge_tol =  1.0e-8;
  const int do_merge = 1;
  const int update_sets= 0; 
  iBase_TagHandle merge_tag = NULL;
  int num1, num2;

  iMesh_getNumOfType(impl, root_set, iBase_VERTEX, &num1, &err);
  ERRORR("Trouble getting number of entities before merge.", err);
  
  err = mm->merge_entities(ents, ents_size, merge_tol,
			   do_merge, update_sets,merge_tag);
  ERRORR("Trouble merging entities.", err);
  
  iMesh_getNumOfType(impl, root_set, iBase_VERTEX, &num2, &err);
  ERRORR("Trouble getting number of entities after merge.", err);

  std::cout << "Merged " << num1 - num2 << " vertices." << std::endl;

    // assign new global ids
  if (global_ids == true){
    std::cout << "Assigning global ids." << std::endl;
    MKUtils mu(impl);
    err = mu.assign_global_ids(root_set, 3, 1, true, false,
                             "GLOBAL_ID");
    ERRORR("Error assigning global ids.", err);
  }

   // export
  iMesh_save(impl, root_set, outfile.c_str(), NULL, &err, 
             strlen(outfile.c_str()), 0);
  ERRORR("Trouble writing output mesh.", err);
  std::cout << "Saved: "<< outfile.c_str() <<std::endl;

  delete cm;
  delete mm;
  
  iMesh_dtor(impl, &err);
  ERRORR("Destructor failed.", 1);

  // get the current date and time
  Timer.GetDateTime (szDateTime);
  std::cout << "\n\nEnding at : " << szDateTime;
 
  // compute the elapsed time
  std::cout << "Elapsed wall clock time: " << Timer.DiffTime ()
	    << " seconds or " << (Timer.DiffTime ())/60.0 << " mins\n";

  return iBase_SUCCESS;
}

int del_orig_mesh(std::vector<iBase_EntitySetHandle> &assys,
                  const bool back_mesh) 
{
  unsigned int end_idx = assys.size();
  if (back_mesh) end_idx--;
  
  iBase_EntityHandle *ents;
  int ents_alloc, ents_size;
  int err;
  
  for (unsigned int i = 0; i < end_idx; i++) {
    ents = NULL; ents_alloc = 0;
    iMesh_getEntities(impl, assys[i], 
                      iBase_ALL_TYPES, iMesh_ALL_TOPOLOGIES,
                      &ents, &ents_alloc, &ents_size, &err);
    ERRORR("Trouble getting entities in original set.", err);
    
    iMesh_deleteEntArr(impl, ents, ents_size, &err);
    ERRORR("Trouble deleting entities in original set.", err);

    free(ents);
  }
  
  return iBase_SUCCESS;
}

int copy_move_hex_assys(CopyMesh *cm,
			const int nrings, const int pack_type,
			const double pitch,
			const int symm,
			std::vector<std::string> &core_alias,
			std::vector<iBase_EntitySetHandle> &assys)
{
  double dx[3] = {0.0, 0.0, 0.0};
  double PI = acos(-1.0);
  iBase_EntityHandle *new_ents;
  int new_ents_alloc, new_ents_size;
  int err; 
  int i = 0, bd = 0;
  int assm_index;
  
  for (int n1 = 0; n1 < nrings; n1++) {
    if(n1%2==0){//check if n1 is even
      for (int n2 = 0; n2 < n1+1; n2++) {

	err = find_assm(i, assm_index);
	if (-1 == assm_index) {
	  i++;
	  if(n2 > n1/2)
	    ++bd; // index for assemblies below diagonal needs updatation
	  continue;
	}
	if(n2 <= n1/2){// before or equal to diagonal
	  dx[0] = n2 * pitch;
	  dx[1] = n1 * pitch * sin(PI/3.0);
	}
	else{//below the diagonal
	  dx[0] = (n1 + 1 + bd) * pitch / 2.0; 
	  dx[1] = (n1 - 1 - bd) * pitch * sin(PI/3.0);
	  ++bd;	  
	}	  
	new_ents = NULL;
	new_ents_alloc = 0;

	err = cm->copy_move_entities(assys[assm_index], dx, 
				     &new_ents, &new_ents_alloc, &new_ents_size,
				     false);
	ERRORR("Failed to copy_move entities.", 1);
	std::cout << "Copy/moved A: " << assm_index 
		  << " n1=" << n1 << ", n2=" << n2 <<" dX = " <<dx[0]<< " dY = " << dx[1] << std::endl;
	free(new_ents);
	i++;
      }
    }
    else{// n1 is odd
      for (int n2 = 0; n2 < n1; n2++) {
	err = find_assm(i, assm_index);
	if (-1 == assm_index){
	  i++;
	  if(n2 > (n1+1)/2)
	    ++bd; // index for assemblies below diagonal needs updatation
	  continue;

	}
	if(n2 <= (n1-1)/2){// before or equal to diagonal
	  dx[0] = (2 * n2 + 1) * pitch / 2.0;
	  dx[1] = n1 * pitch * sin(PI/3.0);

	}
	else{//below the diagonal 
	  dx[0] = (n1 + 1 + bd) * pitch / 2.0; 
	  if (bd == 0) // first n2 = 1 assembly
	    dx[1] = pitch * sin(PI/3.0);
	  dx[1] = (n1 - 1 - bd) * pitch * sin(PI/3.0);
	  ++bd;	
	}
	new_ents = NULL;
	new_ents_alloc = 0;

	err = cm->copy_move_entities(assys[assm_index], dx, 
				     &new_ents, &new_ents_alloc, &new_ents_size,
				     false);
	ERRORR("Failed to copy_move entities.", 1);
	std::cout << "Copy/moved A: " << assm_index 
		  << " n1=" << n1 << ", n2=" << n2 <<" dX = " <<dx[0]<< " dY = " << dx[1] << std::endl;
	free(new_ents);
	i++;
      }
    }
    bd = 0;
  }
  return iBase_SUCCESS;
}  

int copy_move_assys(CopyMesh *cm,
                    const int nrings, const int pack_type,
                    const double pitch,
                    const int symm,
                   	std::vector<std::string> &core_alias, 
                    std::vector<iBase_EntitySetHandle> &assys)
{
  double dx[3] = {0.0, 0.0, 0.0};
  double PI = acos(-1.0);
  iBase_EntityHandle *new_ents;
  int new_ents_alloc, new_ents_size;
  int err, assm_index; 
  int i = 0;
//   // get the copy/expand sets
//   int num_etags = 3, num_ctags = 1;
//   const char *etag_names[] = {"MATERIAL_SET", "DIRICHLET_SET", "NEUMANN_SET"};
//   const char *etag_vals[] = {NULL, NULL, NULL};
//   const char *ctag_names[] = {"GEOM_DIMENSION"};
//   const char *ctag_vals[]={(const char*)&set_DIM};

  for (int n1 = 0; n1 < nrings; n1++) {
    for (int n2 = n1; n2 < nrings; n2++) {
      err = find_assm(i,assm_index);
      if (-1 == assm_index) {
        i++;
        continue;
      }
      dx[1] = n1 * pitch * sin(PI/3.0);
      dx[0] = .5 * n1 * pitch + (n2 - n1) * pitch;
      new_ents = NULL;
      new_ents_alloc = 0;

      err = cm->copy_move_entities(assys[assm_index], dx, 
                                   &new_ents, &new_ents_alloc, &new_ents_size,
                                   false);
      ERRORR("Failed to copy_move entities.", 1);
      std::cout << "Copy/moved A: " << assm_index 
		<< " n1=" << n1 << ", n2=" << n2 <<" dX = " <<dx[0]<< " dY = " << dx[1] << std::endl;
      free(new_ents);
      i++;
    }
  }
  return iBase_SUCCESS;
}


int copy_move_sq_assys(CopyMesh *cm,
		       const int nrings, const int pack_type,
		       const double pitch,
		       const int symm,
		      	std::vector<std::string> &core_alias,
		       std::vector<iBase_EntitySetHandle> &assys) 
{
  double dx[3] = {0.0, 0.0, 0.0};
  iBase_EntityHandle *new_ents;
  int new_ents_alloc, new_ents_size;
  int err; 
  int i = 0, assm_index;

//   // get the copy/expand sets
//   int num_etags = 3, num_ctags = 1;
//   const char *etag_names[] = {"MATERIAL_SET", "DIRICHLET_SET", "NEUMANN_SET"};
//   const char *etag_vals[] = {NULL, NULL, NULL};
//   const char *ctag_names[] = {"GEOM_DIMENSION"};
//   const char *ctag_vals[]={(const char*)&set_DIM};
  for (int n1 = 0; n1 < nringsx; n1++) {
    dx[1] = n1 * pitchy;
    for (int n2 = 0; n2 < nringsy; n2++) {
      
      err = find_assm(i, assm_index);
      if (-1 == assm_index) {
        i++;
        continue;
      }
      dx[0] = n2 * pitchx;
      new_ents = NULL;
      new_ents_alloc = 0;

//       err = get_copy_expand_sets(cm,assys[assy_types[i]], etag_names, etag_vals, num_etags, CopyMesh::EXPAND);
//       ERRORR("Failed to add expand lists.", iBase_FAILURE);
//       err = get_copy_expand_sets(cm,assys[assy_types[i]], ctag_names, ctag_vals, num_ctags, CopyMesh::COPY);
//       ERRORR("Failed to add expand lists.", iBase_FAILURE);
//       err = extend_expand_sets(cm);
//       ERRORR("Failed to extend expand lists.", iBase_FAILURE);

      err = cm->copy_move_entities(assys[assm_index], dx, 
                                   &new_ents, &new_ents_alloc, &new_ents_size,
                                   false);
      ERRORR("Failed to copy_move entities.", 1);
      std::cout << "Copy/moved A: " << assm_index 
		<<  " n1=" << n1 << ", n2=" << n2  <<" dX = " <<dx[0]<< " dY = " << dx[1] << std::endl;
      free(new_ents);
      i++;
    
//       err = cm->tag_copied_sets(ctag_names, ctag_vals, 1);
//       ERRORR("Failed to tag copied sets.", iBase_FAILURE);
    
    }
  }
  return iBase_SUCCESS;
}

int get_copy_expand_sets(CopyMesh *cm,
                         iBase_EntitySetHandle orig_set, 
                         const char **tag_names, const char **tag_vals, 
                         int num_tags, int copy_or_expand) 
{
  int err = iBase_SUCCESS;

  for (int i = 0; i < num_tags; i++) {
    iBase_TagHandle tmp_tag;
    iMesh_getTagHandle(impl, tag_names[i], &tmp_tag, &err, strlen(tag_names[i]));
    ERRORR("Failed to get tag handle.", iBase_FAILURE);
    iBase_EntitySetHandle *tmp_sets = NULL;
    int tmp_alloc = 0, tmp_size;
    iMesh_getEntSetsByTagsRec(impl, orig_set, &tmp_tag, &tag_vals[i], 1, 0,
                              &tmp_sets, &tmp_alloc, &tmp_size, &err);
    ERRORR("Failure getting sets by tags.", err);
    if (0 != tmp_size) {
      err = cm->add_copy_expand_list(tmp_sets, tmp_size, copy_or_expand);
      ERRORR("Failed to add copy/expand lists.", iBase_FAILURE);
    }
    free(tmp_sets);
  }
  return err;
}

int extend_expand_sets(CopyMesh *cm) 
{
    // check expand sets for any contained sets which aren't already copy sets, 
    // and add them to the list
  int err;
  
  for (std::set<iBase_EntitySetHandle>::iterator sit = cm->expand_sets().begin();
       sit != cm->expand_sets().end(); sit++) {
    iBase_EntitySetHandle *sets = NULL;
    int sets_alloc = 0, sets_size;
    iMesh_getEntSets(impl, *sit, 1, &sets, &sets_alloc, &sets_size, &err);
    ERRORR("Failed to get contained sets.", err);
    if (sets_size) {
      err = cm->add_copy_expand_list(sets, sets_size, CopyMesh::COPY);
      ERRORR("Failed to add copy sets for expand extension.", err);
    }
    free(sets);
  }
  
  return iBase_SUCCESS;
}


int  prepareIO (int argc, char *argv[])
{
  bool bDone = false;
  do{
    if (2 == argc) {
      iname = argv[1];
      ifile = iname+".inp";
      outfile = iname+".h5m";
      mfile = iname + ".makefile";
    }
    else { //default case
      std::cerr << "Usage: " << argv[0] << " <input file> WITHOUT EXTENSION"<< std::endl;   
      iname = DEFAULT_TEST_FILE;
      ifile = iname+".inp";
      outfile = iname+".h5m";
      mfile = iname + ".makefile";
      std::cout <<"  No file specified.  Defaulting to: " << DEFAULT_TEST_FILE << std::endl;
    }
    // open the file
    file_input.open (ifile.c_str(), std::ios::in); 
    if (!file_input){
      std::cout << "Unable to open file" << std::endl;
      file_input.clear ();
    }
    else
      bDone = true; // file opened successfully
  } while (!bDone);
  // open Makefile-rgg
  do{
    make_file.open (mfile.c_str(), std::ios::out); 
    if (!make_file){
      std::cout << "Unable to open makefile for writing" << std::endl;
      make_file.clear ();
    }
    else
      bDone = true; // file opened successfully
  } while (!bDone);


  std::cout << "\nEntered input file name: " <<  ifile <<std::endl;
    
  return iBase_SUCCESS;
}

int find_assm(const int i, int &assm_index){
  int flag = 0;
  for(int j=0;j<nassys;j++)
    if(strcmp (core_alias[i].c_str(), assm_alias[j].c_str()) == 0){
      assm_index = j;
      flag = 1;
      break;
    }
  if(flag ==0)//nothing found return -1 or no assembly
    assm_index = -1;
  return iBase_SUCCESS;
}


int banner ()
{ 
  std::cout << '\n';
  std::cout << "\t\t--------------------------------------------------------------------" << '\n';
  std::cout << "\t\tProgram to Assemble Nuclear Reactor Assembly Meshes and Form a Core     " << '\n';
  std::cout << "\t\t\t\t\tArgonne National Laboratory" << '\n';
  std::cout << "\t\t\t\t\t        2009-2010         " << '\n';
  std::cout << "\t\t--------------------------------------------------------------------" << '\n';
 return iBase_SUCCESS;

}


int read_inputs ()
{
  CParser parse;
  int err;
  for(;;){
    if (!parse.ReadNextLine (file_input, linenumber, input_string, 
			     MAXCHARS, comment))
      ERRORR("Reading input file failed",1);

    if (input_string.substr(0,8) == "geometry"){
      std::istringstream formatString (input_string);
      formatString >> card >> geometry;
      if((strcmp(geometry.c_str(), "surface")==0)){
	  set_DIM = 2;
      }
    }
    if (input_string.substr(0,12) == "geometrytype"){
      std::istringstream formatString (input_string);
      formatString >> card >> geom_type;
      if(geom_type == "hexagonal"){

	// reading pitch info
	if (!parse.ReadNextLine (file_input, linenumber, input_string, 
				 MAXCHARS, comment))
	  ERRORR("Reading input file failed",1);
	if (input_string.substr(0,10) == "assemblies"){
	  std::istringstream formatString (input_string);
	  formatString >> card >> nassys >> pitch;
	}
	
	// reading file and alias names
	for(int i=1;i<=nassys;i++){
	  if (!parse.ReadNextLine (file_input, linenumber, input_string, 
				   MAXCHARS, comment))
	    ERRORR("Reading input file failed",1);	
	  std::istringstream formatString (input_string);
	  formatString >> meshfile >> mf_alias;
	  files.push_back(meshfile);
	  assm_alias.push_back(mf_alias);
	}
	
	// reading lattice 
	if (!parse.ReadNextLine (file_input, linenumber, input_string, 
				 MAXCHARS, comment))
	  ERRORR("Reading input file failed",1);
	if (input_string.substr(0,7) == "lattice"){
	  std::istringstream formatString (input_string);
	  formatString >> card >> nrings;
	  if(nrings % 2 == 0)
	    tot_assys = (nrings * (nrings))/2;
	  else
	    tot_assys = ((nrings * (nrings-1))/2) +  (nrings+1)/2;
	}

	// now reading the arrangement
	if (!parse.ReadNextLine (file_input, linenumber, input_string, 
				 MAXCHARS, comment))
	  ERRORR("Reading input file failed",1);	
	std::istringstream formatString (input_string);
	for(int i=1;i<=tot_assys;i++){
	  formatString >> temp_alias;
	  core_alias.push_back(temp_alias);
	}		
      }

      else if(geom_type == "cartesian"){

	// reading pitch info
	if (!parse.ReadNextLine (file_input, linenumber, input_string, 
				 MAXCHARS, comment))
	  ERRORR("Reading input file failed",1);
	if (input_string.substr(0,10) == "assemblies"){
	  std::istringstream formatString (input_string);
	  formatString >> card >> nassys >> pitchx >> pitchy;
	}
	
	// reading file and alias names
	for(int i=1;i<=nassys;i++){
	  if (!parse.ReadNextLine (file_input, linenumber, input_string, 
				   MAXCHARS, comment))
	    ERRORR("Reading input file failed",1);	
	  std::istringstream formatString (input_string);
	  formatString >> meshfile >> mf_alias;
	  files.push_back(meshfile);
	  assm_alias.push_back(mf_alias);
	}
	
	// reading lattice 
	if (!parse.ReadNextLine (file_input, linenumber, input_string, 
				 MAXCHARS, comment))
	  ERRORR("Reading input file failed",1);
	if (input_string.substr(0,7) == "lattice"){
	  std::istringstream formatString (input_string);
	  formatString >> card >> nringsx >> nringsy;
	  tot_assys = nringsx * nringsy;
	}

	// now reading the arrangement
	if (!parse.ReadNextLine (file_input, linenumber, input_string, 
				 MAXCHARS, comment))
	  ERRORR("Reading input file failed",1);	
	std::istringstream formatString (input_string);
	for(int i=1;i<=tot_assys;i++){
	  formatString >> temp_alias;
	  core_alias.push_back(temp_alias);
	}		
      }

      else if (geom_type =="default"){
	// reading pitch info
	if (!parse.ReadNextLine (file_input, linenumber, input_string, 
				 MAXCHARS, comment))
	  ERRORR("Reading input file failed",1);
	if (input_string.substr(0,10) == "assemblies"){
	  std::istringstream formatString (input_string);
	  formatString >> card >> nassys >> pitch;
	}
	
	// reading file and alias names
	for(int i=1;i<=nassys;i++){
	  if (!parse.ReadNextLine (file_input, linenumber, input_string, 
				   MAXCHARS, comment))
	    ERRORR("Reading input file failed",1);	
	  std::istringstream formatString (input_string);
	  formatString >> meshfile >> mf_alias;
	  files.push_back(meshfile);
	  assm_alias.push_back(mf_alias);
	}
	
	// reading lattice 
	if (!parse.ReadNextLine (file_input, linenumber, input_string, 
				 MAXCHARS, comment))
	  ERRORR("Reading input file failed",1);
	if (input_string.substr(0,7) == "lattice"){
	  std::istringstream formatString (input_string);
	  formatString >> card >> nrings;
	  tot_assys = (nrings * (nrings+1))/2;
	}

	// now reading the arrangement
	if (!parse.ReadNextLine (file_input, linenumber, input_string, 
				 MAXCHARS, comment))
	  ERRORR("Reading input file failed",1);	
	std::istringstream formatString (input_string);
	for(int i=1;i<=tot_assys;i++){
	  formatString >> temp_alias;
	  core_alias.push_back(temp_alias);
	}		
      }
      else{
	ERRORR("Invalid geometry type",1);
      }
    }
    // background mesh
    if (input_string.substr(0,10) == "background"){
      std::istringstream formatString (input_string);
      formatString >> card >> back_meshfile;
      files.push_back(back_meshfile);	  
      back_mesh = true;
    }

    // symmetry
    if (input_string.substr(0,8) == "symmetry"){
      std::istringstream formatString (input_string);
      formatString >> card >> symm;	
    }
    // breaking condition
    if(input_string.substr(0,3) == "end"){
      std::istringstream formatstring (input_string);
      break;
    }
  }
  return iBase_SUCCESS;
}


int write_makefile(){
  std::string name;
  std::vector<std::string> f_no_ext, f_sat, f_inp, f_jou, f_injou;
  make_file << "##" << std::endl;
  make_file << "## This makefile is automatically generated by coregen program" << std::endl;
  make_file << "##" << std::endl;
  make_file << "\nCUBIT = /homes/jain/downloads/cubit-10.2-amd64/bin/cubit --no-gui --batch\n" << std::endl;
  make_file << "coregen = /homes/jain/MeshKit/coregen\n" << std::endl;
  make_file << "assygen = /homes/jain/MeshKit-acis/rgg/assygen\n" << std::endl;

  make_file << "MESH_FILES = " ;
  for(unsigned int i=0; i<files.size(); i++)
    make_file << files[i] << "  ";

  
  // get file names without extension
  for(unsigned int i=0; i<files.size(); i++){
    int loc = files[i].find_first_of(".");
    f_no_ext.push_back(files[i].substr(0,loc));
  }
  
  make_file << "\n\nGEOM_FILES = ";
  for(unsigned int i=0; i<files.size(); i++){
    name = f_no_ext[i] + ".sat";
    f_sat.push_back(name);
    make_file << name << "  ";
    name = "";
  }
  
  make_file << "\n\nJOU_FILES = ";
  for(unsigned int i=0; i<files.size(); i++){
    name = f_no_ext[i] + ".jou";
    f_jou.push_back(name);
    make_file << name << "  ";
    name = "";
  }

  make_file << "\n\nINJOU_FILES = ";
  for(unsigned int i=0; i<files.size(); i++){
    name = f_no_ext[i] + ".template.jou";
    f_injou.push_back(name);
    make_file << name << "  ";
    name = "";
  }

  make_file << "\n\nASSYGEN_FILES = ";
  for(unsigned int i=0; i<files.size(); i++){
    name = f_no_ext[i] + ".inp";
    f_inp.push_back(name);
    make_file << name << "  ";
    name = "";
  }

  make_file << "\n\n" << outfile << " : ${MESH_FILES}" << std::endl;
  make_file << "\t" << "${coregen} " << iname << std::endl;
  for(unsigned int i=0; i<files.size(); i++){
    make_file << files[i] << " : " << f_sat[i] << "  " << f_jou[i] << "  " << f_injou[i] << std::endl;
    make_file << "\t" << "${CUBIT} " << f_jou[i] <<"\n" << std::endl;

    make_file << f_sat[i] << " : " << f_inp[i] << std::endl;
    make_file << "\t" << "${assygen} " << f_no_ext[i] << "\n" << std::endl;
  }

  make_file.close();
  return 0;
}
