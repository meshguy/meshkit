/*
 * SIXTH program
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

iMesh_Instance impl;
iBase_EntitySetHandle root_set;

const int UNITCELL_DUCT = 0, ASSY_TYPES = 1;
int set_DIM = 3;
int get_copy_expand_sets(CopyMesh *cm,
                         iBase_EntitySetHandle orig_set, 
                         const char **tag_names, const char **tag_vals, 
                         int num_tags, int copy_or_expand);

int extend_expand_sets(CopyMesh *cm);

int copy_move_assys(CopyMesh *cm,
                    const int nrings, const int pack_type,
                    const double pitch,
                    const int symm,
                    std::vector<int> &assy_types,
                    std::vector<iBase_EntitySetHandle> &assys);

int copy_move_sq_assys(CopyMesh *cm,
		       const int nrings, const int pack_type,
		       const double pitch,
		       const int symm,
		       std::vector<int> &assy_types,
		       std::vector<iBase_EntitySetHandle> &assys);

int copy_move_hex_assys(CopyMesh *cm,
		       const int nrings, const int pack_type,
		       const double pitch,
		       const int symm,
		       std::vector<int> &assy_types,
		       std::vector<iBase_EntitySetHandle> &assys);

int read_input(int &nrings, int &pack_type, double &pitch, int &symm,
               bool &back_mesh,
               std::vector<std::string> &files, 
               std::vector<int> &assy_types,
               std::string &outfile, bool &global_ids);


int read_hex_input(int &nrings, int &pack_type, double &pitch, int &symm,
		   bool &back_mesh,
		   std::vector<std::string> &files, 
		   std::vector<int> &assy_types,
		   std::string &outfile, bool &global_ids);

int read_sq_input(int &nrings, int &pack_type, double &pitch, int &symm,
               bool &back_mesh,
               std::vector<std::string> &files, 
               std::vector<int> &assy_types,
               std::string &outfile, bool &global_ids);

int read_input_defaults(int &nrings, int &pack_type, double &pitch, int &symm,
               bool &back_mesh,
               std::vector<std::string> &files, 
               std::vector<int> &assy_types,
			std::string &outfile, bool &global_ids);

int del_orig_mesh(std::vector<iBase_EntitySetHandle> &assys,
                  const bool back_mesh);


int main(int argc, char **argv) 
{
  // set the dimension based on input mesh - 2D meshes or 3D mesh??
  char* input;
  int test_flag;
    if (1 == argc) {
    test_flag = 1;
  }
    // make a mesh instance
  int err;
  iMesh_newMesh("MOAB", &impl, &err, 4);
  ERRORR("Failed to create instance.", 1);
  
  iMesh_getRootSet(impl, &root_set, &err);

  CopyMesh *cm = new CopyMesh(impl);
  MergeMesh *mm = new MergeMesh(impl);

  int nrings, pack_type, symm;
  double pitch;
  std::vector<std::string> files;
  std::vector<int> assy_types;
  bool back_mesh;
  std::string outfile;
  bool global_ids;

  // reading user input
  if (test_flag == 1){
    input = (char*) "default";
    err = read_input_defaults(nrings, pack_type, pitch, symm, back_mesh,
			      files, assy_types, outfile, global_ids);
  ERRORR("Couldn't parse input.", err);
  }
  else{
    input = argv[1];
    if(!strcmp(input, "hexagonal")){
      err = read_hex_input(nrings, pack_type, pitch, symm, back_mesh,
		       files, assy_types, outfile, global_ids);
      ERRORR("Couldn't parse input.", err);
    }
    else if(!strcmp(input,"square")){
      err = read_sq_input(nrings, pack_type, pitch, symm, back_mesh,
		       files, assy_types, outfile, global_ids);
      ERRORR("Couldn't parse input.", err);
    }
  }

  std::vector<iBase_EntitySetHandle> assys;
  iBase_EntitySetHandle orig_set;
  for (unsigned int i = 0; i < files.size(); i++) {
    iMesh_createEntSet(impl, 0, &orig_set, &err);
    ERRORR("Couldn't create file set.", err);
  
      // read the file
    iMesh_load(impl, orig_set, files[i].c_str(), NULL, &err, 
               strlen(files[i].c_str()), 0);
    ERRORR("Couldn't read mesh file.", err);

    assys.push_back(orig_set);
  }
  std::cout << "Loaded mesh files." << std::endl;

//   // get the copy/expand sets
//   int num_etags = 3, num_ctags = 1;
//   const char *etag_names[] = {"MATERIAL_SET", "DIRICHLET_SET", "NEUMANN_SET"};
//   const char *etag_vals[] = {NULL, NULL, NULL};
//   const char *ctag_names[] = {"GEOM_DIMENSION"};
//   const char *ctag_vals[]={(const char*)&set_DIM};
//   for (unsigned int i = 0; i < files.size(); i++) {
//     err = get_copy_expand_sets(cms[i],assys[assy_types[i]], etag_names, etag_vals, num_etags, CopyMesh::EXPAND);
//     ERRORR("Failed to add expand lists.", iBase_FAILURE);
//     err = get_copy_expand_sets(cms[i],assys[assy_types[i]], ctag_names, ctag_vals, num_ctags, CopyMesh::COPY);
//     ERRORR("Failed to add expand lists.", iBase_FAILURE);
//     err = extend_expand_sets(cms[i]);
//     ERRORR("Failed to extend expand lists.", iBase_FAILURE);
//   }

     // move the assys based on the geometry type
  if(!strcmp(input, "hexagonal")){
    err = copy_move_hex_assys(cm, nrings, pack_type, pitch, symm, 
			  assy_types, assys);
    ERRORR("Failed in copy/move step.", err);
  }
  else if(!strcmp(input,"square")){
    err = copy_move_sq_assys(cm, nrings, pack_type, pitch, symm, 
			  assy_types, assys);
    ERRORR("Failed in copy/move step.", err);
  } 
  else if(!strcmp(input,"default")){
    err = copy_move_assys(cm, nrings, pack_type, pitch, symm, 
			  assy_types, assys);
    ERRORR("Failed in copy/move step.", err);
  }
  
  // delete all original meshes
  err = del_orig_mesh(assys, back_mesh);
  ERRORR("Failed in delete step.", err);

  //get entities for merging
  iBase_EntityHandle *ents = NULL; 
  int ents_alloc = 0, ents_size;

  if(set_DIM ==2){
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
			std::vector<int> &assy_types,
			std::vector<iBase_EntitySetHandle> &assys)
{
  double dx[3] = {0.0, 0.0, 0.0};
  double PI = acos(-1.0);
  iBase_EntityHandle *new_ents;
  int new_ents_alloc, new_ents_size;
  int err; 
  int i = 0, bd = 0;


  for (int n1 = 0; n1 < nrings; n1++) {
    if(n1%2==0){//check if n1 is even
      for (int n2 = 0; n2 < n1+1; n2++) {
	if (-1 == assy_types[i]) {
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

	err = cm->copy_move_entities(assys[assy_types[i]], dx, 
				     &new_ents, &new_ents_alloc, &new_ents_size,
				     false);
	ERRORR("Failed to copy_move entities.", 1);
	std::cout << "Copy/moved n1=" << n1 << ", n2=" << n2 <<" dX = " <<dx[0]<< " dY = " << dx[1] << std::endl;
	free(new_ents);
	i++;
      }
    }
    else{// n1 is odd
      for (int n2 = 0; n2 < n1; n2++) {
	if (-1 == assy_types[i]){
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

	err = cm->copy_move_entities(assys[assy_types[i]], dx, 
				     &new_ents, &new_ents_alloc, &new_ents_size,
				     false);
	ERRORR("Failed to copy_move entities.", 1);
	std::cout << "Copy/moved n1=" << n1 << ", n2=" << n2 <<" dX = " <<dx[0]<< " dY = " << dx[1] << std::endl;
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
                    std::vector<int> &assy_types,
                    std::vector<iBase_EntitySetHandle> &assys)
{
  double dx[3] = {0.0, 0.0, 0.0};
  double PI = acos(-1.0);
  iBase_EntityHandle *new_ents;
  int new_ents_alloc, new_ents_size;
  int err; 
  int i = 0;
  // get the copy/expand sets
  int num_etags = 3, num_ctags = 1;
  const char *etag_names[] = {"MATERIAL_SET", "DIRICHLET_SET", "NEUMANN_SET"};
  const char *etag_vals[] = {NULL, NULL, NULL};
  const char *ctag_names[] = {"GEOM_DIMENSION"};
  const char *ctag_vals[]={(const char*)&set_DIM};






  err = get_copy_expand_sets(cm,assys[assy_types[i]], etag_names, etag_vals, num_etags, CopyMesh::EXPAND);
  ERRORR("Failed to add expand lists.", iBase_FAILURE);
  err = get_copy_expand_sets(cm,assys[assy_types[i]], ctag_names, ctag_vals, num_ctags, CopyMesh::COPY);
  ERRORR("Failed to add expand lists.", iBase_FAILURE);
  err = extend_expand_sets(cm);
  ERRORR("Failed to extend expand lists.", iBase_FAILURE);

  for (int n1 = 0; n1 < nrings; n1++) {
    for (int n2 = n1; n2 < nrings; n2++) {
      if (-1 == assy_types[i]) {
        i++;
        continue;
      }
      dx[1] = n1 * pitch * sin(PI/3.0);
      dx[0] = .5 * n1 * pitch + (n2 - n1) * pitch;
      new_ents = NULL;
      new_ents_alloc = 0;

      err = cm->copy_move_entities(assys[assy_types[i]], dx, 
                                   &new_ents, &new_ents_alloc, &new_ents_size,
                                   false);
      ERRORR("Failed to copy_move entities.", 1);
      std::cout << "Copy/moved n1=" << n1 << ", n2=" << n2 <<" dX = " <<dx[0]<< " dY = " << dx[1] << std::endl;
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
		       std::vector<int> &assy_types,
		       std::vector<iBase_EntitySetHandle> &assys) 
{
  double dx[3] = {0.0, 0.0, 0.0};
  iBase_EntityHandle *new_ents;
  int new_ents_alloc, new_ents_size;
  int err; 
  int i = 0;

  // get the copy/expand sets
  int num_etags = 3, num_ctags = 1;
  const char *etag_names[] = {"MATERIAL_SET", "DIRICHLET_SET", "NEUMANN_SET"};
  const char *etag_vals[] = {NULL, NULL, NULL};
  const char *ctag_names[] = {"GEOM_DIMENSION"};
  const char *ctag_vals[]={(const char*)&set_DIM};
  for (int n1 = 0; n1 < nrings; n1++) {
    dx[1] = n1 * pitch;
    for (int n2 = 0; n2 < nrings; n2++) {
      if (-1 == assy_types[i]) {
        i++;
        continue;
      }
      dx[0] = n2 * pitch;
      new_ents = NULL;
      new_ents_alloc = 0;

      err = get_copy_expand_sets(cm,assys[assy_types[i]], etag_names, etag_vals, num_etags, CopyMesh::EXPAND);
      ERRORR("Failed to add expand lists.", iBase_FAILURE);
      err = get_copy_expand_sets(cm,assys[assy_types[i]], ctag_names, ctag_vals, num_ctags, CopyMesh::COPY);
      ERRORR("Failed to add expand lists.", iBase_FAILURE);
      err = extend_expand_sets(cm);
      ERRORR("Failed to extend expand lists.", iBase_FAILURE);

      err = cm->copy_move_entities(assys[assy_types[i]], dx, 
                                   &new_ents, &new_ents_alloc, &new_ents_size,
                                   false);
      ERRORR("Failed to copy_move entities.", 1);
      std::cout << "Copy/moved n1=" << n1 << ", n2=" << n2  <<" dX = " <<dx[0]<< " dY = " << dx[1] << std::endl;
      free(new_ents);
      i++;
    
//       err = cm->tag_copied_sets(ctag_names, ctag_vals, 1);
//       ERRORR("Failed to tag copied sets.", iBase_FAILURE);
    
    }
  }
  return iBase_SUCCESS;
}


int read_hex_input(int &nrings, int &pack_type, double &pitch, int &symm,
               bool &back_mesh,
               std::vector<std::string> &files, 
               std::vector<int> &assy_types,
               std::string &outfile, bool &global_ids) 
{
  int tot_assys;
  std::cout << "Nrings: ";
  std::cin >> nrings;
  std::cout << "Packing type: ";
  std::cin >> pack_type;
  std::cout << "Pitch: ";
  std::cin >> pitch;
  std::cout << "1/Symmetry: ";
  std::cin >> symm;

  if (UNITCELL_DUCT == pack_type) {
    std::string filename;
    std::cout << "Unit cell file: "; std::cin >> filename; files.push_back(filename);
    std::cout << "Duct file: "; std::cin >> filename; files.push_back(filename);
    back_mesh = true;
  }
  else if (ASSY_TYPES == pack_type) {
      // read n(n+1)/2 assy types, then file names
    if(nrings % 2 == 0)
      tot_assys = (nrings * (nrings))/2;
    else
      tot_assys = ((nrings * (nrings-1))/2) +  (nrings+1)/2;
    int dum;
    for (int i = 0; i < tot_assys; i++) {
      std::cin >> dum;
      assy_types.push_back(dum);
    }
    
      // get the number of assy types then read that many filenames
    std::set<int> all_assys;
    for (std::vector<int>::iterator vit = assy_types.begin();
         vit != assy_types.end(); vit++) all_assys.insert(*vit);
    all_assys.erase(-1);
    char filename[1024];
    for (unsigned int i = 0; i < all_assys.size(); i++) {
      std::cout << "Assy " << i << " file: "; 
      std::cin >> filename;
      files.push_back(std::string(filename));
    }
    std::cout << "Background assembly (y/n)?" << std::endl;
    std::cin >> filename;
    if (filename[0] == 'y' || filename[0] == 'Y') {
      std::cin >> filename;
      files.push_back(std::string(filename));
      back_mesh = true;
    }
    else
      back_mesh = false;

    std::cout << "New global ids (y/n)?" << std::endl;
    std::cin >> filename;
    if (filename[0] == 'y' || filename[0] == 'Y') {
      std::cin >> filename;
      global_ids = true;
    }
    else
      global_ids = false;
  }

  char filename[1024];
  std::cout << "Output file: " << std::endl;
  std::cin >> filename;
  outfile = filename;

  return iBase_SUCCESS;
}



int read_input(int &nrings, int &pack_type, double &pitch, int &symm,
               bool &back_mesh,
               std::vector<std::string> &files, 
               std::vector<int> &assy_types,
               std::string &outfile, bool &global_ids) 
{
  std::cout << "Nrings: ";
  std::cin >> nrings;
  std::cout << "Packing type: ";
  std::cin >> pack_type;
  std::cout << "Pitch: ";
  std::cin >> pitch;
  std::cout << "1/Symmetry: ";
  std::cin >> symm;

  if (UNITCELL_DUCT == pack_type) {
    std::string filename;
    std::cout << "Unit cell file: "; std::cin >> filename; files.push_back(filename);
    std::cout << "Duct file: "; std::cin >> filename; files.push_back(filename);
    back_mesh = true;
  }
  else if (ASSY_TYPES == pack_type) {
      // read n(n+1)/2 assy types, then file names
    int tot_assys = (nrings * (nrings+1))/2;
    int dum;
    for (int i = 0; i < tot_assys; i++) {
      std::cin >> dum;
      assy_types.push_back(dum);
    }
    
      // get the number of assy types then read that many filenames
    std::set<int> all_assys;
    for (std::vector<int>::iterator vit = assy_types.begin();
         vit != assy_types.end(); vit++) all_assys.insert(*vit);
    all_assys.erase(-1);
    char filename[1024];
    for (unsigned int i = 0; i < all_assys.size(); i++) {
      std::cout << "Assy " << i << " file: "; 
      std::cin >> filename;
      files.push_back(std::string(filename));
    }
    std::cout << "Background assembly (y/n)?" << std::endl;
    std::cin >> filename;
    if (filename[0] == 'y' || filename[0] == 'Y') {
      std::cin >> filename;
      files.push_back(std::string(filename));
      back_mesh = true;
    }
    else
      back_mesh = false;

    std::cout << "New global ids (y/n)?" << std::endl;
    std::cin >> filename;
    if (filename[0] == 'y' || filename[0] == 'Y') {
      std::cin >> filename;
      global_ids = true;
    }
    else
      global_ids = false;
  }

  char filename[1024];
  std::cout << "Output file: " << std::endl;
  std::cin >> filename;
  outfile = filename;

  return iBase_SUCCESS;
}


int read_input_defaults(int &nrings, int &pack_type, double &pitch, int &symm,
               bool &back_mesh,
               std::vector<std::string> &files, 
               std::vector<int> &assy_types,
	       std::string &outfile, bool &global_ids) 
{
  std::cout << "Using defaults for running sixth program \n  Usage: sixth <any_character> for user defined inputs";
  nrings = 2;
  std::cout << "\n\n----Inputs to Sixth Program----\n\nNrings: " << nrings;
  pack_type = 1;
  std::cout << "\nPacking type: 1 (for Assembly) 0 (for Unit Cell Duct)"; 
  std::cout << "\nUsing packing type: " << pack_type;
  pitch = 15.685;
  std::cout << "\nPitch: "<< pitch;
  symm = 6;
  std::cout << "\n1/Symmetry: " << symm;


  // ASSY_TYPES is setup as default case
  if (UNITCELL_DUCT == pack_type) {
    std::string filename;
    std::cout << "\nUnit cell file: "; std::cin >> filename; files.push_back(filename);
    std::cout << "\nDuct file: "; std::cin >> filename; files.push_back(filename);
    back_mesh = true;
  }
  else if (ASSY_TYPES == pack_type) {
    std::cout << "\nNo. of assembly types = 2"; //both types are pin1.cub
    //manually enter nrings(nrings+1)/2 assy types
    //  for(int test=1;test<=15;test++){
    assy_types.push_back(0);
    assy_types.push_back(1);
    assy_types.push_back(2);
     //    }
    
    // get the number of assy types then read that many filenames
    std::set<int> all_assys;
    for (std::vector<int>::iterator vit = assy_types.begin();
	 vit != assy_types.end(); vit++) all_assys.insert(*vit);
    all_assys.erase(-1);
   
    //char filename[] = "fuel_assy.cub";
    char filename[] = "pin1.cub";
    for (unsigned int i = 0; i < all_assys.size(); i++) {
      std::cout << "\nAssy " << i  << ": "<< filename;
      files.push_back(std::string(filename));
    }
    char bfilename[] = "y"; 
    std::cout << "\nBackground assembly (y/n)?: " << bfilename << std::endl;
     
    if (bfilename[0] == 'y' || bfilename[0] == 'Y') {
      char bfile[] = "test_sodium_all.cub";
      //char bfile[] = "sodium_all.cub";
      std::cout << "Defaulting to file: " << bfile;

      files.push_back(std::string(bfile));
      back_mesh = true;
    }
    else
      back_mesh = false;
  }
  char g_id[]= "y";
  std::cout << "\nNew global ids (y/n)? " << g_id;
  global_ids = true;
  char filename[] = "sixth_test.h5m";
  std::cout << "\nOutput file: "<< filename << std::endl;
  outfile = filename;

  std::cout << "\n------------------------------\n\n";

  return iBase_SUCCESS;
}

int read_sq_input(int &nrings, int &pack_type, double &pitch, int &symm,
               bool &back_mesh,
               std::vector<std::string> &files, 
               std::vector<int> &assy_types,
	       std::string &outfile, bool &global_ids) 
{ std::cout << "Number of assemblies in X/Y (assuming square): ";
  std::cin >> nrings;
  std::cout << "Packing type: ";
  std::cin >> pack_type;
  std::cout << "Pitch: ";
  std::cin >> pitch;
//   std::cout << "1/Symmetry: ";
//   std::cin >> symm;

  if (UNITCELL_DUCT == pack_type) {
    std::string filename;
    std::cout << "Unit cell file: "; std::cin >> filename; files.push_back(filename);
    std::cout << "Duct file: "; std::cin >> filename; files.push_back(filename);
    back_mesh = true;
  }
  else if (ASSY_TYPES == pack_type) {
      // read n(n+1)/2 assy types, then file names
    int tot_assys = nrings*nrings;
    int dum;
    for (int i = 0; i < tot_assys; i++) {
      std::cin >> dum;
      assy_types.push_back(dum);
    }
    
    // get the number of assy types then read that many filenames
    std::set<int> all_assys;
    for (std::vector<int>::iterator vit = assy_types.begin();
         vit != assy_types.end(); vit++) all_assys.insert(*vit);
    all_assys.erase(-1);
    char filename[1024];
    for (unsigned int i = 0; i < all_assys.size(); i++) {
      std::cout << "Assy " << i << " file: "; 
      std::cin >> filename;
      files.push_back(std::string(filename));
    }
    std::cout << "Background assembly (y/n)?" << std::endl;
    std::cin >> filename;
    if (filename[0] == 'y' || filename[0] == 'Y') {
      std::cin >> filename;
      files.push_back(std::string(filename));
      back_mesh = true;
    }
    else
      back_mesh = false;

    std::cout << "New global ids (y/n)?" << std::endl;
    std::cin >> filename;
    if (filename[0] == 'y' || filename[0] == 'Y') {
      std::cin >> filename;
      global_ids = true;
    }
    else
      global_ids = false;
  }

  char filename[1024];
  std::cout << "Output file: " << std::endl;
  std::cin >> filename;
  outfile = filename;

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


