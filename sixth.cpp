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

int read_input(int &nrings, int &pack_type, double &pitch, int &symm,
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
  if (test_flag == 1){
    err = read_input_defaults(nrings, pack_type, pitch, symm, back_mesh,
			      files, assy_types, outfile, global_ids);
  ERRORR("Couldn't parse input.", err);
  }
  else{
    err = read_input(nrings, pack_type, pitch, symm, back_mesh,
		     files, assy_types, outfile, global_ids);
  ERRORR("Couldn't parse input.", err);
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

  err = copy_move_assys(cm, nrings, pack_type, pitch, symm, 
                        assy_types, assys);
  ERRORR("Failed in copy/move step.", err);

  err = del_orig_mesh(assys, back_mesh);
  ERRORR("Failed in delete step.", err);

  iBase_EntityHandle *ents = NULL; 
  int ents_alloc = 0, ents_size;
  iMesh_getEntities(impl, root_set,
                    iBase_REGION, iMesh_ALL_TOPOLOGIES,
                    &ents, &ents_alloc, &ents_size, &err);
  ERRORR("Trouble getting entities for merge.", err);
  
  const double merge_tol =  1.0e-8;
  const int do_merge = 0;
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
  for (int n1 = 0; n1 < nrings; n1++) {
    dx[1] = n1 * pitch * sin(PI/3.0);
    for (int n2 = n1; n2 < nrings; n2++) {
      if (-1 == assy_types[i]) {
        i++;
        continue;
      }
      dx[0] = .5 * n1 * pitch + (n2 - n1) * pitch;
      new_ents = NULL;
      new_ents_alloc = 0;
      err = cm->copy_move_entities(assys[assy_types[i]], dx, 
                                   &new_ents, &new_ents_alloc, &new_ents_size,
                                   false);
      ERRORR("Failed to copy_move entities.", 1);
      std::cout << "Copy/moved n1=" << n1 << ", n2=" << n2 << std::endl;
      free(new_ents);
      i++;
    }
  }
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
      scanf("%s", filename); 
      files.push_back(std::string(filename));
    }
    std::cout << "Background assembly (y/n)?" << std::endl;
    scanf("%s", filename);
    if (filename[0] == 'y' || filename[0] == 'Y') {
      scanf("%s", filename); 
      files.push_back(std::string(filename));
      back_mesh = true;
    }
    else
      back_mesh = false;

    std::cout << "New global ids (y/n)?" << std::endl;
    scanf("%s", filename);
    if (filename[0] == 'y' || filename[0] == 'Y') {
      scanf("%s", filename); 
      global_ids = true;
    }
    else
      global_ids = false;
  }

  char filename[1024];
  std::cout << "Output file: " << std::endl;
  scanf("%s", filename); 
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
  nrings = 3;
  std::cout << "\n\n----Inputs to Sixth Program----\n\nNrings: " << nrings;
  pack_type = 1;
  std::cout << "\nPacking type: 1 (for Assembly) 0 (for Unit Cell Duct)"; 
  std::cout << "\nUsing packing type: " << pack_type;
  pitch = 14.685;
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
     assy_types.push_back(0);
     assy_types.push_back(1);
     assy_types.push_back(0);
    
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
  char g_id[]= "n";
  std::cout << "\nNew global ids (y/n)? " << g_id;
  global_ids = false;
  char filename[] = "sixth_test.h5m";
  std::cout << "\nOutput file: "<< filename << std::endl;
  outfile = filename;

  std::cout << "\n------------------------------\n\n";

  return iBase_SUCCESS;
}


