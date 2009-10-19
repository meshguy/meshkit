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
               std::string &outfile);

int del_orig_mesh(std::vector<iBase_EntitySetHandle> &assys,
                  const bool back_mesh);

int main(int argc, char **argv) 
{
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
  err = read_input(nrings, pack_type, pitch, symm, back_mesh,
                   files, assy_types, outfile);
  ERRORR("Couldn't parse input.", err);

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

  int num1, num2;
  iMesh_getNumOfType(impl, root_set, iBase_VERTEX, &num1, &err);
  ERRORR("Trouble getting number of entities before merge.", err);
  
//  err = mm->merge_entities(ents, ents_size, 1.0e-7);
//  ERRORR("Trouble merging entities.", err);
  
  iMesh_getNumOfType(impl, root_set, iBase_VERTEX, &num2, &err);
  ERRORR("Trouble getting number of entities after merge.", err);

  std::cout << "Merged " << num1 - num2 << " vertices." << std::endl;

  iMesh_save(impl, root_set, outfile.c_str(), NULL, &err, 
             strlen(outfile.c_str()), 0);
  ERRORR("Trouble writing output mesh.", err);
  
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
               std::string &outfile) 
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
  }

  char filename[1024];
  std::cout << "Output file: " << std::endl;
  scanf("%s", filename); 
  outfile = filename;

  return iBase_SUCCESS;
}
