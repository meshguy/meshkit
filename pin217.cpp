#include <cfloat>
#include <math.h>
#include "MKUtils.hpp"
#include "MergeMesh.hpp"
#include "CopyMesh.hpp"

iMesh_Instance impl;
iBase_EntityHandle root_set;

int get_copy_expand_sets(CopyMesh *cm,
                         iBase_EntitySetHandle orig_set, 
                         const char **tag_names, const char **tag_vals, 
                         int num_tags, int copy_or_expand);

int extend_expand_sets(CopyMesh *cm);

int main(int argc, char **argv) 
{
    // make a mesh instance
  int err;
  iMesh_newMesh("MOAB", &impl, &err, 4);
  ERRORR("Failed to create instance.", 1);
  
  iMesh_getRootSet(impl, &root_set, &err);

  CopyMesh *cm = new CopyMesh(impl);
  MergeMesh *mm = new MergeMesh(impl);

  if (1 == argc) {
    std::cerr << "Usage: " << argv[0] << " <meshfile>" << std::endl;
    return 1;
  }

  std::vector<std::string> ctags, etags, utags;
  std::vector<char*> cvals, evals;
  std::string infile, outfile;
  double dx[3] = {0.0, 0.0, 0.0};
  double ASSY_PITCH=14.685;
  double PI = acos(-1.0);
  double Y_SEP=.5*ASSY_PITCH/sin(PI/3.0) + 
      .5*ASSY_PITCH*sin(PI/6.0)/sin(PI/3.0);

  int NRINGS = 8;
  bool sixth_core[] = {
      false, true, false, false, false, false, false, false,
      false, false, false, false, false, false, false, false,
      false, false, false, false, false, false, false, false,
      false, false, false, false, false, false, false, false,
      false, false, false, false, false, false, false, false,
      false, false, false, false, false, false, false, false,
      false, false, false, false, false, false, false, false,
      false, false, false, false, false, false, false, false 
  };
/*  bool sixth_core[] = {
      false, true, true, true, true, true, true, true, 
      true, true, true, true, true, true, true, false,
      true, true, true, true, true, true, true, false,
      true, true, true, true, true, true, false, false,
      true, true, true, true, true, false, false, false,
      true, true, true, true, false, false, false, false,
      true, true, true, false, false, false, false, false,
      true, false, false, false, false, false, false, false 
  };
*/
  iBase_EntitySetHandle orig_set;
  iMesh_createEntSet(impl, 0, &orig_set, &err);
  ERRORR("Couldn't create orig entity set.", err);
  
    // read the file
  iMesh_load(impl, orig_set, argv[1], NULL, &err, strlen(argv[1]), 0);
  ERRORR("Couldn't read mesh file.", err);

  std::cout << "Loaded initial mesh." << std::endl;

    // get the entities we want to copy, basically all the non-set entities
  iBase_EntityHandle *orig_ents = NULL;
  int orig_ents_alloc = 0, orig_ents_size;
  iMesh_getEntities(impl, orig_set, iBase_ALL_TYPES, iMesh_ALL_TOPOLOGIES,
                    &orig_ents, &orig_ents_alloc, &orig_ents_size, &err);
  ERRORR("Failed to get any entities from original set.", iBase_FAILURE);
  
    // get the copy/expand sets
  int three = 3;
  int num_etags = 3, num_ctags = 1;
  const char *etag_names[] = {"MATERIAL_SET", "DIRICHLET_SET", "NEUMANN_SET"};
  const char *etag_vals[] = {NULL, NULL, NULL};
  const char *ctag_names[] = {"GEOM_DIMENSION"};
  const char *ctag_vals[] = {(const char*)&three};
  err = get_copy_expand_sets(cm, orig_set, etag_names, etag_vals, num_etags, CopyMesh::EXPAND);
  ERRORR("Failed to add expand lists.", iBase_FAILURE);
  err = get_copy_expand_sets(cm, orig_set, ctag_names, ctag_vals, num_ctags, CopyMesh::COPY);
  ERRORR("Failed to add expand lists.", iBase_FAILURE);
  err = extend_expand_sets(cm);
  ERRORR("Failed to extend expand lists.", iBase_FAILURE);
    
    // copy
  iBase_EntityHandle *new_ents;
  int new_ents_alloc, new_ents_size;
  for (int irow = 0; irow < NRINGS; irow++) {
    for (int icol = 0; icol < NRINGS; icol++) {

      if (!sixth_core[irow*NRINGS + icol]) continue;
      dx[0] = .5*irow*ASSY_PITCH + icol * ASSY_PITCH; dx[1] = irow * Y_SEP;
      new_ents = NULL;
      new_ents_alloc = 0;
      
      err = cm->copy_move_entities(orig_ents, orig_ents_size, dx, 
                                   &new_ents, &new_ents_alloc, &new_ents_size,
                                   false);
      ERRORR("Failed to copy_move entities.", 1);
      std::cout << "Copy/moved irow=" << irow << ", icol=" << icol << std::endl;
      free(new_ents);

      err = cm->tag_copied_sets(ctag_names, ctag_vals, 1);
      ERRORR("Failed to tag copied sets.", iBase_FAILURE);
    }
  }
    
    // merge
  std::cout << "Merged entities." << std::endl;
  err = mm->merge_entities(NULL, 0, 1.0e-8);
  if (0 != err) std::cerr << "Failed to merge entities." << std::endl;

    // assign new global ids
  std::cout << "Assigning global ids." << std::endl;
  MKUtils mu(impl);
  err = mu.assign_global_ids(root_set, 3, 1, true, false,
                             "GLOBAL_ID");
  ERRORR("Error assigning global ids.", err);

    // export
  std::cout << "Saving in file " << argv[2] << std::endl;
  iMesh_save(impl, root_set, argv[2], NULL, &err, 
             strlen(argv[2]), 0);
  ERRORR("Failed to save mesh.", 1);
  std::cout << "Wrote " << argv[2] << std::endl;

  delete cm;
  delete mm;
  free(orig_ents);
  
  iMesh_dtor(impl, &err);
  ERRORR("Destructor failed.", 1);
  
  return 0;
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


