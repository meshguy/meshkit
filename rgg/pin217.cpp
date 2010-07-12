#include <cfloat>
#include <cmath>
#include "MKUtils.hpp"
#include "MergeMesh.hpp"
#include "CopyMesh.hpp"
#include "ExtrudeMesh.hpp"
#include "utils.hpp"

#define DEFAULT_TEST_FILE "pin1.cub"
#define DEFAULT_OUTPUT_FILE "pin217-out.h5m"

iMesh_Instance impl;
iBase_EntitySetHandle root_set;

int main(int argc, char **argv) 
{
  //set the dimension based on input mesh - 2D meshes or 3D mesh??
  int set_DIM = 3;
  std::cout <<"\nExecuting Pin217 Program\n\n";
  int err;
  // make a mesh instance
  iMesh_newMesh("MOAB", &impl, &err, 4);
  ERRORR("Failed to create instance.", 1);
  
  iMesh_getRootSet(impl, &root_set, &err);

  CopyMesh *cm = new CopyMesh(impl);
  MergeMesh *mm = new MergeMesh(impl);
  ExtrudeMesh *ext = new ExtrudeMesh(impl);
  char* input_file_name = NULL;
  char* output_file_name = NULL;

  if (3 == argc) {
    input_file_name = argv[1];
    output_file_name = argv[2];
  }
  else {
    if (argc < 1) return 1;
    std::cerr << "Usage: " << argv[0] << " <meshfile>" << std::endl;
    std::cout <<"  No file specified.  Defaulting to: " << DEFAULT_TEST_FILE << std::endl;
    input_file_name = DEFAULT_TEST_FILE;
    output_file_name = DEFAULT_OUTPUT_FILE;
  }

  std::vector<std::string> ctags, etags, utags;
  std::vector<char*> cvals, evals;
  std::string infile, outfile;
  //  double dx[3] = {0.0, 0.0, 0.0};
  double dx[2] = {0.0, 0.0};
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
  /*    bool sixth_core[] = {
      false, true, true, true, true, true, true, true, 
      true, true, true, true, true, true, true, false,
      true, true, true, true, true, true, true, false,
      true, true, true, true, true, true, false, false,
      true, true, true, true, true, false, false, false,
      true, true, true, true, false, false, false, false,
      true, true, true, false, false, false, false, false,
      true, false, false, false, false, false, false, false 
      };*/

  iBase_EntitySetHandle orig_set;
  iMesh_createEntSet(impl, 0, &orig_set, &err);
  ERRORR("Couldn't create orig entity set.", err);
  
    // read the file 
  iMesh_load(impl, orig_set, input_file_name, NULL, &err, strlen(input_file_name), 0);
  ERRORR("Couldn't read mesh file.", err);

  std::cout << "Loaded initial mesh." << std::endl;

    // get the entities we want to copy, basically all the non-set entities
  int orig_ents_alloc = 0, orig_ents_size; 
  iBase_EntityHandle *orig_ents = NULL;

  iMesh_getEntities(impl, orig_set, iBase_ALL_TYPES,iMesh_ALL_TOPOLOGIES,
                 &orig_ents, &orig_ents_alloc, &orig_ents_size, &err);
   ERRORR("Failed to get any entities from original set.", iBase_FAILURE);
  
    // get the copy/expand sets
  int num_etags = 3, num_ctags = 1;
  const char *etag_names[] = {"MATERIAL_SET", "DIRICHLET_SET", "NEUMANN_SET"};
  const char *etag_vals[] = {NULL, NULL, NULL};
  const char *ctag_names[] = {"GEOM_DIMENSION"};
  const char *ctag_vals[]={(const char*)&set_DIM};
  err = get_expand_sets(cm, orig_set, etag_names, etag_vals, num_etags);
  ERRORR("Failed to add expand lists.", iBase_FAILURE);
  err = get_copy_sets(cm, orig_set, ctag_names, ctag_vals, num_ctags);
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
      
      cm->copy_move_entities(orig_ents, orig_ents_size, dx, 
                             &new_ents, &new_ents_alloc, &new_ents_size,
                             false);
      std::cout << "Copy/moved irow=" << irow << ", icol=" << icol << std::endl;
      free(new_ents);
      
      cm->tag_copied_sets(ctag_names, ctag_vals, 1);
    }
  }
       //getting hexahedron elements for merge_entities   
    const double merge_tol =  1.0e-8;
    const int do_merge = 1; //0 or 1 for off or on respectively 
    const int update_sets= 0; 
    iBase_TagHandle merge_tag = NULL;
    iBase_EntityHandle *ents = NULL;
    int ents_alloc = 0, ents_size;
 if(set_DIM ==2){
    iMesh_getEntities(impl, root_set,
                    iBase_FACE, iMesh_ALL_TOPOLOGIES,
                    &ents, &ents_alloc, &ents_size, &err);
     ERRORR("Failed to get entities from set recursively.", err);
 }   
  else{
  iMesh_getEntities(impl, root_set,
                    iBase_REGION, iMesh_ALL_TOPOLOGIES,
                    &ents, &ents_alloc, &ents_size, &err);
     ERRORR("Failed to get entities from set recursively.", err);
  }
 
      // merge  
    int num1, num2;
  
    iMesh_getNumOfType(impl, root_set, iBase_VERTEX, &num1, &err);
    ERRORR("Trouble getting number of entities after merge.", err);

    mm->merge_entities(ents, ents_size, merge_tol,
                       do_merge, update_sets, merge_tag);
    
    iMesh_getNumOfType(impl, root_set, iBase_VERTEX, &num2, &err);
    ERRORR("Trouble getting number of entities after merge.", err);

    std::cout << "Merged " << num1 - num2 << " vertices." << std::endl;

  // extrude the 2D mesh now
  iBase_EntityHandle faces[] = {*ents };
  double v[] = { 0, 0, 1 };
  int steps = 5;
  // err = ext->translate(faces, ents_size, steps, v);
  ERRORR("Couldn't extrude mesh", 1);
 
   // assign new global ids
  std::cout << "Assigning global ids." << std::endl;
  MKUtils mu(impl);
  err = mu.assign_global_ids(root_set, 3, 1, true, false,
                             "GLOBAL_ID");
  ERRORR("Error assigning global ids.", err);

    // export
  std::cout << "Saving in file " << output_file_name << std::endl;
  iMesh_save(impl, root_set, output_file_name, NULL, &err, 
             strlen(output_file_name), 0);
  ERRORR("Failed to save mesh.", 1);
  std::cout << "Wrote " << output_file_name << std::endl;

  delete cm;
  delete mm;
  free(orig_ents);
   
  iMesh_dtor(impl, &err);
  ERRORR("Destructor failed.", 1);
  
  return 0;
}
