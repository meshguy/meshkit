/*********************************************
 June,10
 Reactor Assembly Mesh Assembler
 Argonne National Laboratory

 CCrgen class functions for cop/moving meshes
 based on symmetry and geometry type
*********************************************/

#include "crgen.hpp"
#include <string.h>
int CCrgen::copymove(const int nrank, const int numprocs)
// ---------------------------------------------------------------------------
// Function: copy/move the assemblies based on the geometrytype and symmetry - Assume 1 meshfile in each instance
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
  if (nrank < ((int) core_alias.size() + ((int) files.size() - nassys))){
    err = set_copymove_coords();
    ERRORR("Failed to set cm coords.", err);
  
    // now copy/move
    err = copymove_all(nrank, numprocs);
    ERRORR("Failed to cm hexflat.", err);
  }
  return iBase_SUCCESS;
}

int CCrgen::set_copymove_coords()
// ---------------------------------------------------------------------------
// Function:
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
  x_coord.resize(tot_assys);
  y_coord.resize(tot_assys);  
  int i = 0;
  int assm_index;
  double dx[3] = {0.0, 0.0, 0.0};
  // move the assys based on the geometry type
  if (!strcmp(geom_type.c_str(), "hexflat") && symm == 6) {
    for (int n1 = 0; n1 < nrings; n1++) {
      for (int n2 = 0; n2 < n1 + 1; n2++) {
	err = find_assm(i, assm_index);
	y_coord[i] = n1 * pitch * sin(PII / 3.0) - n2 * pitch * sin(PII / 3.0);
	x_coord[i] = n1 * pitch * cos(PII / 3.0) + n2 * pitch * cos(PII / 3.0);
	i++;
      }
    }
  }
  if (!strcmp(geom_type.c_str(), "rectangular") && symm == 1) {
    for (int n1 = 0; n1 < nringsx; n1++) {
      for (int n2 = 0; n2 < nringsy; n2++) {
	err = find_assm(i, assm_index);
	y_coord[i] = -n1 * pitchy;
	x_coord[i] = n2 * pitchx;
	i++;
      }
    }
  }
  if (!strcmp(geom_type.c_str(), "hexflat") && symm == 1) {
    int t, width = 2 * nrings - 1;
    for (int n1 = 1; n1 <= width; n1++) {
      if(n1 > nrings)
	t = 2 * nrings - n1;
      else
	t = n1;

      for (int n2 = 1; n2 <= (nrings + t - 1); n2++) {
	err = find_assm(i,assm_index);

	if (n1 < nrings){
	  x_coord[i] = (nrings - n2 + 1) * pitch / 2.0 + n2 * pitch / 2.0 + 
	    (n2 - 1) * pitch - (n1 - 1) * pitch / 2.0;
	  y_coord[i] = -((n1 - 1) * (0.5 * pitch / sin(PII/3.0) + 0.5 * pitch * sin(PII/6.0) / sin(PII/3.0)));
	}
	else{
	  x_coord[i] = (nrings - n2 + 1) * pitch / 2.0 + n2 * pitch / 2.0 + (n2 - 1) * pitch -
	    (2 * nrings - n1 -1) * pitch / 2.0;
	  y_coord[i] = -((n1 -1) * (0.5 * pitch / sin(PII/3.0) + 0.5 * pitch * sin(PII/6.0) / sin(PII/3.0))); 
	}
	i++;
      }
    }
  }

  if (!strcmp(geom_type.c_str(), "hexflat") && symm == 12) {
    int flag = 0;
    for (int n1 = 0; n1 < nrings; n1++) {
      int loc = (n1 + 2)/2;    

      if( flag == 0 ){
	dx[0] = (n1 + loc - 1) * pitch / 2.0;
	dx[1] = (n1 - loc + 1) * pitch * sin(PII/3.0);
	flag = 1;
      }
      else{
	dx[0] = (n1 + loc) * pitch / 2.0;
	dx[1] = (n1 - loc) * pitch * sin(PII/3.0);
	flag = 0;
      }

      for (int n2 = 0; n2 < loc; n2++) {
	err = find_assm(i,assm_index);
	y_coord[i] = dx[1] - n2 * pitch * sin(PII/3.0);
	x_coord[i] = dx[0] + n2 * pitch * cos(PII/3.0);
	i++;
      }
    }
  }

  if (!strcmp(geom_type.c_str(), "hexvertex") && symm == 6) {
    int bd = 0;
    for (int n1 = 0; n1 < nrings; n1++) {
      if(n1%2==0){//check if n1 is even
	for (int n2 = 0; n2 < n1+1; n2++) {

	  err = find_assm(i, assm_index);
	  if (-1 == assm_index){
	    i++;
	    if(n2 > (n1+1)/2)
	      ++bd; // index for assemblies below diagonal needs updatation
	    continue;
	  }
	  if(n2 <= n1/2){// before or equal to diagonal
	    dx[0] = n2 * pitch;
	    dx[1] = n1 * pitch * sin(PII/3.0);
	  }
	  else{//below the diagonal
	    dx[0] = (n1 + 1 + bd) * pitch / 2.0; 
	    dx[1] = (n1 - 1 - bd) * pitch * sin(PII/3.0);
	    ++bd;      
	  }      
	  x_coord[i] = (dx[0] * cos(PII/6.0) + dx[1] * sin(PII/6.0));
	  y_coord[i] = (dx[1] * cos(PII/6.0) - dx[0] * sin(PII/6.0));  
	  i++;
	}
      }
      else{//n1 is odd
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
	    dx[1] = n1 * pitch * sin(PII/3.0);

	  }
	  else{//below the diagonal 
	    dx[0] = (n1 + 1 + bd) * pitch / 2.0; 
	    if (bd == 0) // first n2 = 1 assembly
	      dx[1] = pitch * sin(PII/3.0);
	    dx[1] = (n1 - 1 - bd) * pitch * sin(PII/3.0);
	    ++bd;    
	  }


	  // starting from x-axis
	  x_coord[i] = (dx[0] * cos(PII/6.0) + dx[1] * sin(PII/6.0));
	  y_coord[i] = (dx[1] * cos(PII/6.0) - dx[0] * sin(PII/6.0));
	  i++;
	}
      }
      bd = 0;
    }
  }
  return iBase_SUCCESS;
}

int CCrgen::copymove_all(const int nrank, const int numprocs)
// ---------------------------------------------------------------------------
// Function:
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
  if(prob_type == "mesh"){
    // get the copy/expand sets
    int num_etags = 3, num_ctags = 1;
    const char     *etag_names[] = { "MATERIAL_SET", "DIRICHLET_SET", "NEUMANN_SET" };
    const char *etag_vals[] = { NULL, NULL, NULL };
    const char *ctag_names[] = { "GEOM_DIMENSION" };
    const char *ctag_vals[] = { (const char*) &set_DIM };
  
    int flag = 1;
    int assm_index = -1;
    double dx_orig[3], dx[3];
    iBase_EntityHandle *new_ents;
    int new_ents_alloc, new_ents_size;

    if(numprocs <= (int) files.size()){
      // no distribution of task for copy/move; each file loaded only once 
      int flags[assys.size()], move_index;
      double dx[3] = { 0.0, 0.0, 0.0 };
      double dx_move[3] = { 0.0, 0.0, 0.0 };
      CMatrix<double> dx_orig(assys.size(), 3);
      dx_orig.Set(0.0);

      for (unsigned int i = 0; i < assys.size(); i++) {
	flags[i]=0;
	err = get_expand_sets(cm[i], assys[i], etag_names, etag_vals, num_etags);
	ERRORR("Failed to add expand lists.", iBase_FAILURE);
	err = get_copy_sets(cm[i], assys[i], ctag_names, ctag_vals, num_ctags);
	ERRORR("Failed to add expand lists.", iBase_FAILURE);
	err = extend_expand_sets(cm[i]);
	ERRORR("Failed to extend expand lists.", iBase_FAILURE);
      }
      for (int k=0; k< tot_assys; k++){
	err = find_assm(k, assm_index);
	if(assm_index >= 0){
	  // check if this file is with the proc
	  int out = 0;
	  for(int c=0; c < (int) assys.size(); c++){
	    if(assys_index[c] == assm_index){
	      //found assembly
	      move_index=c;
	      out = 0;
	      break;
	    }
	    else{
	      out = 1;
	    }
	  }
	  if(out == 1){
	    continue;
	  }
	  if (-1 == assm_index) {
	    continue;
	  }
	
	  if (flags[move_index] == 0) {
	    dx_move[0] = x_coord[k];
	    dx_move[1] = y_coord[k];
	    dx_move[2] = 0;
	    dx_orig(move_index + 1, 1) = dx_move[0];
	    dx_orig(move_index + 1, 2) = dx_move[1];
	    dx_orig(move_index + 1, 3) = dx_move[2];

	    move_verts(assys[move_index], dx_move);
	    std::cout << "Moved Assembly: " << assm_index  << " dX = " << dx_move[0] << " dY = "
		      << dx_move[1] << " in proc with rank " << nrank << std::endl;
	    if(strcmp(core_info.c_str(),"on") == 0)
	      info_file << assm_index << " \t" << k  << " \t" << dx_move[0] << " \t" << dx_move[1]  << " \t" << dx_move[2]  << " \t" << nrank << std::endl;

	    flags[move_index]=1;
	  }
	  else{

	    dx[0] =  x_coord[k] - dx_orig(move_index+1, 1);
	    dx[1] =  y_coord[k] - dx_orig(move_index+1, 2);
	    dx[2] =  0.0;

	    int orig_ents_alloc = 0, orig_ents_size;
	    iBase_EntityHandle *orig_ents = NULL;
	    new_ents = NULL;
	    new_ents_alloc = 0;

	    iMesh_getEntities(impl, assys[move_index], iBase_ALL_TYPES,iMesh_ALL_TOPOLOGIES,
			      &orig_ents, &orig_ents_alloc, &orig_ents_size, &err);
	    ERRORR("Failed to get any entities from original set.", iBase_FAILURE);

	    cm[move_index]->copy(orig_ents,orig_ents_size, copy::Translate(dx), 
				 &new_ents, &new_ents_alloc, &new_ents_size, false);

	    std::cout << "Copy/moved A: " << assm_index 
		      <<" dX = " <<dx[0]<< " dY = " << dx[1] << " rank " << nrank << std::endl;
	    if(strcmp(core_info.c_str(),"on") == 0)
	      info_file << assm_index << " \t" << k  << " \t" << dx[0] << " \t" << dx[1]  << " \t" << dx[2]  << " \t" << nrank << std::endl;
	    free(new_ents);
	    free(orig_ents);

	    cm[move_index]->tag_copied_sets(ctag_names, ctag_vals, 1);
	  }

	}
      }
    }
    else{
      for(int i =0; i < (int) position_core[nrank].size(); i++){
	assm_index = position_core[nrank][i];
	if(assm_index >= 0){
	  // some files are loaded in two or more processors, copy/move task distribution takes place
	  int j = 0;
	  err = get_expand_sets(cm[j], assys[j], etag_names, etag_vals, num_etags);
	  ERRORR("Failed to add expand lists.", iBase_FAILURE);
	  err = get_copy_sets(cm[j], assys[j], ctag_names, ctag_vals, num_ctags);
	  ERRORR("Failed to add expand lists.", iBase_FAILURE);
	  err = extend_expand_sets(cm[j]);
	  ERRORR("Failed to extend expand lists.", iBase_FAILURE);
	
	  if(flag == 0){
	    dx[0] = x_coord[assm_index] - dx_orig[0];
	    dx[1] = y_coord[assm_index] - dx_orig[1];
	    dx[2] = 0.0;
      
	    new_ents = NULL;
	    new_ents_alloc = 0;
	    new_ents_size = 0;
      
	    int orig_ents_alloc = 0, orig_ents_size;
	    iBase_EntityHandle *orig_ents = NULL;
      
      
	    iMesh_getEntities(impl, assys[0], iBase_ALL_TYPES,
			      iMesh_ALL_TOPOLOGIES, &orig_ents, &orig_ents_alloc,
			      &orig_ents_size, &err);
	    ERRORR("Failed to get any entities from original set.", iBase_FAILURE);
      
	    cm[0]->copy(orig_ents, orig_ents_size, copy::Translate(dx),
			&new_ents, &new_ents_alloc, &new_ents_size, false);
	    std::cout << "Copy/moved Assm: " << assm_index << " dX = " << dx[0] << " dY = "
		      << dx[1]  << " rank " << nrank << std::endl;
	    if(strcmp(core_info.c_str(),"on") == 0)
	      info_file << assm_index << " \t" << i  << " \t" << dx[0] << " \t" << dx[1]  << " \t" << dx[2]  << " \t" << nrank << std::endl;

	    cm[0]->tag_copied_sets(ctag_names, ctag_vals, 1);
      
	    free(new_ents);
	    free(orig_ents);
	  } else {
	    flag = 0;
	    dx_orig[0] = x_coord[assm_index];
	    dx_orig[1] = y_coord[assm_index];
	    dx_orig[2] = 0;
	    move_verts(assys[0], dx_orig);
	    std::cout << "Moved Assm: " << assm_index << " dX = " << dx_orig[0] << " dY = "
		      << dx_orig[1] << " rank " << nrank << std::endl;
	    if(strcmp(core_info.c_str(),"on") == 0)
	      info_file << assm_index << " \t" << i  << " \t" << dx_orig[0] << " \t" << dx_orig[1]  << " \t" << dx_orig[2]  << " \t" << nrank << std::endl;

	  }
	}
      }
    }
  }
  else{ // prob type is geometry
    int assm_index = -1;
    iBase_EntityHandle *new_ents;
    int new_ents_alloc, new_ents_size;

    // no distribution of task for copy/move; each file loaded only once 
    int flags[assys.size()], move_index;
    double dx[3] = { 0.0, 0.0, 0.0 };
    double dx_move[3] = { 0.0, 0.0, 0.0 };
    CMatrix<double> dx_orig(assys.size(), 3);
    dx_orig.Set(0.0);

    for (unsigned int i = 0; i < assys.size(); i++) {
      flags[i]=0;
    }

    for (int k=0; k< tot_assys; k++){
      err = find_assm(k, assm_index);
      if(assm_index >= 0){
	// check if this file is with the proc
	int out = 0;
	for(int c=0; c < (int) assys.size(); c++){
	  if(assys_index[c] == assm_index){
	    //found assembly
	    move_index=c;
	    out = 0;
	    break;
	  }
	  else{
	    out = 1;
	  }
	}
	if(out == 1){
	  continue;
	}
	if (-1 == assm_index) {
	  continue;
	}
	
	if (flags[move_index] == 0) {
	  dx_move[0] = x_coord[k];
	  dx_move[1] = y_coord[k];
	  dx_move[2] = 0;
	  dx_orig(move_index + 1, 1) = dx_move[0];
	  dx_orig(move_index + 1, 2) = dx_move[1];
	  dx_orig(move_index + 1, 3) = dx_move[2];

	  move_geoms(assys[move_index], dx_move);
	  std::cout << "Moved Assembly: " << assm_index  << " dX = " << dx_move[0] << " dY = "
		    << dx_move[1] << " in proc with rank " << nrank << std::endl;
	  flags[move_index]=1;
	}
	else{

	  dx[0] =  x_coord[k] - dx_orig(move_index+1, 1);
	  dx[1] =  y_coord[k] - dx_orig(move_index+1, 2);
	  dx[2] =  0.0;

	  int orig_ents_alloc = 0, orig_ents_size;
	  iBase_EntityHandle *orig_ents = NULL;
	  new_ents = NULL;
	  new_ents_alloc = 0;

	  iGeom_getEntities(geom, assys[move_index], iBase_ALL_TYPES,
			    &orig_ents, &orig_ents_alloc, &orig_ents_size, &err);
	  ERRORR("Failed to get any entities from original set.", iBase_FAILURE);

	  cg[move_index]->copy(orig_ents,orig_ents_size, dx, 
			       &new_ents, &new_ents_alloc, &new_ents_size, false);

	  std::cout << "Copy/moved A: " << assm_index 
		    <<" dX = " <<dx[0]<< " dY = " << dx[1] << " rank " << nrank << std::endl;
	  free(new_ents);
	  free(orig_ents);

	}

      }
    }
      
  }
  return iBase_SUCCESS;
}
