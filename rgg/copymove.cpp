/*********************************************
June,10
Reactor Assembly Mesh Assembler
Argonne National Laboratory

CCrgen class functions for cop/moving meshes
based on symmetry and geometry type
*********************************************/

#include "crgen.hpp"
#include <string.h>
int CCrgen::copy_move()
// ---------------------------------------------------------------------------
// Function: copy/move the assemblies based on the geometrytype and symmetry
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
  // move the assys based on the geometry type
  if(!strcmp(geom_type.c_str(), "hexvertex") && symm == 6){
    if(strcmp(prob_type.c_str(), "mesh")==0){
      err = copy_move_hex_vertex_assys_p1(cm, nrings, pack_type, pitch, symm, 
					  core_alias, assys);
      ERRORR("Failed in copy/move step.", err);

      err = copy_move_hex_vertex_assys(cm, nrings, pack_type, pitch, symm, 
				       core_alias, assys);
      ERRORR("Failed in copy/move step.", err);
    }
    else if(strcmp(prob_type.c_str(), "geometry")==0){
      err = copy_move_hex_vertex_assys_p1(cg, nrings, pack_type, pitch, symm,
					  core_alias, assys);
      ERRORR("Failed in copy/move step.", err);

      err = copy_move_hex_vertex_assys(cg, nrings, pack_type, pitch, symm,
				       core_alias, assys);
      ERRORR("Failed in copy/move step.", err);
    }
  }
  else if(!strcmp(geom_type.c_str(),"rectangular") && symm == 1){

    if(strcmp(prob_type.c_str(), "mesh")==0){
      err = copy_move_sq_assys_p1(cm, nrings, pack_type, pitch, symm,
				  core_alias, assys);
      ERRORR("Failed in copy/move step.", err);

      err = copy_move_sq_assys(cm, nrings, pack_type, pitch, symm, 
			       core_alias, assys);
      ERRORR("Failed in copy/move step.", err);
    }
    else if(strcmp(prob_type.c_str(), "geometry")==0){
      err = copy_move_sq_assys_p1(cg, nrings, pack_type, pitch, symm,
				  core_alias, assys);
      ERRORR("Failed in copy/move step.", err);

      err = copy_move_sq_assys(cg, nrings, pack_type, pitch, symm,
			       core_alias, assys);
      ERRORR("Failed in copy/move step.", err);
    }
  } 
  else if(!strcmp(geom_type.c_str(),"hexflat") && symm == 6){
    if(strcmp(prob_type.c_str(), "mesh")==0){
      err = copy_move_hex_flat_assys_p1(cm, nrings, pack_type, pitch, symm, 
					core_alias, assys);
      ERRORR("Failed in copy/move step.", err);

      err = copy_move_hex_flat_assys(cm, nrings, pack_type, pitch, symm, 
				     core_alias, assys);
      ERRORR("Failed in copy/move step.", err);
    }
    else if(strcmp(prob_type.c_str(), "geometry")==0){
      err = copy_move_hex_flat_assys_p1(cg, nrings, pack_type, pitch, symm,
					core_alias, assys);
      ERRORR("Failed in copy/move step.", err);

      err = copy_move_hex_flat_assys(cg, nrings, pack_type, pitch, symm,
				     core_alias, assys);
      ERRORR("Failed in copy/move step.", err);
    }
  }
  else if(!strcmp(geom_type.c_str(),"hexflat") && symm == 1){
    if(strcmp(prob_type.c_str(), "mesh")==0){
      err = copy_move_hex_full_assys_p1(cm, nrings, pack_type, pitch, symm, 
					core_alias, assys);
      ERRORR("Failed in copy/move step.", err);

      err = copy_move_hex_full_assys(cm, nrings, pack_type, pitch, symm, 
				     core_alias, assys);
      ERRORR("Failed in copy/move step.", err);
    }
    if(strcmp(prob_type.c_str(), "geometry")==0){
      err = copy_move_hex_full_assys_p1(cg, nrings, pack_type, pitch, symm,
					core_alias, assys);
      ERRORR("Failed in copy/move step.", err);

      err = copy_move_hex_full_assys(cg, nrings, pack_type, pitch, symm,
				     core_alias, assys);
      ERRORR("Failed in copy/move step.", err);
    }
  }
  else if(!strcmp(geom_type.c_str(),"hexflat") && symm == 12){
    if(strcmp(prob_type.c_str(), "mesh")==0){
      err = copy_move_one_twelfth_assys_p1(cm, nrings, pack_type, pitch, symm, 
					   core_alias, assys);
      ERRORR("Failed in copy/move step.", err);

      err = copy_move_one_twelfth_assys(cm, nrings, pack_type, pitch, symm, 
					core_alias, assys);
      ERRORR("Failed in copy/move step.", err);
    }
    if(strcmp(prob_type.c_str(), "geometry")==0){
      err = copy_move_one_twelfth_assys_p1(cg, nrings, pack_type, pitch, symm,
					   core_alias, assys);
      ERRORR("Failed in copy/move step.", err);

      err = copy_move_one_twelfth_assys(cg, nrings, pack_type, pitch, symm,
					core_alias, assys);
      ERRORR("Failed in copy/move step.", err);
    }
  }
  return iBase_SUCCESS;
}



int CCrgen::copy_move_hex_vertex_assys_p1(CopyMesh **cm,
					  const int nrings, const int pack_type,
					  const double pitch,
					  const int symm,
					  std::vector<std::string> &core_alias,
					  std::vector<iBase_EntitySetHandle> &assys)
{
  double dx[3] = {0.0, 0.0, 0.0};
  double dxnew[3] = {0.0, 0.0, 0.0};
  int i = 0, bd = 0;
  int assm_index;
  int flags[files.size()];

  for (unsigned int i = 0; i < files.size(); i++) {
    flags[i] = 0;
  }

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

	// starting from x-axis
	dxnew[0] = (dx[0] * cos(PI/6.0) + dx[1] * sin(PI/6.0));
	dxnew[1] = (dx[1] * cos(PI/6.0) - dx[0] * sin(PI/6.0));      

	if(flags[assm_index]==0){
	  move_verts(assys[assm_index], dxnew);

	  std::cout << "Copy/moved A: " << assm_index 
		    <<  " n1=" << n1 << ", n2=" << n2  <<" dX = " <<dx[0]<< " dY = " << dx[1] << std::endl;
	}
	i++;
	flags[assm_index] = 1;
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


	// starting from x-axis
	dxnew[0] = (dx[0] * cos(PI/6.0) + dx[1] * sin(PI/6.0));
	dxnew[1] = (dx[1] * cos(PI/6.0) - dx[0] * sin(PI/6.0));


	if(flags[assm_index]==0){
	  move_verts(assys[assm_index], dxnew);

	  std::cout << "Copy/moved A: " << assm_index 
		    <<  " n1=" << n1 << ", n2=" << n2  <<" dX = " <<dx[0]<< " dY = " << dx[1] << std::endl;
	}
	i++;
	flags[assm_index] = 1;
      }
    }
    bd = 0;
  }
  return iBase_SUCCESS;
}  


int CCrgen::copy_move_hex_vertex_assys(CopyMesh **cm,
				       const int nrings, const int pack_type,
				       const double pitch,
				       const int symm,
				       std::vector<std::string> &core_alias,
				       std::vector<iBase_EntitySetHandle> &assys)
{
  double dx[3] = {0.0, 0.0, 0.0};
  double dxnew[3] = {0.0, 0.0, 0.0};
  iBase_EntityHandle *new_ents;
  int new_ents_alloc, new_ents_size;
  int i = 0, bd = 0;
  int assm_index;

  int flags[files.size()];
  CMatrix<double> dx_orig(files.size(), 3);
  dx_orig.Set(0.0);

  // get the copy/expand sets
  int num_etags = 3, num_ctags = 1;
  const char *etag_names[] = {"MATERIAL_SET", "DIRICHLET_SET", "NEUMANN_SET"};
  const char *etag_vals[] = {NULL, NULL, NULL};
  const char *ctag_names[] = {"GEOM_DIMENSION"};
  const char *ctag_vals[]={(const char*)&set_DIM};

  for (unsigned int i = 0; i < files.size(); i++) {
    flags[i]=0;
    err = get_expand_sets(cm[i],assys[i], etag_names, etag_vals, num_etags);
    ERRORR("Failed to add expand lists.", iBase_FAILURE);
    err = get_copy_sets(cm[i],assys[i], ctag_names, ctag_vals, num_ctags);
    ERRORR("Failed to add expand lists.", iBase_FAILURE);
    err = extend_expand_sets(cm[i]);
    ERRORR("Failed to extend expand lists.", iBase_FAILURE);
  }  

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

	// starting from x-axis
	dxnew[0] = (dx[0] * cos(PI/6.0) + dx[1] * sin(PI/6.0));
	dxnew[1] = (dx[1] * cos(PI/6.0) - dx[0] * sin(PI/6.0));      

	if(flags[assm_index]==1){
	  dxnew[0]-=dx_orig(assm_index+1, 1);
	  dxnew[1]-=dx_orig(assm_index+1, 2);
	  dxnew[2]-=dx_orig(assm_index+1, 3);

	  int orig_ents_alloc = 0, orig_ents_size;
	  iBase_EntityHandle *orig_ents = NULL;

	  iMesh_getEntities(impl, assys[assm_index], iBase_ALL_TYPES,iMesh_ALL_TOPOLOGIES,
			    &orig_ents, &orig_ents_alloc, &orig_ents_size, &err);
	  ERRORR("Failed to get any entities from original set.", iBase_FAILURE);


	  cm[assm_index]->copy(orig_ents,orig_ents_size, copy::Translate(dxnew), 
	                       &new_ents, &new_ents_alloc, &new_ents_size, false);
	  std::cout << "Copy/moved A: " << assm_index 
		    << " n1=" << n1 << ", n2=" << n2 <<" dX = " <<dxnew[0]<< " dY = " << dxnew[1] << std::endl;
	  free(new_ents);
	  free(orig_ents);
	  i++;

	  cm[assm_index]->tag_copied_sets(ctag_names, ctag_vals, 1);
	}
	else{
	  i++;
	  flags[assm_index]=1;
	  dx_orig(assm_index+1, 1)=dxnew[0];
	  dx_orig(assm_index+1, 2)=dxnew[1];
	  dx_orig(assm_index+1, 3)=dxnew[2];
	}
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


	// starting from x-axis
	dxnew[0] = (dx[0] * cos(PI/6.0) + dx[1] * sin(PI/6.0));
	dxnew[1] = (dx[1] * cos(PI/6.0) - dx[0] * sin(PI/6.0));


	if(flags[assm_index]==1){
	  dxnew[0]-=dx_orig(assm_index+1, 1);
	  dxnew[1]-=dx_orig(assm_index+1, 2);
	  dxnew[2]-=dx_orig(assm_index+1, 3);

	  int orig_ents_alloc = 0, orig_ents_size;
	  iBase_EntityHandle *orig_ents = NULL;

	  iMesh_getEntities(impl, assys[assm_index], iBase_ALL_TYPES,iMesh_ALL_TOPOLOGIES,
			    &orig_ents, &orig_ents_alloc, &orig_ents_size, &err);
	  ERRORR("Failed to get any entities from original set.", iBase_FAILURE);

	  cm[assm_index]->copy(orig_ents,orig_ents_size, copy::Translate(dxnew), 
	                       &new_ents, &new_ents_alloc, &new_ents_size, false);
	  std::cout << "Copy/moved A: " << assm_index 
		    << " n1=" << n1 << ", n2=" << n2 <<" dX = " <<dxnew[0]<< " dY = " << dxnew[1] << std::endl;
	  free(new_ents);
	  free(orig_ents);
	  i++;

	  cm[assm_index]->tag_copied_sets(ctag_names, ctag_vals, 1);
	}
	else{
	  i++;
	  flags[assm_index]=1;
	  dx_orig(assm_index+1, 1)=dxnew[0];
	  dx_orig(assm_index+1, 2)=dxnew[1];
	  dx_orig(assm_index+1, 3)=dxnew[2];
	}
    
      }
    }
    bd = 0;
  }
  return iBase_SUCCESS;
}  
int CCrgen::copy_move_one_twelfth_assys_p1(CopyMesh **cm,
					   const int nrings, const int pack_type,
					   const double pitch,
					   const int symm,
					   std::vector<std::string> &core_alias,
					   std::vector<iBase_EntitySetHandle> &assys)
{
  double dx[3] = {0.0, 0.0, 0.0};
  double dxnew[3] = {0.0, 0.0, 0.0};
  int i = 0, flag = 0;
  int flags[files.size()];
  int assm_index;

  for (unsigned int i = 0; i < files.size(); i++) {
    flags[i] = 0;
  }  


  for (int n1 = 0; n1 < nrings; n1++) {
    int loc = (n1 + 2)/2;    

    if( flag == 0 ){
      dx[0] = (n1 + loc - 1) * pitch / 2.0;
      dx[1] = (n1 - loc + 1) * pitch * sin(PI/3.0);
      flag = 1;
    }
    else{
      dx[0] = (n1 + loc) * pitch / 2.0;
      dx[1] = (n1 - loc) * pitch * sin(PI/3.0);
      flag = 0;
    }

    for (int n2 = 0; n2 < loc; n2++) {
      err = find_assm(i,assm_index);
      if (-1 == assm_index) {
        i++;
        continue;
      }  
    
      dxnew[1] = dx[1] - n2 * pitch * sin(PI/3.0);
      dxnew[0] = dx[0] + n2 * pitch * cos(PI/3.0);
      if(flags[assm_index] == 0){

 	move_verts(assys[assm_index], dxnew);

	std::cout << "Copy/moved A: " << assm_index 
		  <<  " n1=" << n1 << ", n2=" << n2  <<" dX = " <<dxnew[0]<< " dY = " << dxnew[1] << std::endl;
      }
      i++;
      flags[assm_index] = 1;
    }
    dx[0] = 0.0;
    dx[1] = 0.0;
  }

  return iBase_SUCCESS;
}  

int CCrgen::copy_move_one_twelfth_assys(CopyMesh **cm,
					const int nrings, const int pack_type,
					const double pitch,
					const int symm,
					std::vector<std::string> &core_alias,
					std::vector<iBase_EntitySetHandle> &assys)
{
  double dx[3] = {0.0, 0.0, 0.0};
  double dxnew[3] = {0.0, 0.0, 0.0};
  iBase_EntityHandle *new_ents;
  int new_ents_alloc, new_ents_size;
  int i = 0, flag = 0;
  int assm_index;
  int flags[files.size()];
  CMatrix<double> dx_orig(files.size(), 3);
  dx_orig.Set(0.0);
  
  // get the copy/expand sets
  int num_etags = 3, num_ctags = 1;
  const char *etag_names[] = {"MATERIAL_SET", "DIRICHLET_SET", "NEUMANN_SET"};
  const char *etag_vals[] = {NULL, NULL, NULL};
  const char *ctag_names[] = {"GEOM_DIMENSION"};
  const char *ctag_vals[]={(const char*)&set_DIM};
 
  for (unsigned int i = 0; i < files.size(); i++) {
    flags[i]=0;
    err = get_expand_sets(cm[i],assys[i], etag_names, etag_vals, num_etags);
    ERRORR("Failed to add expand lists.", iBase_FAILURE);
    err = get_copy_sets(cm[i],assys[i], ctag_names, ctag_vals, num_ctags);
    ERRORR("Failed to add expand lists.", iBase_FAILURE);
    err = extend_expand_sets(cm[i]);
    ERRORR("Failed to extend expand lists.", iBase_FAILURE);
  }  


  for (int n1 = 0; n1 < nrings; n1++) {
    int loc = (n1 + 2)/2;    

    if( flag == 0 ){
      dx[0] = (n1 + loc - 1) * pitch / 2.0;
      dx[1] = (n1 - loc + 1) * pitch * sin(PI/3.0);
      flag = 1;
    }
    else{
      dx[0] = (n1 + loc) * pitch / 2.0;
      dx[1] = (n1 - loc) * pitch * sin(PI/3.0);
      flag = 0;
    }

    for (int n2 = 0; n2 < loc; n2++) {
      err = find_assm(i,assm_index);
      if (-1 == assm_index) {
        i++;
        continue;
      }  
    
      dxnew[1] = dx[1] - n2 * pitch * sin(PI/3.0);
      dxnew[0] = dx[0] + n2 * pitch * cos(PI/3.0);

      if(flags[assm_index]==1){
	dxnew[0]-=dx_orig(assm_index+1, 1);
	dxnew[1]-=dx_orig(assm_index+1, 2);
	dxnew[2]-=dx_orig(assm_index+1, 3);
	new_ents = NULL;
	new_ents_alloc = 0;
	new_ents_size = 0;

	int orig_ents_alloc = 0, orig_ents_size;
	iBase_EntityHandle *orig_ents = NULL;

	iMesh_getEntities(impl, assys[assm_index], iBase_ALL_TYPES,iMesh_ALL_TOPOLOGIES,
			  &orig_ents, &orig_ents_alloc, &orig_ents_size, &err);
	ERRORR("Failed to get any entities from original set.", iBase_FAILURE);

	cm[assm_index]->copy(orig_ents,orig_ents_size, copy::Translate(dxnew),
	                     &new_ents, &new_ents_alloc, &new_ents_size, false);
	std::cout << "Copy/moved A: " << assm_index 
		  << " n1=" << n1 << ", n2=" << n2 <<" dX = " <<dxnew[0]<< " dY = " << dxnew[1] << std::endl;

	free(new_ents);
	free(orig_ents);
	i++;
	dxnew[0] = 0.0;
	dxnew[1] = 0.0;
	cm[assm_index]->tag_copied_sets(ctag_names, ctag_vals, 1);
      }
      else{
	i++;
	flags[assm_index]=1;
	dx_orig(assm_index+1, 1)=dxnew[0];
	dx_orig(assm_index+1, 2)=dxnew[1];
	dx_orig(assm_index+1, 3)=dxnew[2];
      }
    }
    dx[0] = 0.0;
    dx[1] = 0.0;
  }

  return iBase_SUCCESS;
}  

int CCrgen::copy_move_hex_flat_assys_p1(CopyMesh **cm,
					const int nrings, const int pack_type,
					const double pitch,
					const int symm,
					std::vector<std::string> &core_alias, 
					std::vector<iBase_EntitySetHandle> &assys)
{
  double dx[3] = {0.0, 0.0, 0.0};
  int assm_index; 
  int i = 0;
  int flags[files.size()];



  for (unsigned int i = 0; i < files.size(); i++) {
    flags[i] = 0;
  }

  for (int n1 = 0; n1 < nrings; n1++) {
    for (int n2 = 0; n2 < n1+1; n2++) {
      err = find_assm(i,assm_index);
      if (-1 == assm_index) {
        i++;
        continue;
      }


      dx[1] = n1 * pitch * sin(PI/3.0)  - n2 * pitch * sin(PI/3.0);
      dx[0] = n1 * pitch *  cos(PI/3.0)  + n2 * pitch * cos(PI/3.0);
      if(flags[assm_index] == 0){

 	move_verts(assys[assm_index], dx);

	std::cout << "Copy/moved A: " << assm_index 
		  <<  " n1=" << n1 << ", n2=" << n2  <<" dX = " <<dx[0]<< " dY = " << dx[1] << std::endl;
      }
      i++;
      flags[assm_index] = 1;
    }
  }
  return iBase_SUCCESS;
}


int CCrgen::copy_move_hex_flat_assys(CopyMesh **cm,
				     const int nrings, const int pack_type,
				     const double pitch,
				     const int symm,
				     std::vector<std::string> &core_alias, 
				     std::vector<iBase_EntitySetHandle> &assys)
{
  double dx[3] = {0.0, 0.0, 0.0};
  iBase_EntityHandle *new_ents;
  int new_ents_alloc, new_ents_size;
  int assm_index; 
  int flags[files.size()];
  CMatrix<double> dx_orig(files.size(), 3);
  dx_orig.Set(0.0);
  int i = 0;

  // get the copy/expand sets
  int num_etags = 3, num_ctags = 1;
  const char *etag_names[] = {"MATERIAL_SET", "DIRICHLET_SET", "NEUMANN_SET"};
  const char *etag_vals[] = {NULL, NULL, NULL};
  const char *ctag_names[] = {"GEOM_DIMENSION"};
  const char *ctag_vals[]={(const char*)&set_DIM};

  for (unsigned int i = 0; i < files.size(); i++) {
    flags[i]=0;
    err = get_expand_sets(cm[i],assys[i], etag_names, etag_vals, num_etags);
    ERRORR("Failed to add expand lists.", iBase_FAILURE);
    err = get_copy_sets(cm[i],assys[i], ctag_names, ctag_vals, num_ctags);
    ERRORR("Failed to add expand lists.", iBase_FAILURE);
    err = extend_expand_sets(cm[i]);
    ERRORR("Failed to extend expand lists.", iBase_FAILURE);
  }

  for (int n1 = 0; n1 < nrings; n1++) {
    for (int n2 = 0; n2 < n1+1; n2++) {
      err = find_assm(i,assm_index);
      if (-1 == assm_index) {
        i++;
        continue;
      }


      dx[1] = n1 * pitch * sin(PI/3.0)  - n2 * pitch * sin(PI/3.0);
      dx[0] = n1 * pitch *  cos(PI/3.0)  + n2 * pitch * cos(PI/3.0);


      if(flags[assm_index]==1){
	dx[0]-=dx_orig(assm_index+1, 1);
	dx[1]-=dx_orig(assm_index+1, 2);
	dx[2]-=dx_orig(assm_index+1, 3);	

	new_ents = NULL;
	new_ents_alloc = 0;
	new_ents_size = 0;

	int orig_ents_alloc = 0, orig_ents_size;
	iBase_EntityHandle *orig_ents = NULL;

	iMesh_getEntities(impl, assys[assm_index], iBase_ALL_TYPES,iMesh_ALL_TOPOLOGIES,
			  &orig_ents, &orig_ents_alloc, &orig_ents_size, &err);
	ERRORR("Failed to get any entities from original set.", iBase_FAILURE);

	cm[assm_index]->copy(orig_ents,orig_ents_size, copy::Translate(dx),
	                     &new_ents, &new_ents_alloc, &new_ents_size, false);
	std::cout << "Copy/moved A: " << assm_index 
		  << " n1=" << n1 << ", n2=" << n2 <<" dX = " <<dx[0]<< " dY = " << dx[1] << std::endl;

	free(new_ents);
	free(orig_ents);
	i++;
    
	cm[assm_index]->tag_copied_sets(ctag_names, ctag_vals, 1);
      }
      else{
	i++;
	flags[assm_index]=1;
	dx_orig(assm_index+1, 1)=dx[0];
	dx_orig(assm_index+1, 2)=dx[1];
	dx_orig(assm_index+1, 3)=dx[2];
      }
    }
  }
  return iBase_SUCCESS;
}


int CCrgen::copy_move_hex_full_assys_p1(CopyMesh **cm,
					const int nrings, const int pack_type,
					const double pitch,
					const int symm,
					std::vector<std::string> &core_alias, 
					std::vector<iBase_EntitySetHandle> &assys)
{
  double dx[3] = {0.0, 0.0, 0.0};
  int  assm_index; 
  int flags[files.size()];
  int i = 0;

  for (unsigned int i = 0; i < files.size(); i++) {
    flags[i] = 0;
  }

  int t, width = 2 * nrings - 1;
  for (int n1 = 1; n1 <= width; n1++) {
    if(n1 > nrings)
      t = 2 * nrings - n1;
    else
      t = n1;

    for (int n2 = 1; n2 <= (nrings + t - 1); n2++) {
      err = find_assm(i,assm_index);
      if (-1 == assm_index) {
        i++;
        continue;
      }

      if (n1 < nrings){
	dx[0] = (nrings - n2 + 1) * pitch / 2.0 + n2 * pitch / 2.0 + 
	  (n2 - 1) * pitch - (n1 - 1) * pitch / 2.0;
	dx[1] = (n1 - 1) * (0.5 * pitch / sin(PI/3.0) + 0.5 * pitch * sin(PI/6.0) / sin(PI/3.0));
      }
      else{
	dx[0] = (nrings - n2 + 1) * pitch / 2.0 + n2 * pitch / 2.0 + (n2 - 1) * pitch -
	  (2 * nrings - n1 -1) * pitch / 2.0;
	dx[1] = (n1 -1) * (0.5 * pitch / sin(PI/3.0) + 0.5 * pitch * sin(PI/6.0) / sin(PI/3.0)); 
      } 
      if(flags[assm_index] == 0){

 	move_verts(assys[assm_index], dx);

	std::cout << "Copy/moved A: " << assm_index 
		  <<  " n1=" << n1 << ", n2=" << n2  <<" dX = " <<dx[0]<< " dY = " << dx[1] << std::endl;
      }
      i++;
      flags[assm_index] = 1;
    }
  }
  return iBase_SUCCESS;
}

int CCrgen::copy_move_hex_full_assys(CopyMesh **cm,
				     const int nrings, const int pack_type,
				     const double pitch,
				     const int symm,
				     std::vector<std::string> &core_alias, 
				     std::vector<iBase_EntitySetHandle> &assys)
{
  double dx[3] = {0.0, 0.0, 0.0};
  iBase_EntityHandle *new_ents;
  int new_ents_alloc, new_ents_size;
  int assm_index; 
  int flags[files.size()];
  CMatrix<double> dx_orig(files.size(), 3);
  dx_orig.Set(0.0);
  int i = 0;

  // get the copy/expand sets
  int num_etags = 3, num_ctags = 1;
  const char *etag_names[] = {"MATERIAL_SET", "DIRICHLET_SET", "NEUMANN_SET"};
  const char *etag_vals[] = {NULL, NULL, NULL};
  const char *ctag_names[] = {"GEOM_DIMENSION"};
  const char *ctag_vals[]={(const char*)&set_DIM};

  for (unsigned int i = 0; i < files.size(); i++) {
    flags[i]=0;
    err = get_expand_sets(cm[i],assys[i], etag_names, etag_vals, num_etags);
    ERRORR("Failed to add expand lists.", iBase_FAILURE);
    err = get_copy_sets(cm[i],assys[i], ctag_names, ctag_vals, num_ctags);
    ERRORR("Failed to add expand lists.", iBase_FAILURE);
    err = extend_expand_sets(cm[i]);
    ERRORR("Failed to extend expand lists.", iBase_FAILURE);
  }

  int t, width = 2 * nrings - 1;
  for (int n1 = 1; n1 <= width; n1++) {
    if(n1 > nrings)
      t = 2 * nrings - n1;
    else
      t = n1;

    for (int n2 = 1; n2 <= (nrings + t - 1); n2++) {
      err = find_assm(i,assm_index);
      if (-1 == assm_index) {
        i++;
        continue;
      }

      if (n1 < nrings){
	dx[0] = (nrings - n2 + 1) * pitch / 2.0 + n2 * pitch / 2.0 + 
	  (n2 - 1) * pitch - (n1 - 1) * pitch / 2.0;
	dx[1] = (n1 - 1) * (0.5 * pitch / sin(PI/3.0) + 0.5 * pitch * sin(PI/6.0) / sin(PI/3.0));
      }
      else{
	dx[0] = (nrings - n2 + 1) * pitch / 2.0 + n2 * pitch / 2.0 + (n2 - 1) * pitch -
	  (2 * nrings - n1 -1) * pitch / 2.0;
	dx[1] = (n1 -1) * (0.5 * pitch / sin(PI/3.0) + 0.5 * pitch * sin(PI/6.0) / sin(PI/3.0)); 
      }     

      if(flags[assm_index]==1){
	dx[0]-=dx_orig(assm_index+1, 1);
	dx[1]-=dx_orig(assm_index+1, 2);
	dx[2]-=dx_orig(assm_index+1, 3);

	new_ents = NULL;
	new_ents_alloc = 0;
	new_ents_size = 0;

	int orig_ents_alloc = 0, orig_ents_size;
	iBase_EntityHandle *orig_ents = NULL;

	iMesh_getEntities(impl, assys[assm_index], iBase_ALL_TYPES,iMesh_ALL_TOPOLOGIES,
			  &orig_ents, &orig_ents_alloc, &orig_ents_size, &err);
	ERRORR("Failed to get any entities from original set.", iBase_FAILURE);

	cm[assm_index]->copy(orig_ents,orig_ents_size, copy::Translate(dx),
	                     &new_ents, &new_ents_alloc, &new_ents_size, false);
	std::cout << "Copy/moved A: " << assm_index 
		  << " n1=" << n1 << ", n2=" << n2 <<" dX = " <<dx[0]<< " dY = " << dx[1] << std::endl;

	free(new_ents);
	free(orig_ents);
	i++;
    
	cm[assm_index]->tag_copied_sets(ctag_names, ctag_vals, 1);
      }
      else{
	i++;
	flags[assm_index]=1;
	dx_orig(assm_index+1, 1)=dx[0];
	dx_orig(assm_index+1, 2)=dx[1];
	dx_orig(assm_index+1, 3)=dx[2];
      }
    }
  }
  return iBase_SUCCESS;
}

int CCrgen::copy_move_sq_assys_p1(CopyMesh **cm,
				  const int nrings, const int pack_type,
				  const double pitch,
				  const int symm,
				  std::vector<std::string> &core_alias,
				  std::vector<iBase_EntitySetHandle> &assys) 
{
  double dx[3] = {0.0, 0.0, 0.0};
  int i = 0, assm_index;

  int flags[files.size()];

  for (unsigned int i = 0; i < files.size(); i++) {
    flags[i] = 0;
  }

  for (int n1 = 0; n1 < nringsx; n1++) {
    dx[1] = n1 * pitchy;
    for (int n2 = 0; n2 < nringsy; n2++) {
      
      err = find_assm(i, assm_index);
      if (-1 == assm_index) {
        i++;
        continue;
      }

      if(flags[assm_index] == 0){

	dx[0] = n2 * pitchx;
 	move_verts(assys[assm_index], dx);

	std::cout << "Copy/moved A: " << assm_index 
		  <<  " n1=" << n1 << ", n2=" << n2  <<" dX = " <<dx[0]<< " dY = " << dx[1] << std::endl;
      }
      i++;
      flags[assm_index] = 1;
    }
  }
  return iBase_SUCCESS;
}


int CCrgen::copy_move_sq_assys(CopyMesh **cm,
			       const int nrings, const int pack_type,
			       const double pitch,
			       const int symm,
			       std::vector<std::string> &core_alias,
			       std::vector<iBase_EntitySetHandle> &assys) 
{

  double dx[3] = {0.0, 0.0, 0.0};
  iBase_EntityHandle *new_ents;
  int new_ents_alloc, new_ents_size;
  int i = 0, assm_index;
  int flags[files.size()];
  CMatrix<double> dx_orig(files.size(), 3);
  dx_orig.Set(0.0);

  // get the copy/expand sets
  int num_etags = 3, num_ctags = 1;
  const char *etag_names[] = {"MATERIAL_SET", "DIRICHLET_SET", "NEUMANN_SET"};
  const char *etag_vals[] = {NULL, NULL, NULL};
  const char *ctag_names[] = {"GEOM_DIMENSION"};
  const char *ctag_vals[]={(const char*)&set_DIM};

  for (unsigned int i = 0; i < files.size(); i++) {
    flags[i]=0;
    err = get_expand_sets(cm[i],assys[i], etag_names, etag_vals, num_etags);
    ERRORR("Failed to add expand lists.", iBase_FAILURE);
    err = get_copy_sets(cm[i],assys[i], ctag_names, ctag_vals, num_ctags);
    ERRORR("Failed to add expand lists.", iBase_FAILURE);
    err = extend_expand_sets(cm[i]);
    ERRORR("Failed to extend expand lists.", iBase_FAILURE);
  }

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

      if(flags[assm_index]==1){
	dx[0]-=dx_orig(assm_index+1, 1);
	dx[1]-=dx_orig(assm_index+1, 2);
	dx[2]-=dx_orig(assm_index+1, 3);
	int orig_ents_alloc = 0, orig_ents_size;
	iBase_EntityHandle *orig_ents = NULL;

	iMesh_getEntities(impl, assys[assm_index], iBase_ALL_TYPES,iMesh_ALL_TOPOLOGIES,
			  &orig_ents, &orig_ents_alloc, &orig_ents_size, &err);
	ERRORR("Failed to get any entities from original set.", iBase_FAILURE);

	cm[assm_index]->copy(orig_ents,orig_ents_size, copy::Translate(dx), 
	                     &new_ents, &new_ents_alloc, &new_ents_size, false);
	std::cout << "Copy/moved A: " << assm_index 
		  <<  " n1=" << n1 << ", n2=" << n2  <<" dX = " <<dx[0]<< " dY = " << dx[1] << std::endl;
	free(new_ents);
	i++;
    
	cm[assm_index]->tag_copied_sets(ctag_names, ctag_vals, 1);
      }
      else{
	i++;
	flags[assm_index]=1;
	dx_orig(assm_index+1, 1)=dx[0];
	dx_orig(assm_index+1, 2)=dx[1];
	dx_orig(assm_index+1, 3)=dx[2];
      }    
    }
  }
  return iBase_SUCCESS;
}

//////////////////////////////////////////////////////////////////////////////////
//////////////// GEOMETRY VERSIONS////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////


int CCrgen::copy_move_hex_vertex_assys_p1(CopyGeom **cg,
					  const int nrings, const int pack_type,
					  const double pitch,
					  const int symm,
					  std::vector<std::string> &core_alias,
					  std::vector<iBase_EntitySetHandle> &assys)
{
  double dx[3] = {0.0, 0.0, 0.0};
  double dxnew[3] = {0.0, 0.0, 0.0};
  int i = 0, bd = 0;
  int assm_index;
  int flags[files.size()];

  for (unsigned int i = 0; i < files.size(); i++) {
    flags[i] = 0;
  }

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

	// starting from x-axis
	dxnew[0] = (dx[0] * cos(PI/6.0) + dx[1] * sin(PI/6.0));
	dxnew[1] = (dx[1] * cos(PI/6.0) - dx[0] * sin(PI/6.0));

	if(flags[assm_index]==0){
	  move_verts(assys[assm_index], dxnew);

	  std::cout << "Copy/moved A: " << assm_index
		    <<  " n1=" << n1 << ", n2=" << n2  <<" dX = " <<dx[0]<< " dY = " << dx[1] << std::endl;
	}
	i++;
	flags[assm_index] = 1;
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


	// starting from x-axis
	dxnew[0] = (dx[0] * cos(PI/6.0) + dx[1] * sin(PI/6.0));
	dxnew[1] = (dx[1] * cos(PI/6.0) - dx[0] * sin(PI/6.0));


	if(flags[assm_index]==0){
    	  move_geoms(assys[assm_index], dxnew);

	  std::cout << "Copy/moved A: " << assm_index
		    <<  " n1=" << n1 << ", n2=" << n2  <<" dX = " <<dx[0]<< " dY = " << dx[1] << std::endl;
	}
	i++;
	flags[assm_index] = 1;
      }
    }
    bd = 0;
  }
  return iBase_SUCCESS;
}


int CCrgen::copy_move_hex_vertex_assys(CopyGeom **cg,
				       const int nrings, const int pack_type,
				       const double pitch,
				       const int symm,
				       std::vector<std::string> &core_alias,
				       std::vector<iBase_EntitySetHandle> &assys)
{
  double dx[3] = {0.0, 0.0, 0.0};
  double dxnew[3] = {0.0, 0.0, 0.0};
  iBase_EntityHandle *new_ents;
  int new_ents_alloc, new_ents_size;
  int i = 0, bd = 0;
  int assm_index;

  int flags[files.size()];
  CMatrix<double> dx_orig(files.size(), 3);
  dx_orig.Set(0.0);

  // get the copy/expand sets
  int num_etags = 3, num_ctags = 1;
  const char *etag_names[] = {"MATERIAL_SET", "DIRICHLET_SET", "NEUMANN_SET"};
  const char *etag_vals[] = {NULL, NULL, NULL};
  const char *ctag_names[] = {"GEOM_DIMENSION"};
  const char *ctag_vals[]={(const char*)&set_DIM};

  for (unsigned int i = 0; i < files.size(); i++) {
    flags[i]=0;
    err = get_expand_sets(cg[i],assys[i], etag_names, etag_vals, num_etags);
    ERRORR("Failed to add expand lists.", iBase_FAILURE);
    err = get_copy_sets(cg[i],assys[i], ctag_names, ctag_vals, num_ctags);
    ERRORR("Failed to add expand lists.", iBase_FAILURE);
    err = extend_expand_sets(cg[i]);
    ERRORR("Failed to extend expand lists.", iBase_FAILURE);
  }

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

 	// starting from x-axis
 	dxnew[0] = (dx[0] * cos(PI/6.0) + dx[1] * sin(PI/6.0));
 	dxnew[1] = (dx[1] * cos(PI/6.0) - dx[0] * sin(PI/6.0));

 	if(flags[assm_index]==1){
 	  dxnew[0]-=dx_orig(assm_index+1, 1);
 	  dxnew[1]-=dx_orig(assm_index+1, 2);
 	  dxnew[2]-=dx_orig(assm_index+1, 3);

 	  int orig_ents_alloc = 0, orig_ents_size;
 	  iBase_EntityHandle *orig_ents = NULL;

 	  iGeom_getEntities(geom, assys[assm_index], iBase_ALL_TYPES, &orig_ents, &orig_ents_alloc, &orig_ents_size, &err);
 	  ERRORR("Failed to get any entities from original set.", iBase_FAILURE);


 	  cg[assm_index]->copy(orig_ents,orig_ents_size, dxnew,
 	                       &new_ents, &new_ents_alloc, &new_ents_size, false);
 	  std::cout << "Copy/moved A: " << assm_index
 		    << " n1=" << n1 << ", n2=" << n2 <<" dX = " <<dxnew[0]<< " dY = " << dxnew[1] << std::endl;
 	  free(new_ents);
 	  free(orig_ents);
 	  i++;

 	  cg[assm_index]->tag_copied_sets(ctag_names, ctag_vals, 1);
 	}
 	else{
 	  i++;
 	  flags[assm_index]=1;
 	  dx_orig(assm_index+1, 1)=dxnew[0];
 	  dx_orig(assm_index+1, 2)=dxnew[1];
 	  dx_orig(assm_index+1, 3)=dxnew[2];
 	}
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


 	// starting from x-axis
 	dxnew[0] = (dx[0] * cos(PI/6.0) + dx[1] * sin(PI/6.0));
 	dxnew[1] = (dx[1] * cos(PI/6.0) - dx[0] * sin(PI/6.0));


 	if(flags[assm_index]==1){
 	  dxnew[0]-=dx_orig(assm_index+1, 1);
 	  dxnew[1]-=dx_orig(assm_index+1, 2);
 	  dxnew[2]-=dx_orig(assm_index+1, 3);

 	  int orig_ents_alloc = 0, orig_ents_size;
 	  iBase_EntityHandle *orig_ents = NULL;

 	  iGeom_getEntities(geom, assys[assm_index], iBase_ALL_TYPES, &orig_ents, &orig_ents_alloc, &orig_ents_size, &err);
 	  ERRORR("Failed to get any entities from original set.", iBase_FAILURE);

 	  cg[assm_index]->copy(orig_ents,orig_ents_size, dxnew,
 	                       &new_ents, &new_ents_alloc, &new_ents_size, false);
 	  std::cout << "Copy/moved A: " << assm_index
 		    << " n1=" << n1 << ", n2=" << n2 <<" dX = " <<dxnew[0]<< " dY = " << dxnew[1] << std::endl;
 	  free(new_ents);
 	  free(orig_ents);
 	  i++;

 	  cg[assm_index]->tag_copied_sets(ctag_names, ctag_vals, 1);
 	}
 	else{
 	  i++;
 	  flags[assm_index]=1;
 	  dx_orig(assm_index+1, 1)=dxnew[0];
 	  dx_orig(assm_index+1, 2)=dxnew[1];
 	  dx_orig(assm_index+1, 3)=dxnew[2];
 	}

      }
    }
    bd = 0;
  }
  return iBase_SUCCESS;
}
int CCrgen::copy_move_one_twelfth_assys_p1(CopyGeom **cg,
					   const int nrings, const int pack_type,
					   const double pitch,
					   const int symm,
					   std::vector<std::string> &core_alias,
					   std::vector<iBase_EntitySetHandle> &assys)
{
  double dx[3] = {0.0, 0.0, 0.0};
  double dxnew[3] = {0.0, 0.0, 0.0};
  int i = 0, flag = 0;
  int flags[files.size()];
  int assm_index;

  for (unsigned int i = 0; i < files.size(); i++) {
    flags[i] = 0;
  }


  for (int n1 = 0; n1 < nrings; n1++) {
    int loc = (n1 + 2)/2;

    if( flag == 0 ){
      dx[0] = (n1 + loc - 1) * pitch / 2.0;
      dx[1] = (n1 - loc + 1) * pitch * sin(PI/3.0);
      flag = 1;
    }
    else{
      dx[0] = (n1 + loc) * pitch / 2.0;
      dx[1] = (n1 - loc) * pitch * sin(PI/3.0);
      flag = 0;
    }

    for (int n2 = 0; n2 < loc; n2++) {
      err = find_assm(i,assm_index);
      if (-1 == assm_index) {
	i++;
	continue;
      }

      dxnew[1] = dx[1] - n2 * pitch * sin(PI/3.0);
      dxnew[0] = dx[0] + n2 * pitch * cos(PI/3.0);
      if(flags[assm_index] == 0){

	move_geoms(assys[assm_index], dxnew);

	std::cout << "Copy/moved A: " << assm_index
		  <<  " n1=" << n1 << ", n2=" << n2  <<" dX = " <<dxnew[0]<< " dY = " << dxnew[1] << std::endl;
      }
      i++;
      flags[assm_index] = 1;
    }
    dx[0] = 0.0;
    dx[1] = 0.0;
  }

  return iBase_SUCCESS;
}

int CCrgen::copy_move_one_twelfth_assys(CopyGeom **cg,
					const int nrings, const int pack_type,
					const double pitch,
					const int symm,
					std::vector<std::string> &core_alias,
					std::vector<iBase_EntitySetHandle> &assys)
{
  double dx[3] = {0.0, 0.0, 0.0};
  double dxnew[3] = {0.0, 0.0, 0.0};
  iBase_EntityHandle *new_ents;
  int new_ents_alloc, new_ents_size;
  int i = 0, flag = 0;
  int assm_index;
  int flags[files.size()];
  CMatrix<double> dx_orig(files.size(), 3);
  dx_orig.Set(0.0);

  // get the copy/expand sets
  int num_etags = 3, num_ctags = 1;
  const char *etag_names[] = {"MATERIAL_SET", "DIRICHLET_SET", "NEUMANN_SET"};
  const char *etag_vals[] = {NULL, NULL, NULL};
  const char *ctag_names[] = {"GEOM_DIMENSION"};
  const char *ctag_vals[]={(const char*)&set_DIM};

  for (unsigned int i = 0; i < files.size(); i++) {
    flags[i]=0;
    err = get_expand_sets(cg[i],assys[i], etag_names, etag_vals, num_etags);
    ERRORR("Failed to add expand lists.", iBase_FAILURE);
    err = get_copy_sets(cg[i],assys[i], ctag_names, ctag_vals, num_ctags);
    ERRORR("Failed to add expand lists.", iBase_FAILURE);
    err = extend_expand_sets(cg[i]);
    ERRORR("Failed to extend expand lists.", iBase_FAILURE);
  }


  for (int n1 = 0; n1 < nrings; n1++) {
    int loc = (n1 + 2)/2;

    if( flag == 0 ){
      dx[0] = (n1 + loc - 1) * pitch / 2.0;
      dx[1] = (n1 - loc + 1) * pitch * sin(PI/3.0);
      flag = 1;
    }
    else{
      dx[0] = (n1 + loc) * pitch / 2.0;
      dx[1] = (n1 - loc) * pitch * sin(PI/3.0);
      flag = 0;
    }

    for (int n2 = 0; n2 < loc; n2++) {
      err = find_assm(i,assm_index);
      if (-1 == assm_index) {
	i++;
	continue;
      }

      dxnew[1] = dx[1] - n2 * pitch * sin(PI/3.0);
      dxnew[0] = dx[0] + n2 * pitch * cos(PI/3.0);

      if(flags[assm_index]==1){
	dxnew[0]-=dx_orig(assm_index+1, 1);
	dxnew[1]-=dx_orig(assm_index+1, 2);
	dxnew[2]-=dx_orig(assm_index+1, 3);
	new_ents = NULL;
	new_ents_alloc = 0;
	new_ents_size = 0;

	int orig_ents_alloc = 0, orig_ents_size;
	iBase_EntityHandle *orig_ents = NULL;

	iGeom_getEntities(geom, assys[assm_index], iBase_ALL_TYPES, &orig_ents, &orig_ents_alloc, &orig_ents_size, &err);
	ERRORR("Failed to get any entities from original set.", iBase_FAILURE);

	cg[assm_index]->copy(orig_ents,orig_ents_size, dxnew,
	                     &new_ents, &new_ents_alloc, &new_ents_size, false);
	std::cout << "Copy/moved A: " << assm_index
		  << " n1=" << n1 << ", n2=" << n2 <<" dX = " <<dxnew[0]<< " dY = " << dxnew[1] << std::endl;

	free(new_ents);
	free(orig_ents);
	i++;
	dxnew[0] = 0.0;
	dxnew[1] = 0.0;
	cg[assm_index]->tag_copied_sets(ctag_names, ctag_vals, 1);
      }
      else{
	i++;
	flags[assm_index]=1;
	dx_orig(assm_index+1, 1)=dxnew[0];
	dx_orig(assm_index+1, 2)=dxnew[1];
	dx_orig(assm_index+1, 3)=dxnew[2];
      }
    }
    dx[0] = 0.0;
    dx[1] = 0.0;
  }

  return iBase_SUCCESS;
}

int CCrgen::copy_move_hex_flat_assys_p1(CopyGeom **cg,
					const int nrings, const int pack_type,
					const double pitch,
					const int symm,
					std::vector<std::string> &core_alias,
					std::vector<iBase_EntitySetHandle> &assys)
{
  double dx[3] = {0.0, 0.0, 0.0};
  int assm_index;
  int i = 0;
  int flags[files.size()];



  for (unsigned int i = 0; i < files.size(); i++) {
    flags[i] = 0;
  }

  for (int n1 = 0; n1 < nrings; n1++) {
    for (int n2 = 0; n2 < n1+1; n2++) {
      err = find_assm(i,assm_index);
      if (-1 == assm_index) {
        i++;
        continue;
      }


      dx[1] = n1 * pitch * sin(PI/3.0)  - n2 * pitch * sin(PI/3.0);
      dx[0] = n1 * pitch *  cos(PI/3.0)  + n2 * pitch * cos(PI/3.0);
      if(flags[assm_index] == 0){

	move_geoms(assys[assm_index], dx);
	std::cout << "Copy/moved A: " << assm_index
		  <<  " n1=" << n1 << ", n2=" << n2  <<" dX = " <<dx[0]<< " dY = " << dx[1] << std::endl;
      }
      i++;
      flags[assm_index] = 1;
    }
  }
  return iBase_SUCCESS;
}


int CCrgen::copy_move_hex_flat_assys(CopyGeom **cg,
				     const int nrings, const int pack_type,
				     const double pitch,
				     const int symm,
				     std::vector<std::string> &core_alias,
				     std::vector<iBase_EntitySetHandle> &assys)
{
  double dx[3] = {0.0, 0.0, 0.0};
  iBase_EntityHandle *new_ents;
  int new_ents_alloc, new_ents_size;
  int assm_index;
  int flags[files.size()];
  CMatrix<double> dx_orig(files.size(), 3);
  dx_orig.Set(0.0);
  int i = 0;

  // get the copy/expand sets
  int num_etags = 3, num_ctags = 1;
  const char *etag_names[] = {"MATERIAL_SET", "DIRICHLET_SET", "NEUMANN_SET"};
  const char *etag_vals[] = {NULL, NULL, NULL};
  const char *ctag_names[] = {"GEOM_DIMENSION"};
  const char *ctag_vals[]={(const char*)&set_DIM};

  for (unsigned int i = 0; i < files.size(); i++) {
    flags[i]=0;
    err = get_expand_sets(cg[i],assys[i], etag_names, etag_vals, num_etags);
    ERRORR("Failed to add expand lists.", iBase_FAILURE);
    err = get_copy_sets(cg[i],assys[i], ctag_names, ctag_vals, num_ctags);
    ERRORR("Failed to add expand lists.", iBase_FAILURE);
    err = extend_expand_sets(cg[i]);
  }


  for (int n1 = 0; n1 < nrings; n1++) {
    for (int n2 = 0; n2 < n1+1; n2++) {
      err = find_assm(i,assm_index);
      if (-1 == assm_index) {
        i++;
        continue;
      }


      dx[1] = n1 * pitch * sin(PI/3.0)  - n2 * pitch * sin(PI/3.0);
      dx[0] = n1 * pitch *  cos(PI/3.0)  + n2 * pitch * cos(PI/3.0);


      if(flags[assm_index]==1){
	dx[0]-=dx_orig(assm_index+1, 1);
	dx[1]-=dx_orig(assm_index+1, 2);
	dx[2]-=dx_orig(assm_index+1, 3);

	new_ents = NULL;
	new_ents_alloc = 0;
	new_ents_size = 0;

	// copy/move the geometry
	int entities_ehsize =0, entities_ehallocated =0;
	iBase_EntityHandle *entities = NULL;

	iGeom_getEntities( geom, assys[assm_index], iBase_ALL_TYPES,
			   &entities, &entities_ehallocated, &entities_ehsize, &err );
	ERRORR("Failed to get entities from set.", iBase_FAILURE);

	cg[assm_index]->copy(entities, entities_ehsize, dx, &new_ents, &new_ents_alloc, &new_ents_size);
	ERRORR("Failed to get entities from set.", iBase_FAILURE);

	std::cout << "Copy/moved A: " << assm_index
		  << " n1=" << n1 << ", n2=" << n2 <<" dX = " <<dx[0]<< " dY = " << dx[1] << std::endl;

	free(new_ents);
	free(entities);
	i++;

	cg[assm_index]->tag_copied_sets(ctag_names, ctag_vals, 1);
	ERRORR("Failed to get entities from set.", iBase_FAILURE);

      }
      else{
	i++;
	flags[assm_index]=1;
	dx_orig(assm_index+1, 1)=dx[0];
	dx_orig(assm_index+1, 2)=dx[1];
	dx_orig(assm_index+1, 3)=dx[2];
      }
    }
  }
  return iBase_SUCCESS;
}


int CCrgen::copy_move_hex_full_assys_p1(CopyGeom **cg,
					const int nrings, const int pack_type,
					const double pitch,
					const int symm,
					std::vector<std::string> &core_alias,
					std::vector<iBase_EntitySetHandle> &assys)
{
  double dx[3] = {0.0, 0.0, 0.0};
  int  assm_index;
  int flags[files.size()];
  int i = 0;

  for (unsigned int i = 0; i < files.size(); i++) {
    flags[i] = 0;
  }

  int t, width = 2 * nrings - 1;
  for (int n1 = 1; n1 <= width; n1++) {
    if(n1 > nrings)
      t = 2 * nrings - n1;
    else
      t = n1;

    for (int n2 = 1; n2 <= (nrings + t - 1); n2++) {
      err = find_assm(i,assm_index);
      if (-1 == assm_index) {
	i++;
	continue;
      }

      if (n1 < nrings){
	dx[0] = (nrings - n2 + 1) * pitch / 2.0 + n2 * pitch / 2.0 +
	  (n2 - 1) * pitch - (n1 - 1) * pitch / 2.0;
	dx[1] = (n1 - 1) * (0.5 * pitch / sin(PI/3.0) + 0.5 * pitch * sin(PI/6.0) / sin(PI/3.0));
      }
      else{
	dx[0] = (nrings - n2 + 1) * pitch / 2.0 + n2 * pitch / 2.0 + (n2 - 1) * pitch -
	  (2 * nrings - n1 -1) * pitch / 2.0;
	dx[1] = (n1 -1) * (0.5 * pitch / sin(PI/3.0) + 0.5 * pitch * sin(PI/6.0) / sin(PI/3.0));
      }
      if(flags[assm_index] == 0){

	move_geoms(assys[assm_index], dx);

	std::cout << "Copy/moved A: " << assm_index
		  <<  " n1=" << n1 << ", n2=" << n2  <<" dX = " <<dx[0]<< " dY = " << dx[1] << std::endl;
      }
      i++;
      flags[assm_index] = 1;
    }
  }
  return iBase_SUCCESS;
}

int CCrgen::copy_move_hex_full_assys(CopyGeom **cg,
				     const int nrings, const int pack_type,
				     const double pitch,
				     const int symm,
				     std::vector<std::string> &core_alias,
				     std::vector<iBase_EntitySetHandle> &assys)
{
  double dx[3] = {0.0, 0.0, 0.0};
  iBase_EntityHandle *new_ents;
  int new_ents_alloc, new_ents_size;
  int assm_index;
  int flags[files.size()];
  CMatrix<double> dx_orig(files.size(), 3);
  dx_orig.Set(0.0);
  int i = 0;

  // get the copy/expand sets
  int num_etags = 3, num_ctags = 1;
  const char *etag_names[] = {"MATERIAL_SET", "DIRICHLET_SET", "NEUMANN_SET"};
  const char *etag_vals[] = {NULL, NULL, NULL};
  const char *ctag_names[] = {"GEOM_DIMENSION"};
  const char *ctag_vals[]={(const char*)&set_DIM};

  for (unsigned int i = 0; i < files.size(); i++) {
    flags[i]=0;
    err = get_expand_sets(cg[i],assys[i], etag_names, etag_vals, num_etags);
    ERRORR("Failed to add expand lists.", iBase_FAILURE);
    err = get_copy_sets(cg[i],assys[i], ctag_names, ctag_vals, num_ctags);
    ERRORR("Failed to add expand lists.", iBase_FAILURE);
    err = extend_expand_sets(cg[i]);
    ERRORR("Failed to extend expand lists.", iBase_FAILURE);
  }

  int t, width = 2 * nrings - 1;
  for (int n1 = 1; n1 <= width; n1++) {
    if(n1 > nrings)
      t = 2 * nrings - n1;
    else
      t = n1;

    for (int n2 = 1; n2 <= (nrings + t - 1); n2++) {
      err = find_assm(i,assm_index);
      if (-1 == assm_index) {
	i++;
	continue;
      }

      if (n1 < nrings){
	dx[0] = (nrings - n2 + 1) * pitch / 2.0 + n2 * pitch / 2.0 +
	  (n2 - 1) * pitch - (n1 - 1) * pitch / 2.0;
	dx[1] = (n1 - 1) * (0.5 * pitch / sin(PI/3.0) + 0.5 * pitch * sin(PI/6.0) / sin(PI/3.0));
      }
      else{
	dx[0] = (nrings - n2 + 1) * pitch / 2.0 + n2 * pitch / 2.0 + (n2 - 1) * pitch -
	  (2 * nrings - n1 -1) * pitch / 2.0;
	dx[1] = (n1 -1) * (0.5 * pitch / sin(PI/3.0) + 0.5 * pitch * sin(PI/6.0) / sin(PI/3.0));
      }

      if(flags[assm_index]==1){
	dx[0]-=dx_orig(assm_index+1, 1);
	dx[1]-=dx_orig(assm_index+1, 2);
	dx[2]-=dx_orig(assm_index+1, 3);

	new_ents = NULL;
	new_ents_alloc = 0;
	new_ents_size = 0;

	int orig_ents_alloc = 0, orig_ents_size;
	iBase_EntityHandle *orig_ents = NULL;

	iGeom_getEntities(geom, assys[assm_index], iBase_ALL_TYPES, &orig_ents, &orig_ents_alloc, &orig_ents_size, &err);
	ERRORR("Failed to get any entities from original set.", iBase_FAILURE);

	cg[assm_index]->copy(orig_ents,orig_ents_size, dx,
	                     &new_ents, &new_ents_alloc, &new_ents_size, false);
	std::cout << "Copy/moved A: " << assm_index
		  << " n1=" << n1 << ", n2=" << n2 <<" dX = " <<dx[0]<< " dY = " << dx[1] << std::endl;

	free(new_ents);
	free(orig_ents);
	i++;

	cg[assm_index]->tag_copied_sets(ctag_names, ctag_vals, 1);
      }
      else{
	i++;
	flags[assm_index]=1;
	dx_orig(assm_index+1, 1)=dx[0];
	dx_orig(assm_index+1, 2)=dx[1];
	dx_orig(assm_index+1, 3)=dx[2];
      }
    }
  }
  return iBase_SUCCESS;
}

int CCrgen::copy_move_sq_assys_p1(CopyGeom **cg,
				  const int nrings, const int pack_type,
				  const double pitch,
				  const int symm,
				  std::vector<std::string> &core_alias,
				  std::vector<iBase_EntitySetHandle> &assys)
{
  double dx[3] = {0.0, 0.0, 0.0};
  int i = 0, assm_index;

  int flags[files.size()];

  for (unsigned int i = 0; i < files.size(); i++) {
    flags[i] = 0;
  }

  for (int n1 = 0; n1 < nringsx; n1++) {
    dx[1] = n1 * pitchy;
    for (int n2 = 0; n2 < nringsy; n2++) {

      err = find_assm(i, assm_index);
      if (-1 == assm_index) {
	i++;
	continue;
      }

      if(flags[assm_index] == 0){

	dx[0] = n2 * pitchx;
	move_geoms(assys[assm_index], dx);

	std::cout << "Copy/moved A: " << assm_index
		  <<  " n1=" << n1 << ", n2=" << n2  <<" dX = " <<dx[0]<< " dY = " << dx[1] << std::endl;
      }
      i++;
      flags[assm_index] = 1;
    }
  }
  return iBase_SUCCESS;
}


int CCrgen::copy_move_sq_assys(CopyGeom **cg,
			       const int nrings, const int pack_type,
			       const double pitch,
			       const int symm,
			       std::vector<std::string> &core_alias,
			       std::vector<iBase_EntitySetHandle> &assys)
{

  double dx[3] = {0.0, 0.0, 0.0};
  iBase_EntityHandle *new_ents;
  int new_ents_alloc, new_ents_size;
  int i = 0, assm_index;
  int flags[files.size()];
  CMatrix<double> dx_orig(files.size(), 3);
  dx_orig.Set(0.0);

  // get the copy/expand sets
  int num_etags = 3, num_ctags = 1;
  const char *etag_names[] = {"MATERIAL_SET", "DIRICHLET_SET", "NEUMANN_SET"};
  const char *etag_vals[] = {NULL, NULL, NULL};
  const char *ctag_names[] = {"GEOM_DIMENSION"};
  const char *ctag_vals[]={(const char*)&set_DIM};

  for (unsigned int i = 0; i < files.size(); i++) {
    flags[i]=0;
    err = get_expand_sets(cg[i],assys[i], etag_names, etag_vals, num_etags);
    ERRORR("Failed to add expand lists.", iBase_FAILURE);
    err = get_copy_sets(cg[i],assys[i], ctag_names, ctag_vals, num_ctags);
    ERRORR("Failed to add expand lists.", iBase_FAILURE);
    err = extend_expand_sets(cg[i]);
    ERRORR("Failed to extend expand lists.", iBase_FAILURE);
  }

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

      if(flags[assm_index]==1){
	dx[0]-=dx_orig(assm_index+1, 1);
	dx[1]-=dx_orig(assm_index+1, 2);
	dx[2]-=dx_orig(assm_index+1, 3);
	int orig_ents_alloc = 0, orig_ents_size;
	iBase_EntityHandle *orig_ents = NULL;

	iGeom_getEntities(geom, assys[assm_index], iBase_ALL_TYPES, &orig_ents, &orig_ents_alloc, &orig_ents_size, &err);
	ERRORR("Failed to get any entities from original set.", iBase_FAILURE);

	cg[assm_index]->copy(orig_ents,orig_ents_size, dx,
	                     &new_ents, &new_ents_alloc, &new_ents_size, false);
	std::cout << "Copy/moved A: " << assm_index
		  <<  " n1=" << n1 << ", n2=" << n2  <<" dX = " <<dx[0]<< " dY = " << dx[1] << std::endl;
	free(new_ents);
	i++;

	cg[assm_index]->tag_copied_sets(ctag_names, ctag_vals, 1);
      }
      else{
	i++;
	flags[assm_index]=1;
	dx_orig(assm_index+1, 1)=dx[0];
	dx_orig(assm_index+1, 2)=dx[1];
	dx_orig(assm_index+1, 3)=dx[2];
      }
    }
  }
  return iBase_SUCCESS;
}


