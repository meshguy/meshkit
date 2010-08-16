#include "CopyMesh.hpp"
#include "MergeMesh.hpp"
#include "TestFramework.hpp"
#include <cfloat>

iMesh_Instance impl;
iBase_EntitySetHandle root_set;

int check_num_ents(int ent_type, int expected) 
{
  int num_type, err;
  iMesh_getNumOfType(impl, root_set, ent_type, &num_type, &err);
  CHECK_ERR("Failed to get # entities");
  if (num_type != expected) {
    std::cerr << "Didn't get right # entities of dimension " << ent_type
              << "; got " << num_type
              << ", expected " << expected << std::endl;
    return 1;
  }
  return 0;
}

int parse_input(int argc, char **argv, 
                std::string &infile, std::string &outfile, 
                double *dx, 
                std::vector<std::string> &ctags, 
                std::vector<std::string> &etags, 
                std::vector<std::string> &utags,
                std::vector<char*> &cvals, 
                std::vector<char*> &evals);

int make_ents(CopyMesh *cm, iBase_EntityHandle *&ents, 
              int &ents_alloc, int &ents_size, 
              iBase_EntityHandle *&verts, 
              int &verts_alloc, int &verts_size);

int main(int argc, char **argv) 
{
  // make a mesh instance
  int err;
  iMesh_newMesh("MOAB", &impl, &err, 4);
  CHECK_ERR("Failed to create instance.");
  
  iMesh_getRootSet(impl, &root_set, &err);

  CopyMesh *cm = new CopyMesh(impl);
  MergeMesh *mm = new MergeMesh(impl);

  // some entity tag types are always copy or expand
  CHECK_THROW( cm->expand_sets().add_tag("MATERIAL_SET") );
  CHECK_THROW( cm->expand_sets().add_tag("DIRICHLET_SET") );
  CHECK_THROW( cm->expand_sets().add_tag("NEUMANN_SET") );
  
  if (1 == argc) {
    int tag_val = 3;
    CHECK_THROW( cm->copy_sets().add_tag(std::string("GEOM_DIMENSION"),
                                         (const char*) &tag_val) );

    // create vertices
    iBase_EntityHandle *ents = NULL, *verts = NULL;
    int ents_alloc = 0, ents_size, verts_alloc = 0, verts_size, 
      new_ents_alloc, new_ents_size;

    err = make_ents(cm, ents, ents_alloc, ents_size, 
                    verts, verts_alloc, verts_size);
    CHECK_ERR("Making entities failed.");
    
    double dx[] = {1.0,0.0,0.0};
    iBase_EntityHandle *new_ents = ents + ents_size;
    new_ents_alloc = ents_size;
    
    CHECK_THROW( cm->copy(ents, ents_size, copy::Translate(dx), &new_ents,
                          &new_ents_alloc, &new_ents_size, false) );

    // check # entities
    if (check_num_ents(iBase_VERTEX, 32)) return 1;
    if (check_num_ents(iBase_REGION, 6)) return 1;
  
    CHECK_THROW( mm->merge_entities(ents, ents_size+new_ents_size, 1.0e-8) );

    // now get all vertices, put in new verts array
    iBase_EntityHandle *nverts = NULL;
    int nverts_alloc = 0, nverts_size;
    iMesh_getEntities(impl, root_set, iBase_VERTEX, iMesh_ALL_TOPOLOGIES, 
                      &nverts, &nverts_alloc, &nverts_size, &err);
    CHECK_ERR("Didn't get all vertices.");

    // and their coords
    double *vert_coords = NULL;
    int vert_coords_alloc = 0, vert_coords_size;
    int sorder = iBase_INTERLEAVED;
    iMesh_getVtxArrCoords (impl, nverts, nverts_size, sorder, 
                           &vert_coords, &vert_coords_alloc, &vert_coords_size,
                           &err);
    CHECK_ERR("Didn't get vtx coords.");

    double xmin[3] = {DBL_MAX, DBL_MAX, DBL_MAX};
    double xmax[3] = {DBL_MIN, DBL_MIN, DBL_MIN};

    for (int i = 0; i < nverts_size; i++) {
      for (int j = 0; j < 3; j++) {
        if (vert_coords[3*i+j] < xmin[j]) xmin[j] = vert_coords[3*i+j];
        if (vert_coords[3*i+j] > xmax[j]) xmax[j] = vert_coords[3*i+j];
      }
    }
    if (xmin[0] != 0.0 || xmin[1] != 0.0 || xmin[2] != 0.0 ||
        xmax[0] != 2.0 || xmax[1] != 1.0 || xmax[2] != 3.0) {
      std::cerr << "Didn't get correct min/max; output values = ("
                << xmin[0] << ", " << xmin[1] << ", " << xmin[2] << "), ("
                << xmax[0] << ", " << xmax[1] << ", " << xmax[2] << ")" << std::endl;
      return 1;
    }

    int expected = 2*verts_size - 8;
    if (nverts_size != expected) {
      std::cerr << "Didn't get right # vertices; got " << nverts_size
                << ", expected " << expected << std::endl;
    
      return 1;
    }
  
    free(verts);
    free(nverts);
    free(ents);
  }
  
  else {
    std::vector<std::string> ctags, etags, utags;
    std::vector<char*> cvals, evals;
    std::string infile, outfile;
    double dx[3] = {0.0, 0.0, 0.0};

    err = parse_input(argc, argv, 
                      infile, outfile, dx, 
                      ctags, etags, utags,
                      cvals, evals);
    if (-1 == err) return err;
    CHECK_ERR("Couldn't parse input, exiting.");
    
    // read the file
    iMesh_load(impl, 0, infile.c_str(), NULL, &err, infile.length(), 0);
    CHECK_ERR("Couldn't read mesh file.");
    
    // set copy/expand/unique tags
    for (unsigned int i = 0; i < ctags.size(); i++)
      CHECK_THROW( cm->copy_sets().add_tag(ctags[i], cvals[i]) );
    for (unsigned int i = 0; i < etags.size(); i++)
      CHECK_THROW( cm->expand_sets().add_tag(etags[i], evals[i]) );
    for (unsigned int i = 0; i < utags.size(); i++)
      CHECK_THROW( cm->add_unique_tag(utags[i]) );
      
    // copy
    iBase_EntityHandle* new_ents = NULL;
    int new_ents_alloc = 0;
    int new_ents_size = 0;
    CHECK_THROW( cm->copy(root_set, copy::Translate(dx), &new_ents,
                          &new_ents_alloc, &new_ents_size,
                          false) );
    //getting elements for merge mesh   
    iBase_EntityHandle *ents = NULL;
    int ents_alloc = 0, ents_size;
    iMesh_getEntitiesRec(impl, root_set, 
			 iBase_ALL_TYPES, iMesh_HEXAHEDRON, true,
			 &ents, &ents_alloc, &ents_size, &err);
    CHECK_ERR("Failed to get entities from set recursively.");
    // merge 
    const double merge_tol =  1.0e-8;
    const int do_merge = 0;
    const int update_sets= 0;
    iBase_TagHandle merge_tag = NULL;
    CHECK_THROW( mm->merge_entities(ents, ents_size, merge_tol,
                                    do_merge, update_sets, merge_tag) );

    // export
    iMesh_save(impl, root_set, outfile.c_str(), NULL, &err, 
               outfile.length(), 0);
    CHECK_ERR("Failed to save mesh.");
    std::cout << "Wrote " << outfile << std::endl;
  }

  delete cm;
  delete mm;
  
  iMesh_dtor(impl, &err);
  CHECK_ERR("Destructor failed.");
  
  return 0;
}

int parse_input(int argc, char **argv, 
                std::string &infile, std::string &outfile, 
                double *dx, 
                std::vector<std::string> &ctags, 
                std::vector<std::string> &etags, 
                std::vector<std::string> &utags,
                std::vector<char*> &cvals, 
                std::vector<char*> &evals) 
{
  if (!strcmp(argv[1], "-h")) {
    std::cout << "Usage: " << argv[0] << " _OR_ " << std::endl
              << argv[0] << " [-c copy_tag [-cv copy_val]]* "
              << "[-e expand_tag [-ev copy_val]]* "
              << "[-u unique_tag]* " 
              << "[-x <dx> <dy> <dz>] <infile> <outfile>" << std::endl;
    std::cout << "* - repeated zero or more times" << std::endl;
    return -1;
  }
  
  int pos = 1;
  for (; pos < argc; pos++)
    {
      if (!strcmp(argv[pos],"-c")) {
	pos++;
	ctags.push_back(argv[pos]);
	if (!strcmp(argv[pos],"-cv")) {
	  pos++;
	  cvals.push_back(argv[pos]);
	}
	else cvals.push_back(NULL);
      }
      else if (!strcmp(argv[pos],"-e")) {
	pos++;
	etags.push_back(argv[pos]);
	if (!strcmp(argv[pos],"-ev")) {
	  pos++;
	  evals.push_back(argv[pos]);
	}
	else evals.push_back(NULL);
      }
      else if (!strcmp(argv[pos],"-u")) {
	pos++;
	utags.push_back(argv[pos]);
      }
      else if (!strcmp(argv[pos],"-x")) {
	pos++;
	dx[0] = atof(argv[pos++]);
	dx[1] = atof(argv[pos++]);
	dx[2] = atof(argv[pos]);
      }
      else if (argv[pos][0] != '-') {
	break;
      }
      else {
	std::cerr << "Invalid option: \"" << argv[pos] << '"' << std::endl;
      }
    }

  if (pos == argc) {
    std::cerr << "Not enough arguments, infile missing." << std::endl;
    return 1;
  }
  else
    infile = argv[pos++];

  if (pos == argc) {
    std::cerr << "Not enough arguments, outfile missing." << std::endl;
    return 1;
  }
  
  outfile = argv[pos++];
  
  return 0;
}

int make_ents(CopyMesh *cm, iBase_EntityHandle *&ents, 
              int &ents_alloc, int &ents_size, 
              iBase_EntityHandle *&verts, 
              int &verts_alloc, int &verts_size) 
{
  // make an arbitrary mesh to copy
  double vert_pos[] = {
    0.0, 0.0, 0.0,
    1.0, 0.0, 0.0,
    1.0, 1.0, 0.0,
    0.0, 1.0, 0.0,
    0.0, 0.0, 1.0,
    1.0, 0.0, 1.0,
    1.0, 1.0, 1.0,
    0.0, 1.0, 1.0,
    0.0, 0.0, 2.0,
    1.0, 0.0, 2.0,
    1.0, 1.0, 2.0,
    0.0, 1.0, 2.0,
    0.0, 0.0, 3.0,
    1.0, 0.0, 3.0,
    1.0, 1.0, 3.0,
    0.0, 1.0, 3.0
  };
	 
  int connect_idx[] = {
    0, 1, 2, 3, 4, 5, 6, 7,
    4, 5, 6, 7, 8, 9, 10, 11,
    8, 9, 10, 11, 12, 13, 14, 15
  };

  int adj_table[] = {
    1, 0, 0, 1,
    1, 0, 0, 0,
    1, 0, 1, 0, 
    1, 0, 0, 1
  };
  int err;

  // set the adjacency table to request faces
  iMesh_setAdjTable(cm->impl(), adj_table, 16, &err);
  CHECK_ERR("Failed to set adj table.");
	
  // create vertices
  verts = NULL; verts_alloc = 0;
  iMesh_createVtxArr(cm->impl(), 16, iBase_INTERLEAVED, vert_pos, 48,
                     &verts, &verts_alloc, &verts_size, &err);
  CHECK_ERR("Failed to create 16 vertices.");
	
  // pre-allocate ents space
  ents = (iBase_EntityHandle*) malloc(6*sizeof(iBase_EntityHandle));
  ents_alloc = 6;
	 
  // create hexes
  iBase_EntityHandle connect[sizeof(connect_idx)/sizeof(int)];
  int num_connect = sizeof(connect_idx)/sizeof(int);
  for (int i = 0; i < num_connect; i++)
    connect[i] = verts[connect_idx[i]];
	 
  int *status = NULL, status_alloc = 0, status_size;
  iMesh_createEntArr(impl,  iMesh_HEXAHEDRON, connect, 24,
                     &ents, &ents_alloc, &ents_size,
                     &status, &status_alloc, &status_size, &err);
  CHECK_ERR("Failed to create 3 hexes.");

  // create entity sets and add things to it
  iBase_EntitySetHandle esets[7];
  for (int i = 0; i < 7; i++) {
    iMesh_createEntSet(cm->impl(), false, esets+i, &err);
    CHECK_ERR("Failed to create entity set.");
  }
  
  iBase_TagHandle gtag, mtag, ntag, dtag;
  iMesh_getTagHandle(cm->impl(), "GEOM_DIMENSION", &gtag, &err, 14);
  if (iBase_TAG_NOT_FOUND == err) {
    iMesh_createTag(cm->impl(), "GEOM_DIMENSION", 1, iBase_INTEGER,
                    &gtag, &err, 14);
    CHECK_ERR("Failed to create tag handle for geom dimension.");
  }
  
  iMesh_getTagHandle(cm->impl(), "MATERIAL_SET", &mtag, &err, 12);
  if (iBase_TAG_NOT_FOUND == err) {
    iMesh_createTag(cm->impl(), "MATERIAL_SET", 1, iBase_INTEGER,
                    &mtag, &err, 12);
    CHECK_ERR("Failed to create tag handle for material set.");
  }
  
  iMesh_getTagHandle(cm->impl(), "NEUMANN_SET", &ntag, &err, 11);
  if (iBase_TAG_NOT_FOUND == err) {
    iMesh_createTag(cm->impl(), "NEUMANN_SET", 1, iBase_INTEGER,
                    &ntag, &err, 11);
    CHECK_ERR("Failed to create tag handle for neumann set.");
  }
  
  iMesh_getTagHandle(cm->impl(), "DIRICHLET_SET", &dtag, &err, 13);
  if (iBase_TAG_NOT_FOUND == err) {
    iMesh_createTag(cm->impl(), "DIRICHLET_SET", 1, iBase_INTEGER,
                    &dtag, &err, 13);
    CHECK_ERR("Failed to create tag handle for dirichlet set.");
  }
  
  // ES1: geom dim = 3, hexes 1-3, C
  iMesh_addEntArrToSet(cm->impl(), ents, 3, *esets, &err);
  CHECK_ERR("Failed to add ents to set 1.");
  int dum = 3;
  iMesh_setEntSetIntData(cm->impl(),esets[0], gtag, dum, &err);
  CHECK_ERR("Failed to set geom set tag on set 1.");
  
  // ES2: geom dim = 2, faces all, -
  iBase_EntityHandle *faces = NULL;
  int faces_alloc = 0, faces_size;
  iMesh_getEntities(cm->impl(), root_set, iBase_FACE, iMesh_ALL_TOPOLOGIES,
                    &faces, &faces_alloc, &faces_size, &err);
  CHECK_ERR("Failed to get faces.");
  iMesh_addEntArrToSet(cm->impl(), faces, faces_size, *(esets+1), &err);
  CHECK_ERR("Failed to add ents to set 2.");
  dum = 2;
  iMesh_setEntSetIntData(cm->impl(), esets[1], gtag, dum, &err);
  CHECK_ERR("Failed to set geom set tag on set 2.");
  free(faces);

  // ES3: - , skin nodes, E
  iMesh_addEntArrToSet(cm->impl(), verts, verts_size, *(esets+2), &err);
  CHECK_ERR("Failed to add ents to set 3.");
  cm->expand_sets().add_set(esets[2]);

  // ES4: mat set, set 1, E
  iMesh_addEntSet(cm->impl(), esets[0],*(esets+3), &err);
  CHECK_ERR("Failed to add set to set 4.");
  dum = 100;
  iMesh_setEntSetIntData(cm->impl(), esets[3], mtag, dum, &err);
  CHECK_ERR("Failed to set mat set tag on set 4.");

  // ES5: neu set, set 2, E
  iMesh_addEntSet(cm->impl(), esets[1],*(esets+4), &err);
  CHECK_ERR("Failed to add set to set 5.");
  dum = 101;
  iMesh_setEntSetIntData(cm->impl(), esets[4], ntag, dum, &err);
  CHECK_ERR("Failed to set neu set tag on set 5.");

  // ES6: dir set, set 3, E
  iMesh_addEntSet(cm->impl(), esets[2],*(esets+5), &err);
  CHECK_ERR("Failed to add set to set 6.");
  dum = 102;
  iMesh_setEntSetIntData(cm->impl(), esets[5], dtag, dum, &err);
  CHECK_ERR("Failed to set dir set tag on set 6.");

  // ES7: -, set 4, C
  iMesh_addEntSet(cm->impl(), esets[3],*(esets+6), &err);
  CHECK_ERR("Failed to add set to set 7.");
  cm->copy_sets().add_set(esets[6]);
	
  // check # entities
  if (check_num_ents(iBase_VERTEX, 16)) return 1;
  if (check_num_ents(iBase_REGION, 3)) return 1;

  return 0;
}
