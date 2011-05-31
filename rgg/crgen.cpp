/*********************************************
 June,10
 Reactor Assembly Mesh Assembler
 Argonne National Laboratory

 CCrgen class definition.
******************************************/
#define STRINGIFY_(X) #X
#define STRINGIFY(X) STRINGIFY_(X)
#define DIR STRINGIFY(SRCDIR) "/"
#define TEST_FILE_NAME "twoassm"
#define DEFAULT_TEST_FILE STRINGIFY(SRCDIR) "/twoassm"

#include "crgen.hpp"
#include "utils.hpp"
#include <string.h>
/* ==================================================================
   ======================= CCrgen class =============================
   ================================================================== */

int CCrgen::assign_gids_parallel(const int nrank, const int numprocs) {
#ifdef USE_MPI   
  // assign new global ids
  if (global_ids == true) {
    if(nrank==0)
      std::cout << "Assigning global ids in parallel" << std::endl;
    err = pc->assign_global_ids(0, 3, 1, false, true, false);
    ERRORR("Error assigning global ids.", err);
  }
  // // assign global ids for all entity sets
  // SimpleArray<iBase_EntitySetHandle> sets;
  // iBase_TagHandle gid_tag;
  // const char *tag_name = "GLOBAL_ID"; 
  // iMesh_getTagHandle(impl, tag_name, &gid_tag, &err, 9);
  // if (iBase_TAG_NOT_FOUND == err) {
  //   iMesh_createTag(impl, tag_name, 1, iBase_INTEGER, 
  //                   &gid_tag, &err, strlen(tag_name));
  //   ERRORR("Couldn't create global id tag", err);
  // }



  // iMesh_getEntSets(impl,root_set, 1, ARRAY_INOUT(sets), &err);
  // ERRORR("Failed to get contained sets.", err);

  // for(int id = 1; id <= sets.size(); id++){
  //   iMesh_setEntSetIntData(impl, sets[id-1], gid_tag, id, &err);
  //   ERRORR("Failed to set tags on sets.", err);
  // }
#endif
  return iBase_SUCCESS;
}

int CCrgen::close_parallel(const int nrank, const int numprocs)
// ---------------------------------------------------------------------------
// Function: dellocating
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
#ifdef USE_MPI
  // deallocate ... deallocate ... deallocate
  if (prob_type == "mesh") {
//    for (unsigned int i = 0; i < assys.size(); i++) {
//      delete cm[i];
//    }

    iMesh_dtor(impl, &err);
    ERRORR("Failed in call iMesh_dtor", err);

  }

  if (prob_type == "geometry") {
    for (unsigned int i = 0; i < 1; i++) {
      delete cg[i];
    }
    iGeom_dtor(geom, &err);
    ERRORR("Failed in call iGeom_dtor", err);
  }
#endif
  return 0;
}

int CCrgen::save_mesh_parallel(const int nrank, const int numprocs) {
  #ifdef USE_MPI
  // write file
  if (nrank == 0) {
    std::cout << "Saving mesh file in parallel. " << std::endl;
  }
#ifdef HAVE_MOAB
  int rval = mbImpl()->write_file(outfile.c_str() , 0,"PARALLEL=WRITE_PART");
  if(rval != moab::MB_SUCCESS) {
    std::cerr<<"Writing output file failed Code:";
    std::string foo = ""; mbImpl()->get_last_error(foo);
    std::cerr<<"File Error: "<<foo<<std::endl;
    return 1;
  }
  if (nrank == 0) {
    std::cout << "Done saving mesh file: " << outfile << std::endl;
  }
#endif
#endif
  return iBase_SUCCESS;

}

int CCrgen::merge_nodes_parallel(const int nrank, const int numprocs)
// -------------------------------------------------------------------------------------------
// Function: merge the nodes within a set tolerance in the model
// Input:    none
// Output:   none
// -------------------------------------------------------------------------------------------
{
#ifdef USE_MPI
  //Call the resolve parallel function
  moab::ParallelMergeMesh pm(pc, merge_tol);
  err = pm.merge();
  if (err != moab::MB_SUCCESS) {
    std::cerr << "Merge Failed" << std::endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
#endif
  return iBase_SUCCESS;
}

int CCrgen::distribute_mesh(const int nrank, const int numprocs)
// -------------------------------------------------------------------------------------------
// Function: merge the nodes within a set tolerance in the model
// Input:    none
// Output:   none
// -------------------------------------------------------------------------------------------
{
#ifdef USE_MPI 
  std::vector<int> rank_load;
  rank_load.resize(numprocs);
  int extra_procs = numprocs - files.size();
  if(numprocs > (int) files.size() && numprocs <= tot_assys){
    // again fill assm_meshfiles
    for(int p=0; p<nassys; p++)
      assm_meshfiles[p]=0;
   
     for(int p=0; p<tot_assys; p++){
       for(int q=0; q<nassys; q++){
	if(strcmp(core_alias[p].c_str(), assm_alias[q].c_str()) ==0) {
	 assm_meshfiles[q]+=1;
//	 assm_location[q].push_back(p);
         }
       }
    }
    //distribute
    for(int i=0; i< (int) files.size(); i++){
      rank_load[i] = i;
    }
    
    int temp = 0;
    int assm_load = - 1;
    int e= 0;
    // compute the rank, mf vector for extra procs
      while(e < extra_procs){
  	for(int i = 0; i < nassys; i++){
  	  if (assm_meshfiles[i] > temp){
  	    temp = assm_meshfiles[i];
  	    assm_load = i;
  	  }
  	  else if (assm_load == -1){
  	    std::cout << "No assemblies mesh files used in core" << std::endl;
  	    exit(0);
  	  }
  	}
  	assm_meshfiles[assm_load]-=1;
	int temp_rank = files.size()+ e;
  	rank_load[temp_rank] = assm_load;
	e++;
  	temp = 0;	  
      }
  }
  else{
    // bailout
    std::cout << "Invalid number of procs"<< std::endl;
    exit(0);
}
  std::vector<int> times_loaded(nassys);
  std::vector<std::vector<int> > meshfiles_rank (files.size());
  for(int i=0; i < (int) files.size(); i++){
    for(int j=0; j < (int) rank_load.size(); j++){
      if(rank_load[j]==i){
	meshfiles_rank[i].push_back(j);
	times_loaded[i]+=1;
      }
    }
  }

  position_core.resize(numprocs);
  for(int i=0; i < (int) files.size(); i++){
    int k = 0;
    for(int j=0; j < (int) assm_location[i].size(); j++){
      if (k >= (int) meshfiles_rank[i].size()){
	k = 0;
      }
      int p = meshfiles_rank[i][k];
      int q = assm_location[i][j];
      position_core[p].push_back(q);

      ++k;
    }
  }
  // for(int i =0; i< numprocs; i++){
  //   for(int j=0; j<position_core[i].size(); j++){
  //     std::cout << position_core[i][j] << " ";
  //   }
  //   std::cout << "\n" << std::endl;
  // }
#endif
  return 0;
}


int CCrgen::load_meshes_parallel(const int nrank, const int numprocs)
// ---------------------------------------------------------------------------
// Function: loads all the meshes and initializes copymesh object
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
#ifdef USE_MPI
  if(numprocs > (int) core_alias.size()){
    std::cout << "Error: #procs should be <= #assys in core" << std::endl;
    exit(0);
    //MPI::COMM_WORLD.Abort(0);
  }
  iBase_EntitySetHandle orig_set;
  // make a mesh instance, get root set, create ent set for loading meshfiles
  iMesh_newMesh("MOAB", &impl, &err, 4);
  ERRORR("Failed to create instance.", 1);

  iMesh_getRootSet(impl, &root_set, &err);
  ERRORR("Couldn't get the root set", err);
 
  int temp_index;
  int extra_procs = numprocs - files.size();
  std::vector<int> rank_load;
  for(int i = 0; i< (int) files.size(); i++){
    temp_index = numprocs*i + nrank;
    if(temp_index >= (int) files.size()){      
      if (nrank >= (int) files.size() && nrank <= tot_assys){
  	int temp = 0;
  	int assm_load = - 1;
	int p = extra_procs;
	// compute the rank, mf vector for extra procs
	while(p !=0){
	  for(int i = 0; i < nassys; i++){
	    if ((int) assm_meshfiles[i] > temp){
	      temp = assm_meshfiles[i];
	      assm_load = i;
	    }
	    else if (assm_load == -1){
	      std::cout << "No assemblies mesh files used in core" << std::endl;
	      exit(0);
	    }
	  }
	  assm_meshfiles[assm_load]-=1;
	  rank_load.push_back(assm_load);
	  --p;
	  temp = 0;	  
	}
	temp_index = nrank - files.size();
	iMesh_createEntSet(impl, 0, &orig_set, &err);
	ERRORR("Couldn't create file set.", err);
	// load this file
	iMesh_load(impl, orig_set, files[rank_load[temp_index]].c_str(), NULL, &err, strlen(files[rank_load[temp_index]].c_str()), 0);
	ERRORR("Couldn't read mesh file.", err);
	std::cout << "Loaded mesh file " << rank_load[temp_index] << " in processor: " << nrank << std::endl;
	assys.push_back(orig_set);
	assys_index.push_back(rank_load[temp_index]);
	break;
      }
    }
    else{
      iMesh_createEntSet(impl, 0, &orig_set, &err);
      ERRORR("Couldn't create file set.", err);
      // load this file
      iMesh_load(impl, orig_set, files[temp_index].c_str(), NULL, &err, strlen(files[temp_index].c_str()), 0);
      ERRORR("Couldn't read mesh file.", err);
      std::cout << "Loaded mesh file " << temp_index << " in processor: " << nrank << std::endl;
      assys.push_back(orig_set);
      assys_index.push_back(temp_index);
    }
  }
  
  // create cm instances for each mesh file
  cm = new CopyMesh*[assys.size()];
  for (unsigned int i = 0; i < assys.size(); i++) {
    cm[i] = new CopyMesh(impl);
  }
#endif
  return iBase_SUCCESS;
}


CCrgen::CCrgen()
// ---------------------------------------------------------------------------
// Function: default constructor
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
  err = 0;
  UNITCELL_DUCT = 0;
  ASSY_TYPES = 1;
  pack_type = 1;
  symm = 1;
  z_height = 1;
  z_divisions = 2;
  set_DIM = 3; // default is 3D
  PII = acos(-1.0);
  comment = "!";
  MAXCHARS = 300;
  extrude_flag = false;
  mem_tflag = false;
  global_ids = true;
  merge_tol = 1.0e-9;
  do_merge = 1;
  update_sets = 0;
  merge_tag = NULL;
  nst_flag = false;
  nsb_flag = false;
  nss_flag = false;
  nst_Id = 9997;
  nsb_Id = 9998;
  nss_Id = 9999;
  prob_type = "mesh";
}

CCrgen::~CCrgen()
// ---------------------------------------------------------------------------
// Function: destructor
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
}

int CCrgen::close()
// ---------------------------------------------------------------------------
// Function: dellocating 
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
  // deallocate ... deallocate ... deallocate
  if (prob_type == "mesh") {
    for (unsigned int i = 0; i < files.size(); i++) {
      delete cm[i];
    }

    iMesh_dtor(impl, &err);
    ERRORR("Failed in call iMesh_dtor", err);

  }
  if (prob_type == "geometry") {
    for (unsigned int i = 0; i < files.size(); i++) {
      delete cg[i];
    }
    iGeom_dtor(geom, &err);
    ERRORR("Failed in call iGeom_dtor", err);
  }
  return 0;
}

int CCrgen::prepareIO(int argc, char *argv[], int myID, int nprocs)
// -----------------------------------------------------------------------------------
// Function: Obtains file names and opens input/output files and then read/write them
// Input:    command line arguments
// Output:   none
// -----------------------------------------------------------------------------------
{
  bool bDone = false;
  do {
    if (2 == argc) {
      if (argv[1][0] == '-' && myID == 0) {
	if (argv[1][1] == 'h') {
	  std::cout
	    << "Usage: coregen [-t -m -h] <coregen input file>"
	    << std::endl;
	  std::cout
	    << "        -t print timing and memory usage info in each step"
	    << std::endl;
	  std::cout << "        -m create makefile only" << std::endl;
	  std::cout << "        -h print help" << std::endl;
	  std::cout
	    << "\nInstruction on writing coregen input file can be found at: "
	    << std::endl;
	  std::cout
	    << "        https://trac.mcs.anl.gov/projects/fathom/browser/MeshKit/trunk/rgg/README"
	    << std::endl;
	  exit(0);
	}
      }

      iname = argv[1];
      ifile = iname + ".inp";
      outfile = iname + ".h5m";
      mfile = iname + ".makefile";
    } else if (3 == argc) {
      int i = 1;// will loop through arguments, and process them
      for (i = 1; i < argc - 1; i++) {
	if (argv[i][0] == '-') {
	  switch (argv[i][1]) {
	  case 'm': {
	    if (myID == 0) {
	      std::cout << "Creating Makefile Only" << std::endl;
	    }
	    // only makefile creation specified
	    iname = argv[2];
	    ifile = iname + ".inp";
	    outfile = iname + ".h5m";
	    mfile = iname + ".makefile";
	    break;
	  }
	  case 't': {
	    mem_tflag = true;
	    iname = argv[2];
	    ifile = iname + ".inp";
	    outfile = iname + ".h5m";
	    mfile = iname + ".makefile";
	    break;
	  }
	  case 'h': {
	    if (myID == 0) {
	      std::cout
		<< "Usage: coregen [-t -m -h] <coregen input file>"
		<< std::endl;
	      std::cout
		<< "        -t print timing and memory usage info in each step"
		<< std::endl;
	      std::cout << "        -m create makefile only"
			<< std::endl;
	      std::cout << "        -h print help" << std::endl;
	      std::cout
		<< "\nInstruction on writing coregen input file can also be found at: "
		<< std::endl;
	      std::cout
		<< "        https://trac.mcs.anl.gov/projects/fathom/browser/MeshKit/trunk/rgg/README"
		<< std::endl;
	      exit(0);
	      break;
	    }
	  }
	  }
	}
      }
    } else { //default case
      if (myID == 0) {
	std::cerr << "Usage: " << argv[0]
		  << " <input file> WITHOUT EXTENSION" << std::endl;
	std::cout << "  No file specified.  Defaulting to: "
		  << DEFAULT_TEST_FILE << std::endl;
      }
      iname = DEFAULT_TEST_FILE;
      ifile = iname + ".inp";
      std::string temp = TEST_FILE_NAME;
      outfile = temp + ".h5m";
      mfile = temp + ".makefile";
    }

    // open the file
    file_input.open(ifile.c_str(), std::ios::in);
    if (!file_input) {
      if (myID == 0) {
	std::cout << "Unable to open file" << std::endl;
	std::cout << "Usage: coregen [-t -m -h] <coregen input file>"
		  << std::endl;
	std::cout
	  << "        -t print timing and memory usage info in each step"
	  << std::endl;
	std::cout << "        -m create makefile only" << std::endl;
	std::cout << "        -h print help" << std::endl;
	std::cout
	  << "\nInstruction on writing coregen input file can be found at: "
	  << std::endl;
	std::cout
	  << "        https://trac.mcs.anl.gov/projects/fathom/browser/MeshKit/trunk/rgg/README"
	  << std::endl;
      }
      file_input.clear();
      exit(1);
    } else
      bDone = true; // file opened successfully
  } while (!bDone);

  // open Makefile-rgg
  do {
    make_file.open(mfile.c_str(), std::ios::out);
    if (!make_file) {
      if (myID == 0) {
	std::cout << "Unable to open makefile for writing" << std::endl;
      }
      make_file.clear();
    } else
      bDone = true; // file opened successfully
  } while (!bDone);

  if (myID == 0) {
    std::cout << "\nEntered input file name: " << ifile << std::endl;
  }

  // now call the functions to read and write
  err = read_inputs_phase1();
  ERRORR("Failed to read inputs in phase1.", 1);

  err = read_inputs_phase2();
  ERRORR("Failed to read inputs in phase2.", 1);

  if (myID == 0) {
    err = write_makefile();
    ERRORR("Failed to write a makefile.", 1);
  }
  return 0;
}

int CCrgen::load_meshes()
// ---------------------------------------------------------------------------
// Function: loads all the meshes and initializes copymesh and merge mesh objects
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
  // make a mesh instance
  iMesh_newMesh("MOAB", &impl, &err, 4);
  ERRORR("Failed to create instance.", 1);

  iMesh_getRootSet(impl, &root_set, &err);
  ERRORR("Couldn't get the root set", err);

  // create cm instances for each mesh file
  cm = new CopyMesh*[files.size()];
  for (unsigned int i = 0; i < files.size(); i++) {
    cm[i] = new CopyMesh(impl);
  }

  iBase_EntitySetHandle orig_set;

  // loop reading all mesh files
  for (unsigned int i = 0; i < files.size(); i++) {
    iMesh_createEntSet(impl, 0, &orig_set, &err);
    ERRORR("Couldn't create file set.", err);

    iMesh_load(impl, orig_set, files[i].c_str(), NULL, &err, strlen(
								    files[i].c_str()), 0);
    ERRORR("Couldn't read mesh file.", err);

    assys.push_back(orig_set);
    assys_index.push_back(i);
  }
  std::cout << "Loaded mesh files." << std::endl;

  return iBase_SUCCESS;
}

int CCrgen::load_geometries()
// ---------------------------------------------------------------------------
// Function: loads all the meshes and initializes copymesh and merge mesh objects
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
  std::cout << "\n--Loading geometry files." << std::endl;
  // make a mesh instance
  iGeom_newGeom("GEOM", &geom, &err, 4);
  ERRORR("Failed to create instance.", 1);

  iGeom_getRootSet(geom, &root_set, &err);
  ERRORR("Couldn't get the root set", err);

  // create cm instances for each mesh file
  cg = new CopyGeom*[files.size()];
  for (unsigned int i = 0; i < files.size(); i++) {
    cg[i] = new CopyGeom(geom);
  }

  iBase_EntitySetHandle orig_set, temp_set, temp_set1;

  iGeom_createEntSet(geom, 0, &temp_set1, &err);
  ERRORR( "Problem creating entity set.",err );

  // loop reading all geom files
  for (unsigned int i = 0; i < files.size(); i++) {
    iGeom_createEntSet(geom, 0, &orig_set, &err);
    ERRORR( "Problem creating entity set.", err);

    iGeom_createEntSet(geom, 0, &temp_set, &err);
    ERRORR( "Problem creating entity set.",err );

    iGeom_load(geom, files[i].c_str(), NULL, &err,
	       strlen(files[i].c_str()), 0);
    ERRORR("Couldn't read mesh file.", err);

    iBase_EntityHandle *entities = NULL;
    int entities_ehsize = 0, entities_ehallocated = 0;

    iGeom_getEntities(geom, root_set, iBase_REGION, &entities,
		      &entities_ehsize, &entities_ehallocated, &err);
    ERRORR( "Problem getting entities." , err);

    // add the entity
    for (int j = 0; j < entities_ehsize; j++) {
      iGeom_addEntToSet(geom, entities[j], temp_set, &err);
      ERRORR( "Problem adding to set.", err );
    }

    iGeom_subtract(geom, temp_set, temp_set1, &orig_set, &err);
    ERRORR( "Unable to subtract entity sets.", err );

    assys.push_back(orig_set);
    // store this set for subtraction with next entity set
    temp_set1 = temp_set;
  }
  std::cout << "\n--Loaded geometry files.\n" << std::endl;

  return iBase_SUCCESS;
}

int CCrgen::read_inputs_phase1() {
  // ---------------------------------------------------------------------------
  // Function: Reads the dimension and symmetry of the problem
  // Input:    none
  // Output:   none
  // ---------------------------------------------------------------------------
  CParser parse;
  for (;;) {
    if (!parse.ReadNextLine(file_input, linenumber, input_string, MAXCHARS,
			    comment))
      ERRORR("Reading input file failed",1);

    if (input_string.substr(0, 11) == "problemtype") {
      std::istringstream formatString(input_string);
      formatString >> card >> prob_type;
      if ((strcmp(prob_type.c_str(), "geometry") == 0)) {
	prob_type = "geometry";
      }
    }
    if (input_string.substr(0, 8) == "geometry") {
      std::istringstream formatString(input_string);
      formatString >> card >> geometry;
      if ((strcmp(geometry.c_str(), "surface") == 0)) {
	set_DIM = 2;
      }
    }

    // geom engine
    if (input_string.substr(0, 10) == "geomengine") {
      std::istringstream formatString(input_string);
      formatString >> card >> geom_engine;
    }

    // symmetry
    if (input_string.substr(0, 8) == "symmetry") {
      std::istringstream formatString(input_string);
      formatString >> card >> symm;
    }

    // merge tolerance
    if (input_string.substr(0, 14) == "mergetolerance") {
      std::istringstream formatString(input_string);
      formatString >> card >> merge_tol;
    }

    // neumannset card
    if (input_string.substr(0, 10) == "neumannset") {
      std::istringstream formatString(input_string);
      std::string nsLoc = "";
      int nsId = 0;
      formatString >> card >> nsLoc >> nsId;
      if ((strcmp(nsLoc.c_str(), "top") == 0)) {
	nst_flag = true;
	nst_Id = nsId;
      } else if ((strcmp(nsLoc.c_str(), "bot") == 0)) {
	nsb_flag = true;
	nsb_Id = nsId;
      } else if ((strcmp(nsLoc.c_str(), "side") == 0)) {
	nss_flag = true;
	nss_Id = nsId;
      } else {
	std::cout << "Invalid Neumann set specification" << std::endl;
      }
    }

    // breaking condition
    if (input_string.substr(0, 3) == "end") {
      std::istringstream formatstring(input_string);
      break;
    }
  }
  return iBase_SUCCESS;
}

int CCrgen::read_inputs_phase2()
// ---------------------------------------------------------------------------
// Function: read all the inputs 
// Input:    command line arguments
// Output:   none
// ---------------------------------------------------------------------------
{
  //Rewind the input file
  file_input.clear(std::ios_base::goodbit);
  file_input.seekg(0L, std::ios::beg);
  linenumber = 0;

  CParser parse;
  for (;;) {
    if (!parse.ReadNextLine(file_input, linenumber, input_string, MAXCHARS,
			    comment))
      ERRORR("Reading input file failed",1);

    if (input_string.substr(0, 12) == "geometrytype") {
      std::istringstream formatString(input_string);
      formatString >> card >> geom_type;
      if (geom_type == "hexvertex" && symm == 6) {

	// reading pitch info
	if (!parse.ReadNextLine(file_input, linenumber, input_string,
				MAXCHARS, comment))
	  ERRORR("Reading input file failed",1);
	if (input_string.substr(0, 10) == "assemblies") {
	  std::istringstream formatString(input_string);
	  formatString >> card >> nassys >> pitch;
	}

	// reading file and alias names
	for (int i = 1; i <= nassys; i++) {
	  if (!parse.ReadNextLine(file_input, linenumber,
				  input_string, MAXCHARS, comment))
	    ERRORR("Reading input file failed",1);
	  std::istringstream formatString(input_string);
	  formatString >> meshfile >> mf_alias;
	  meshfile = DIR + meshfile;
	  files.push_back(meshfile);
	  assm_alias.push_back(mf_alias);
	}
	// reading lattice
	if (!parse.ReadNextLine(file_input, linenumber, input_string,
				MAXCHARS, comment))
	  ERRORR("Reading input file failed",1);
	if (input_string.substr(0, 7) == "lattice") {
	  std::istringstream formatString(input_string);
	  formatString >> card >> nrings;
	  if (nrings % 2 == 0)
	    tot_assys = (nrings * (nrings)) / 2;
	  else
	    tot_assys = ((nrings * (nrings - 1)) / 2)
	      + (nrings + 1) / 2;
	}

	// now reading the arrangement
	if (!parse.ReadNextLine(file_input, linenumber, input_string,
				MAXCHARS, comment))
	  ERRORR("Reading input file failed",1);
	std::istringstream formatString(input_string);
	for (int i = 1; i <= tot_assys; i++) {
	  formatString >> temp_alias;
	  core_alias.push_back(temp_alias);
	}
      }

      else if (geom_type == "rectangular" && symm == 1) {

	// reading pitch info
	if (!parse.ReadNextLine(file_input, linenumber, input_string,
				MAXCHARS, comment))
	  ERRORR("Reading input file failed",1);
	if (input_string.substr(0, 10) == "assemblies") {
	  std::istringstream formatString(input_string);
	  formatString >> card >> nassys >> pitchx >> pitchy;
	}

	// reading file and alias names
	for (int i = 1; i <= nassys; i++) {
	  if (!parse.ReadNextLine(file_input, linenumber,
				  input_string, MAXCHARS, comment))
	    ERRORR("Reading input file failed",1);
	  std::istringstream formatString(input_string);
	  formatString >> meshfile >> mf_alias;
	  meshfile = DIR + meshfile;
	  files.push_back(meshfile);
	  assm_alias.push_back(mf_alias);
	}

	// reading lattice
	if (!parse.ReadNextLine(file_input, linenumber, input_string,
				MAXCHARS, comment))
	  ERRORR("Reading input file failed",1);
	if (input_string.substr(0, 7) == "lattice") {
	  std::istringstream formatString(input_string);
	  formatString >> card >> nringsx >> nringsy;
	  tot_assys = nringsx * nringsy;
	}

	// now reading the arrangement
	if (!parse.ReadNextLine(file_input, linenumber, input_string,
				MAXCHARS, comment))
	  ERRORR("Reading input file failed",1);
	std::istringstream formatString(input_string);
	for (int i = 1; i <= tot_assys; i++) {
	  formatString >> temp_alias;
	  core_alias.push_back(temp_alias);
	}
      }

      else if (geom_type == "hexflat" && symm == 6) {
	// reading pitch info
	if (!parse.ReadNextLine(file_input, linenumber, input_string,
				MAXCHARS, comment))
	  ERRORR("Reading input file failed",1);
	if (input_string.substr(0, 10) == "assemblies") {
	  std::istringstream formatString(input_string);
	  formatString >> card >> nassys >> pitch;
	}

	// reading file and alias names
	for (int i = 1; i <= nassys; i++) {
	  if (!parse.ReadNextLine(file_input, linenumber,
				  input_string, MAXCHARS, comment))
	    ERRORR("Reading input file failed",1);
	  std::istringstream formatString(input_string);
	  formatString >> meshfile >> mf_alias;
	  meshfile = DIR + meshfile;
	  files.push_back(meshfile);
	  assm_alias.push_back(mf_alias);
	}

	// reading lattice
	if (!parse.ReadNextLine(file_input, linenumber, input_string,
				MAXCHARS, comment))
	  ERRORR("Reading input file failed",1);
	if (input_string.substr(0, 7) == "lattice") {
	  std::istringstream formatString(input_string);
	  formatString >> card >> nrings;
	  tot_assys = (nrings * (nrings + 1)) / 2;
	}

	// now reading the arrangement
	if (!parse.ReadNextLine(file_input, linenumber, input_string,
				MAXCHARS, comment))
	  ERRORR("Reading input file failed",1);
	std::istringstream formatString(input_string);
	for (int i = 1; i <= tot_assys; i++) {
	  formatString >> temp_alias;
	  core_alias.push_back(temp_alias);
	}
      } else if (geom_type == "hexflat" && symm == 1) {
	// reading pitch info
	if (!parse.ReadNextLine(file_input, linenumber, input_string,
				MAXCHARS, comment))
	  ERRORR("Reading input file failed",1);
	if (input_string.substr(0, 10) == "assemblies") {
	  std::istringstream formatString(input_string);
	  formatString >> card >> nassys >> pitch;
	}

	// reading file and alias names
	for (int i = 1; i <= nassys; i++) {
	  if (!parse.ReadNextLine(file_input, linenumber,
				  input_string, MAXCHARS, comment))
	    ERRORR("Reading input file failed",1);
	  std::istringstream formatString(input_string);
	  formatString >> meshfile >> mf_alias;
	  meshfile = DIR + meshfile;
	  files.push_back(meshfile);
	  assm_alias.push_back(mf_alias);
	}

	// reading lattice
	if (!parse.ReadNextLine(file_input, linenumber, input_string,
				MAXCHARS, comment))
	  ERRORR("Reading input file failed",1);
	if (input_string.substr(0, 7) == "lattice") {
	  std::istringstream formatString(input_string);
	  formatString >> card >> nrings;
	  tot_assys = 3 * (nrings * (nrings - 1)) + 1;
	}

	// now reading the arrangement
	if (!parse.ReadNextLine(file_input, linenumber, input_string,
				MAXCHARS, comment))
	  ERRORR("Reading input file failed",1);
	std::istringstream formatString(input_string);
	for (int i = 1; i <= tot_assys; i++) {
	  formatString >> temp_alias;
	  core_alias.push_back(temp_alias);
	}
      } else if (geom_type == "hexflat" && symm == 12) {

	// reading pitch info
	if (!parse.ReadNextLine(file_input, linenumber, input_string,
				MAXCHARS, comment))
	  ERRORR("Reading input file failed",1);
	if (input_string.substr(0, 10) == "assemblies") {
	  std::istringstream formatString(input_string);
	  formatString >> card >> nassys >> pitch;
	}

	// reading file and alias names
	for (int i = 1; i <= nassys; i++) {
	  if (!parse.ReadNextLine(file_input, linenumber,
				  input_string, MAXCHARS, comment))
	    ERRORR("Reading input file failed",1);
	  std::istringstream formatString(input_string);
	  formatString >> meshfile >> mf_alias;
	  meshfile = DIR + meshfile;
	  files.push_back(meshfile);
	  assm_alias.push_back(mf_alias);
	}

	// reading lattice
	if (!parse.ReadNextLine(file_input, linenumber, input_string,
				MAXCHARS, comment))
	  ERRORR("Reading input file failed",1);
	if (input_string.substr(0, 7) == "lattice") {
	  std::istringstream formatString(input_string);
	  formatString >> card >> nrings;
	  if (nrings % 2 == 0)
	    tot_assys = (nrings * (nrings + 2)) / 4;
	  else
	    tot_assys = ((nrings + 1) * (nrings + 1)) / 4;
	}

	// now reading the arrangement
	if (!parse.ReadNextLine(file_input, linenumber, input_string,
				MAXCHARS, comment))
	  ERRORR("Reading input file failed",1);
	std::istringstream formatString(input_string);
	for (int i = 1; i <= tot_assys; i++) {
	  formatString >> temp_alias;
	  core_alias.push_back(temp_alias);
	}
      }

      else {
	ERRORR("Invalid geometry type",1);
      }
    }
    // background mesh
    if (input_string.substr(0, 10) == "background") {
      std::istringstream formatString(input_string);
      formatString >> card >> back_meshfile;
      back_meshfile = DIR + back_meshfile;
      files.push_back(back_meshfile);
      back_mesh = true;
    }
    // z-height and z-divisions
    if (input_string.substr(0, 7) == "extrude") {
      extrude_flag = true;
      std::istringstream formatString(input_string);
      formatString >> card >> z_height >> z_divisions;
    }

    // OutputFileName
    if (input_string.substr(0, 14) == "outputfilename") {
      std::istringstream formatString(input_string);
      formatString >> card >> outfile;
    }

    // breaking condition
    if (input_string.substr(0, 3) == "end") {
      std::istringstream formatstring(input_string);
      break;
    }
  }
  // set some variables
  assm_meshfiles.resize(nassys);

  int flag = 0;
  assm_location.resize(nassys);
  for(int i = 0; i < tot_assys; i++){
    for (int j = 0; j < nassys; j++){
      if (strcmp(core_alias[i].c_str(), assm_alias[j].c_str()) == 0) {
	assm_meshfiles[j]+=1;
	assm_location[j].push_back(i);
	flag = 1;
	break;
      }
    }
  }
  // for(int i = 0; i < nassys; i++){
  //   for (int j = 0; j < (int) assm_location[i].size(); j++){
  //     //      std::cout << assm_location[i][j]<<std::endl;
  //   }
  // }

  return iBase_SUCCESS;
}

int CCrgen::find_assm(const int i, int &assm_index)
// ---------------------------------------------------------------------------
// Function: find the assembly index (0 to n) for n assemblies for core alias i
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
  int flag = 0;
  for (int j = 0; j < nassys; j++)
    if (strcmp(core_alias[i].c_str(), assm_alias[j].c_str()) == 0) {
      assm_index = j;
      flag = 1;
      break;
    }
  if (flag == 0)//nothing found return -1 or no assembly
    assm_index = -1;
  return iBase_SUCCESS;
}

int CCrgen::banner()
// ---------------------------------------------------------------------------
// Function: display the program banner
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
  std::cout << '\n';
  std::cout
    << "\t\t--------------------------------------------------------------------"
    << '\n';
  std::cout
    << "\t\tProgram to Assemble Nuclear Reactor Assembly Meshes and Form a Core     "
    << '\n';
  std::cout << "\t\t\t\t\tArgonne National Laboratory" << '\n';
  std::cout << "\t\t\t\t\t        2009-2010         " << '\n';
  std::cout
    << "\t\t--------------------------------------------------------------------"
    << '\n';
  return iBase_SUCCESS;
}

int CCrgen::write_makefile()
// ---------------------------------------------------------------------------
// Function: write the makefile based on inputs read from input file
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
  std::string name;
  std::vector<std::string> f_no_ext, f_sat, f_inp, f_jou, f_injou;
  make_file << "##" << std::endl;
  make_file
    << "## This makefile is automatically generated by coregen program"
    << std::endl;
  make_file << "##" << std::endl;
  make_file << "## Check your coregen, assygen and cubit location"
	    << std::endl;
  make_file << "##" << std::endl;
  make_file << "\nCUBIT = cubit\n" << std::endl;
  make_file << "COREGEN = ../../coregen\n" << std::endl;
  make_file << "ASSYGEN = ../../assygen\n" << std::endl;

  // remove the ./ if run from the current working directory
  make_file << "MESH_FILES = ";
  std::string filename;
  for (unsigned int i = 0; i < files.size(); i++) {
    if (files[i][0] == '.' && files[i][1] == '/') {
      filename = files[i].substr(2, files[i].length());
    } else if (files[i][0] != '.' || files[i][1] != '/') {
      int loc1 = files[i].find_last_of(".");
      int loc2 = files[i].find_last_of("/");
      filename = files[i].substr(loc2 + 1, loc1);
    } else {
      filename = files[i];
    }
    mk_files.push_back(filename);
    make_file << mk_files[i] << "  ";
  }

  // get file names without extension
  for (unsigned int i = 0; i < mk_files.size(); i++) {
    int loc = mk_files[i].find_first_of(".");
    f_no_ext.push_back(mk_files[i].substr(0, loc));
  }

  make_file << "\n\nGEOM_FILES = ";
  for (unsigned int i = 0; i < mk_files.size(); i++) {
    if (geom_engine == "occ")
      name = f_no_ext[i] + ".stp";
    else
      name = f_no_ext[i] + ".sat";
    f_sat.push_back(name);
    make_file << name << "  ";
    name = "";
  }

  make_file << "\n\nJOU_FILES = ";
  for (unsigned int i = 0; i < mk_files.size(); i++) {
    name = f_no_ext[i] + ".jou";
    f_jou.push_back(name);
    make_file << name << "  ";
    name = "";
  }

  make_file << "\n\nINJOU_FILES = ";
  for (unsigned int i = 0; i < mk_files.size(); i++) {
    name = f_no_ext[i] + ".template.jou";
    f_injou.push_back(name);
    make_file << name << "  ";
    name = "";
  }

  make_file << "\n\nASSYGEN_FILES = ";
  for (unsigned int i = 0; i < mk_files.size(); i++) {
    name = f_no_ext[i] + ".inp";
    f_inp.push_back(name);
    make_file << name << "  ";
    name = "";
  }

  make_file << "\n\n" << outfile << " : ${MESH_FILES} " << ifile << std::endl;
  make_file << "\t" << "${COREGEN} " << iname << std::endl;
  for (unsigned int i = 0; i < mk_files.size(); i++) {
    make_file << mk_files[i] << " : " << f_sat[i] << "  " << f_jou[i]
	      << "  " << f_injou[i] << std::endl;
    make_file << "\t" << "${CUBIT} -batch " << f_jou[i] << "\n"
	      << std::endl;

    make_file << f_sat[i] << " " << f_jou[i] << " " << f_injou[i] << " : "
	      << f_inp[i] << std::endl;
    make_file << "\t" << "${ASSYGEN} " << f_no_ext[i] << "\n" << std::endl;
  }

  make_file.close();
  std::cout << "Created makefile: " << mfile << std::endl;
  return 0;
}

int CCrgen::move_verts(iBase_EntitySetHandle set, const double *dx)
// ---------------------------------------------------------------------------
// Function: Change the coordinates for moving the assembly to its first loc.
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{

  int verts_ents_alloc = 0, verts_ents_size = 0;
  iBase_EntityHandle *verts_ents = NULL;

  iMesh_getEntities(impl, set, iBase_VERTEX, iMesh_ALL_TOPOLOGIES,
		    &verts_ents, &verts_ents_alloc, &verts_ents_size, &err);
  ERRORR("Failed to get any entities from original set.", iBase_FAILURE);

  double *coords = 0;
  int coords_alloc = 0, coords_size = 0;

  iMesh_getVtxArrCoords(impl, verts_ents, verts_ents_size, iBase_INTERLEAVED,
			&coords, &coords_alloc, &coords_size, &err);
  ERRORR("Failed to get vtx coords from set.", iBase_FAILURE);

  for (int i = 0; i < verts_ents_size; i++) {
    coords[3 * i] += dx[0];
    coords[3 * i + 1] += dx[1];
    coords[3 * i + 2] += dx[2];
  }

  iMesh_setVtxArrCoords(impl, verts_ents, verts_ents_size, iBase_INTERLEAVED,
			coords, coords_size, &err);
  ERRORR("Failed to set vtx coords.", iBase_FAILURE);

  return iBase_SUCCESS;
}

int CCrgen::move_geoms(iBase_EntitySetHandle set, const double *dx)
// ---------------------------------------------------------------------------
// Function: Change the coordinates for moving the assembly to its first loc.
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
  int entities_ehsize = 0, entities_ehallocated = 0;
  iBase_EntityHandle *entities = NULL;

  iGeom_getEntities(geom, set, iBase_ALL_TYPES, &entities, &entities_ehsize,
		    &entities_ehallocated, &err);
  ERRORR("Failed to get entities from set.", iBase_FAILURE);

  for (int i = 0; i < entities_ehsize; i++) {
    iGeom_moveEnt(geom, entities[i], dx[0], dx[1], dx[2], &err);
    ERRORR("Failed to move geometries.", iBase_FAILURE);
  }
  return iBase_SUCCESS;
}

int CCrgen::merge_nodes()
// -------------------------------------------------------------------------------------------
// Function: merge the nodes within a set tolerance in the model
// Input:    none
// Output:   none
// -------------------------------------------------------------------------------------------
{

  //get entities for merging
  iBase_EntityHandle *ents = NULL;
  int ents_alloc = 0, ents_size;
  mm = new MergeMesh(impl);
  if (set_DIM == 2) { // if surface geometry specified
    iMesh_getEntities(impl, root_set, iBase_FACE, iMesh_ALL_TOPOLOGIES,
		      &ents, &ents_alloc, &ents_size, &err);
  } else {
    iMesh_getEntities(impl, root_set, iBase_REGION, iMesh_ALL_TOPOLOGIES,
		      &ents, &ents_alloc, &ents_size, &err);
  }
  ERRORR("Trouble getting entities for merge.", err);

  int num1, num2;

  iMesh_getNumOfType(impl, root_set, iBase_VERTEX, &num1, &err);
  ERRORR("Trouble getting number of entities before merge.", err);

  // merge now
  mm->merge_entities(ents, ents_size, merge_tol, do_merge, update_sets,
		     merge_tag);

  iMesh_getNumOfType(impl, root_set, iBase_VERTEX, &num2, &err);
  ERRORR("Trouble getting number of entities after merge.", err);

  std::cout << "Merged " << num1 - num2 << " vertices." << std::endl;
  delete mm;
  return iBase_SUCCESS;
}

int CCrgen::assign_gids() {
  // assign new global ids
  if (global_ids == true){
    std::cout << "Assigning global ids." << std::endl;
    mu = new MKUtils(impl);
    err = mu->assign_global_ids(root_set, 3, 1, true, false,
  				"GLOBAL_ID");
    ERRORR("Error assigning global ids.", err);
  }
  delete mu;
  return iBase_SUCCESS;


  // // assign new global ids
  // if (global_ids == true) {
  //   std::cout << "Assigning global ids." << std::endl;
  //   err = pc->assign_global_ids(0, 3, 1, false, false, false);
  //   ERRORR("Error assigning global ids.", err);
  // }

  // // assign global ids for all entity sets
  // SimpleArray<iBase_EntitySetHandle> sets;
  // iBase_TagHandle gid_tag;
  // const char *tag_name = "GLOBAL_ID"; 
  // iMesh_getTagHandle(impl, tag_name, &gid_tag, &err, 9);
  // if (iBase_TAG_NOT_FOUND == err) {
  //   iMesh_createTag(impl, tag_name, 1, iBase_INTEGER, 
  //                   &gid_tag, &err, strlen(tag_name));
  //   ERRORR("Couldn't create global id tag", err);
  // }

  // iMesh_getEntSets(impl,root_set, 1, ARRAY_INOUT(sets), &err);
  // ERRORR("Failed to get contained sets.", err);

  // for(int id = 1; id <= sets.size(); id++){
  //   iMesh_setEntSetIntData(impl, sets[id-1], gid_tag, id, &err);
  //   ERRORR("Failed to set tags on sets.", err);
  // }
  // return iBase_SUCCESS;
}

int CCrgen::save_mesh() {
  // export
  std::cout << "Saving mesh file." << std::endl;
  iMesh_save(impl, root_set, outfile.c_str(), NULL, &err, strlen(
								 outfile.c_str()), 0);
  ERRORR("Trouble writing output mesh.", err);
  std::cout << "Saved mesh file: " << outfile.c_str() << std::endl;

  return iBase_SUCCESS;
}

int CCrgen::save_geometry() {
  double dTol = 1e-3;

  // getting all entities for merge and imprint
  SimpleArray<iBase_EntityHandle> entities_merge, entities_imprint;
  iGeom_getEntities(geom, root_set, iBase_REGION,
		    ARRAY_INOUT(entities_merge), &err );
  ERRORR("Trouble writing output geometry.", err);

  // merge and imprint before save
  std::cout << "Merging.." << std::endl;
  iGeom_mergeEnts(geom, ARRAY_IN(entities_merge), dTol, &err);
  ERRORR("Trouble writing output geometry.", err);

  iGeom_getEntities( geom, root_set, iBase_REGION, ARRAY_INOUT(entities_imprint),&err );
  ERRORR("Trouble writing output geometry.", err);
 
  // std::cout << "Imprinting.." << std::endl;
  // iGeom_imprintEnts(geom, ARRAY_IN(entities_imprint),&err); 
  // ERRORR("Trouble writing output geometry.", err);
  // export
  std::cout << "Saving geometry file: " <<  outfile << std::endl;

  iGeom_save(geom, outfile.c_str(), NULL, &err,
             strlen(outfile.c_str()), 0);
  ERRORR("Trouble writing output geometry.", err);
  std::cout << "Saved geometry file: "<< outfile.c_str() <<std::endl;

  return iBase_SUCCESS;
}

int CCrgen::extrude() {
  // extrude if this is a surface mesh
  if (set_DIM == 2 && extrude_flag == true) { // if surface geometry and extrude
    std::cout << "Extruding surface mesh." << std::endl;

    ExtrudeMesh *ext = new ExtrudeMesh(impl);

    //get entities for extrusion
    iBase_EntityHandle *ents = NULL;
    int ents_alloc = 0, ents_size;
    iMesh_getEntities(impl, root_set, iBase_FACE, iMesh_ALL_TOPOLOGIES,
		      &ents, &ents_alloc, &ents_size, &err);
    ERRORR("Trouble getting face mesh.", err);

    // add entities for extrusion to a set
    iBase_EntitySetHandle set;
    iMesh_createEntSet(impl, false, &set, &err);
    ERRORR("Trouble getting face mesh.", err);

    iMesh_addEntArrToSet(impl, ents, ents_size, set, &err);
    ERRORR("Trouble getting face mesh.", err);

    // This tag needs to be set to the newly created extrude sets
    const char *tag_g1 = "GEOM_DIMENSION";
    iBase_TagHandle gtag;
    iMesh_getTagHandle(impl, tag_g1, &gtag, &err, 14);
    ERRORR("Trouble getting geom dimension set.", err);

    // This tag needs to be set to the newly created extrude sets
    const char *tag_m1 = "MATERIAL_SET";
    iBase_TagHandle mtag;
    iMesh_getTagHandle(impl, tag_m1, &mtag, &err, 12);
    ERRORR("Trouble getting material set.", err);

    // This tag needs to be set to the newly created extrude sets
    const char *tag_n1 = "NEUMANN_SET";
    iBase_TagHandle ntag;
    iMesh_getTagHandle(impl, tag_n1, &ntag, &err, 11);
    ERRORR("Trouble getting neumann set.", err);

    // add only material set tag as extrude set
    ext->extrude_sets().add_tag(mtag, NULL);

    // add 2d neumann set tag as copy tag to create-
    // -new neumann set for destination surfaces formed by extrusion
    ext->copy_sets().add_tag(ntag, NULL);

    double v[] = { 0, 0, z_height };
    int steps = z_divisions;

    // now extrude
    ext->extrude(set, extrude::Translate(v, steps));

    // update the copy sets, this creates the ent set for destination 2D quads.
    ext->copy_sets().update_tagged_sets();

    iMesh_destroyEntSet(impl, set, &err);
    ERRORR("Error in destroying ent set of faces after extrusion is done.", err);

    // Do things to get the metadata right after extrusion
    // Step 1: get all the material sets, remove old one's and add GD=3 on the new one's.

    SimpleArray<iBase_EntitySetHandle> msets;
    iMesh_getEntSetsByTagsRec(impl, root_set, &mtag, NULL, 1, 0,
			      ARRAY_INOUT(msets), &err);
    ERRORR("Trouble getting entity set.", err);

    for (int i = 0; i < msets.size(); i++) {
      int num =0;
      SimpleArray<iBase_EntitySetHandle> in_msets;

      iMesh_getNumOfType(impl, msets[i], iBase_REGION, &num, &err);
      ERRORR("Trouble getting num entities.", err);
      if(num ==0){
	iMesh_destroyEntSet(impl, msets[i], &err); 
	ERRORR("Trouble destroying set.", err);
      }
      else{
	const int gd = 3;
	iMesh_setEntSetIntData(impl, msets[i], gtag, gd, &err);
	ERRORR("Trouble setting tag data.", err);
      }
    }

    // Step 2: get all max. value of neumann sets, then, for newly created NS set a new value and GD =2 tag.

    iBase_EntityHandle *ents1 = NULL;
    int ents_alloc1 = 0, ents_size1;

    SimpleArray<iBase_EntitySetHandle> nsets;
    iMesh_getEntSetsByTagsRec(impl, root_set, &ntag, NULL,
			      1, 0, ARRAY_INOUT(nsets), &err);
    ERRORR("Trouble getting entity set.", err);

    int max_nset_value = 0;
    for (int i = 0; i < nsets.size(); i++) {
      int nvalue;
      iMesh_getEntSetIntData(impl, nsets[i], ntag, &nvalue, &err);
      ERRORR("Trouble getting entity set.", err);
      if (nvalue > max_nset_value)
	max_nset_value = nvalue;
    }
    
    for (int i = 0; i < nsets.size(); i++) {

      iMesh_getEntities(impl, nsets[i],
			iBase_FACE, iMesh_ALL_TOPOLOGIES,
			&ents1, &ents_alloc1, &ents_size1, &err);
      ERRORR("Trouble getting face mesh.", err);

      if(ents_size1 > 0) {
	// set GEOM_DIMENSION tag = 2 and renumber the neumann set
	const int gd = 2;
	const int nvalue = max_nset_value + i;

	iMesh_setEntSetIntData(impl, nsets[i], ntag, nvalue, &err);
	ERRORR("Trouble getting entity set.", err);

	iMesh_setEntSetIntData(impl,nsets[i], gtag, gd, &err);
	ERRORR("Trouble getting entity set.", err);
      }
    }
    delete ext;    
  }
  return iBase_SUCCESS;
}

int CCrgen::create_neumannset() {
  if (nss_flag == true || nsb_flag == true || nst_flag == true) {
    std::cout << "Creating NeumannSet." << std::endl;

    if (extrude_flag == true)
      set_DIM = 3;

#ifdef HAVE_MOAB  
    int err = 0, z_flag = 0, i, ents_alloc = 0, ents_size;
    double z1 = 0.0, z2 = 0.0;
    iBase_TagHandle ntag1, gtag1;
    iBase_EntityHandle *ents = NULL;
    iBase_EntitySetHandle set = NULL, set_z1 = NULL, set_z2 = NULL;

    //get entities for skinner
    if(set_DIM ==2) { // if surface geometry specified
      iMesh_getEntities(impl, root_set,
			iBase_FACE, iMesh_ALL_TOPOLOGIES,
			&ents, &ents_alloc, &ents_size, &err);
    }
    else {
      iMesh_getEntities(impl, root_set,
			iBase_REGION, iMesh_ALL_TOPOLOGIES,
			&ents, &ents_alloc, &ents_size, &err);
    }
    ERRORR("Trouble getting entities for specifying neumannsets via skinner.", err);

    // get tag handle
    const char *tag_neumann1 = "NEUMANN_SET";
    const char *global_id1 = "GLOBAL_ID";

    iMesh_getTagHandle(impl, tag_neumann1, &ntag1, &err, 12);
    ERRORR("Trouble getting handle.", err);

    iMesh_getTagHandle(impl, global_id1, &gtag1, &err, 9);
    ERRORR("Trouble getting handle.", err);

    iMesh_createEntSet(impl,0, &set, &err);
    ERRORR("Trouble creating set handle.", err);

    if (set_DIM == 3) { // sets for collecting top and bottom surface
      iMesh_createEntSet(impl,0, &set_z1, &err);
      ERRORR("Trouble creating set handle.", err);

      iMesh_createEntSet(impl,0, &set_z2, &err);
      ERRORR("Trouble creating set handle.", err);
    }

    MBRange tmp_elems;
    tmp_elems.insert((MBEntityHandle*)ents, (MBEntityHandle*)ents + ents_size);

    // get the skin of the entities
    MBSkinner skinner(mbImpl());
    MBRange skin_range;
    MBErrorCode result;
    MBRange::iterator rit;

    result = skinner.find_skin(tmp_elems, set_DIM-1, skin_range);
    if (MB_SUCCESS != result) return result;

    for (rit = skin_range.begin(), i = 0; rit != skin_range.end(); rit++, i++) {
      if(set_DIM == 3) { // filter top and bottom
	int num_vertex=0, size_vertex =0;
	iBase_EntityHandle *vertex = NULL;
	iMesh_getEntAdj(impl, (iBase_EntityHandle)(*rit), iBase_VERTEX, &vertex,
			&num_vertex, &size_vertex, &err);
	ERRORR("Trouble getting number of entities after merge.", err);

	double *coords = NULL;
	int coords_alloc = 0, coords_size=0;
	iMesh_getVtxArrCoords(impl, vertex, size_vertex, iBase_INTERLEAVED,
			      &coords, &coords_alloc, &coords_size, &err);
	ERRORR("Trouble getting number of entities after merge.", err);

	double ztemp = coords[2];
	int flag = 0;
	for (int p=1; p<num_vertex; p++) {
	  double z1 = coords[3*p+2];
	  if( ztemp != z1) {
	    flag = 1;
	    continue;
	  }
	}
	if(flag == 0) { // this is top or bottom surface
	  if (z_flag == 0) { // store z-coord this is the first call
	    z_flag = 1;
	    z1 = ztemp;
	  }
	  if( ztemp == z1) {
	    iMesh_addEntToSet(impl, (iBase_EntityHandle)(*rit), set_z1, &err);
	    ERRORR("Trouble getting number of entities after merge.", err);
	  }
	  else {
	    z2 = ztemp;
	    iMesh_addEntToSet(impl, (iBase_EntityHandle)(*rit), set_z2, &err);
	    ERRORR("Trouble getting number of entities after merge.", err);
	  }
	}
	else if (flag == 1) { // add the faces that are not top or bottom surface
	  iMesh_addEntToSet(impl, (iBase_EntityHandle)(*rit), set, &err);
	  ERRORR("Trouble getting number of entities after merge.", err);
	}
      }
      else if(set_DIM == 2) { // edges add all for sideset
	iMesh_addEntToSet(impl, (iBase_EntityHandle)(*rit), set, &err);
	ERRORR("Trouble getting number of entities after merge.", err);
      }
    }

    iMesh_setEntSetIntData( impl, set, ntag1, nss_Id, &err);
    ERRORR("Trouble getting handle.", err);

    iMesh_setEntSetIntData( impl, set, gtag1, nss_Id, &err);
    ERRORR("Trouble getting handle.", err);

    if (set_DIM == 3) {
      if (nst_flag == true || nsb_flag == true) {

	iMesh_setEntSetIntData( impl, set_z1, ntag1, nst_Id, &err);
	ERRORR("Trouble getting handle.", err);

	iMesh_setEntSetIntData( impl, set_z1, gtag1, nst_Id, &err);
	ERRORR("Trouble getting handle.", err);

	iMesh_setEntSetIntData( impl, set_z2, ntag1, nsb_Id, &err);
	ERRORR("Trouble getting handle.", err);

	iMesh_setEntSetIntData( impl, set_z2, gtag1, nsb_Id, &err);
	ERRORR("Trouble getting handle.", err);
      }
    }
#endif
  }
  return iBase_SUCCESS;
}

