/*********************************************
June,10
Reactor Assembly Mesh Assembler
Argonne National Laboratory

CCrgen class definition.
*********************************************/
#define DEFAULT_TEST_FILE "twoassm"
#include "crgen.hpp"

/* ==================================================================
   ======================= CCrgen class =============================
   ================================================================== */

CCrgen::CCrgen ()
// ---------------------------------------------------------------------------
// Function: default constructor
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
  err = 0;
  UNITCELL_DUCT = 0;
  ASSY_TYPES = 1;
  pack_type = 1; symm = 1;
  z_height = 1;
  z_divisions = 2;
  set_DIM = 3; // default is 3D
  PI = acos(-1.0);
  comment = "!";
  MAXCHARS = 300;
  global_ids = true;
  merge_tol =  1.0e-3;
  do_merge = 1;
  update_sets= 0;
  merge_tag = NULL;

}

CCrgen::~CCrgen ()
// ---------------------------------------------------------------------------
// Function: destructor
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
}

int CCrgen::close ()
// ---------------------------------------------------------------------------
// Function: dellocating 
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
  // deallocate ... deallocate ... deallocate
  for (unsigned int i = 0; i < files.size(); i++) {
    delete cm[i];
  } 

  iMesh_dtor(impl, &err);
  ERRORR("Failed in call iMesh_dtor", err);

  return 0;
}



int CCrgen::prepareIO (int argc, char *argv[])
// -----------------------------------------------------------------------------------
// Function: Obtains file names and opens input/output files and then read/write them
// Input:    command line arguments
// Output:   none
// -----------------------------------------------------------------------------------
{
  bool bDone = false;
  do{
    if (2 == argc) {
      iname = argv[1];
      ifile = iname+".inp";
      outfile = iname+".h5m";
      mfile = iname + ".makefile";
    }
    else if (3 == argc) {
      int i=1;// will loop through arguments, and process them
      for (i=1; i<argc-1 ; i++) {
	if (argv[i][0]=='-') {
	  switch (argv[i][1]) {
	  case 'm': {
	    std::cout << "Creating Makfile Only" << std::endl;
	    // only makefile creation specified
	    iname = argv[2];
	    ifile = iname+".inp";
	    outfile = iname+".h5m";
	    mfile = iname + ".makefile";
	  }
	  }
	}
      }
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
      std::cout << "Usage: coregen <input file WITHOUT EXTENSION>"<< std::endl; 
      file_input.clear ();
      exit(1);
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


  // now call the functions to read and write
  err = read_inputs_phase1();
  ERRORR("Failed to read inputs in phase1.", 1);

  err = read_inputs_phase2();
  ERRORR("Failed to read inputs in phase2.", 1);
 
  err = write_makefile();
  ERRORR("Failed to write a makefile.", 1);

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

  mm = new MergeMesh(impl);
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
 
  return iBase_SUCCESS;
}

int CCrgen::read_inputs_phase1 (){
  // ---------------------------------------------------------------------------
  // Function: Reads the dimension and symmetry of the problem
  // Input:    none
  // Output:   none
  // ---------------------------------------------------------------------------
  CParser parse;
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

    // geom engine
    if (input_string.substr(0,10) == "geomengine"){
      std::istringstream formatString (input_string);
      formatString >> card >> geom_engine;    
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

int CCrgen::read_inputs_phase2 ()
// ---------------------------------------------------------------------------
// Function: read all the inputs 
// Input:    command line arguments
// Output:   none
// ---------------------------------------------------------------------------
{
  //Rewind the input file
  file_input.clear (std::ios_base::goodbit);
  file_input.seekg (0L, std::ios::beg);
  linenumber = 0;

  CParser parse;
  for(;;){
    if (!parse.ReadNextLine (file_input, linenumber, input_string, 
			     MAXCHARS, comment))
      ERRORR("Reading input file failed",1);

    if (input_string.substr(0,12) == "geometrytype"){
      std::istringstream formatString (input_string);
      formatString >> card >> geom_type;
      if(geom_type == "hexvertex" && symm == 6){

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

      else if(geom_type == "rectangular" && symm ==1){

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

      else if (geom_type =="hexflat" && symm == 6){
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
      else if (geom_type =="hexflat" && symm == 1){
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
	  tot_assys = 3 * (nrings * (nrings - 1)) + 1;
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
      else if(geom_type == "hexflat" && symm == 12){

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
	    tot_assys = (nrings*(nrings + 2))/4;
	  else
	    tot_assys = ((nrings+1) * (nrings+1))/4;
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
    // z-height
    if (input_string.substr(0,8) == "z-height"){
      std::istringstream formatString (input_string);
      formatString >> card >> z_height;
    }

    // z-divisions
    if (input_string.substr(0,11) == "z-divisions"){
      std::istringstream formatString (input_string);
      formatString >> card >> z_divisions;
    }

    // OutputFileName
    if (input_string.substr(0,14) == "outputfilename"){
      std::istringstream formatString (input_string);
      formatString >> card >> outfile;
    }

    // breaking condition
    if(input_string.substr(0,3) == "end"){
      std::istringstream formatstring (input_string);
      break;
    }
  }
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


int CCrgen::banner ()
// ---------------------------------------------------------------------------
// Function: display the program banner
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{ 
  std::cout << '\n';
  std::cout << "\t\t--------------------------------------------------------------------" << '\n';
  std::cout << "\t\tProgram to Assemble Nuclear Reactor Assembly Meshes and Form a Core     " << '\n';
  std::cout << "\t\t\t\t\tArgonne National Laboratory" << '\n';
  std::cout << "\t\t\t\t\t        2009-2010         " << '\n';
  std::cout << "\t\t--------------------------------------------------------------------" << '\n';
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
  make_file << "## This makefile is automatically generated by coregen program" << std::endl;
  make_file << "##" << std::endl;
  make_file << "## Check your coregen, assygen and cubit location" << std::endl;
  make_file << "##" << std::endl;
  make_file << "\nCUBIT = cubit\n" << std::endl;
  make_file << "COREGEN = ../../coregen\n" << std::endl;
  make_file << "ASSYGEN = ../../assygen\n" << std::endl;

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
    if(geom_engine == "occ")
      name = f_no_ext[i] + ".brep";
    else
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

  make_file << "\n\n" << outfile << " : ${MESH_FILES} " << ifile <<  std::endl;
  make_file << "\t" << "${COREGEN} " << iname << std::endl;
  for(unsigned int i=0; i<files.size(); i++){
    make_file << files[i] << " : " << f_sat[i] << "  " << f_jou[i] << "  " << f_injou[i] << std::endl;
    make_file << "\t" << "${CUBIT} -batch" << f_jou[i] <<"\n" << std::endl;

    make_file << f_sat[i] << " " << f_jou[i] << " " << f_injou[i] << " : " << f_inp[i] << std::endl;
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

  int verts_ents_alloc = 0, verts_ents_size;
  iBase_EntityHandle *verts_ents = NULL;

  iMesh_getEntities(impl, set, iBase_VERTEX,iMesh_ALL_TOPOLOGIES,
		    &verts_ents, &verts_ents_alloc, &verts_ents_size, &err);
  ERRORR("Failed to get any entities from original set.", iBase_FAILURE);


  double *coords=0;
  int coords_alloc=0,coords_size=0;

  iMesh_getVtxArrCoords(impl, verts_ents, verts_ents_size, iBase_INTERLEAVED, &coords,
			&coords_alloc, &coords_size, &err); 
  ERRORR("Failed to get vtx coords from set.", iBase_FAILURE);

  for(int i = 0; i < verts_ents_size; i++){
    coords[3*i]   += dx[0];
    coords[3*i+1] += dx[1];
    coords[3*i+2] += dx[2];
  }
  
  iMesh_setVtxArrCoords(impl, verts_ents, verts_ents_size, iBase_INTERLEAVED, coords,
		        coords_size, &err);  
  ERRORR("Failed to set vtx coords.", iBase_FAILURE);

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
  
  int num1, num2;

  iMesh_getNumOfType(impl, root_set, iBase_VERTEX, &num1, &err);
  ERRORR("Trouble getting number of entities before merge.", err);

  // merge now
  mm->merge_entities(ents, ents_size, merge_tol, do_merge, update_sets,
                     merge_tag);
  
  iMesh_getNumOfType(impl, root_set, iBase_VERTEX, &num2, &err);
  ERRORR("Trouble getting number of entities after merge.", err);

  std::cout << "Merged " << num1 - num2 << " vertices." << std::endl;
  return iBase_SUCCESS;
}

int CCrgen::assign_gids()
{
  // assign new global ids
  if (global_ids == true){
    std::cout << "Assigning global ids." << std::endl;
    mu = new MKUtils(impl);
    err = mu->assign_global_ids(root_set, 3, 1, true, false,
				"GLOBAL_ID");
    ERRORR("Error assigning global ids.", err);
  }
  return iBase_SUCCESS;
}


int CCrgen::save()
{
  // export
  iMesh_save(impl, root_set, outfile.c_str(), NULL, &err, 
             strlen(outfile.c_str()), 0);
  ERRORR("Trouble writing output mesh.", err);
  std::cout << "Saved: "<< outfile.c_str() <<std::endl;

  return iBase_SUCCESS;
}



int CCrgen::extrude()
{
  // extrude if this is a surface mesh
  if(set_DIM ==2){ // if surface geometry specified
    std::cout << "Extruding surface mesh." << std::endl;

    //get entities for merging
    iBase_EntityHandle *ents = NULL; 
    int ents_alloc = 0, ents_size;
    int err = 0;


    iMesh_getEntities(impl, root_set,
		      iBase_FACE, iMesh_ALL_TOPOLOGIES,
		      &ents, &ents_alloc, &ents_size, &err);
    ERRORR("Trouble getting face mesh.", err);

    double v[] = { 0, 0, z_height };
    int steps = z_divisions;

    ExtrudeMesh *ext = new ExtrudeMesh(impl);
    ext->extrude(ents, ents_size, extrude::Translate(v, steps));
    ERRORR("Trouble extruding mesh.", err);
  }
  return iBase_SUCCESS;
}
