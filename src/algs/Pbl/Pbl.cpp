#include "../meshkit/Pbl.hpp"
#ifdef HAVE_MOAB
#include "MBSkinner.hpp"
#include "MBAdaptiveKDTree.hpp"
#include "MBRange.hpp"
#include "MBCartVect.hpp"
#endif
namespace MeshKit
{
  bool debug =true;
  // static registration of this mesh scheme
  moab::EntityType Pbl_tps[] = { moab::MBTRI,
				 moab::MBHEX,
				 moab::MBMAXTYPE};
  const moab::EntityType* Pbl::output_types()
  { return Pbl_tps; }

  Pbl::Pbl( MKCore *mk, const MEntVector &me_vec)
    : MeshScheme( mk, me_vec),
      igeom(mk->igeom_instance()), imesh(mk->imesh_instance()), 
      mb (mk->moab_instance())
      // ---------------------------------------------------------------------------
      // Function: Obtains parameters for post meshing boundary layer and load mesh file
      // Input:    command line arguments
      // Output:   none
      // ---------------------------------------------------------------------------
  {
    m_SurfId = -1;
    m_NeumannSet = -1;
    m_nLineNumber = 0;
    szComment = "!";
    MAXCHARS = 300;
    err = 0;
  }

  Pbl::~Pbl()
  // ---------------------------------------------------------------------------
  // Function: Obtains parameters for post meshing boundary layer and load mesh file
  // Input:    command line arguments
  // Output:   none
  // ---------------------------------------------------------------------------
  {}

  bool Pbl::add_modelent(ModelEnt *model_ent)
  // ---------------------------------------------------------------------------
  // Function: Obtains parameters for post meshing boundary layer and load mesh file
  // Input:    command line arguments
  // Output:   none
  // ---------------------------------------------------------------------------
  {
    return MeshOp::add_modelent(model_ent);
  }

  void Pbl::setup_this()
  // ---------------------------------------------------------------------------
  // Function: Obtains parameters for post meshing boundary layer and load mesh file
  // Input:    command line arguments
  // Output:   none
  // ---------------------------------------------------------------------------
  {

    if (debug) {
      std::cout << "\nIn setup this : " <<  "\n";
    }

  }

  void Pbl::execute_this()
  // ---------------------------------------------------------------------------
  // Function: Obtains parameters for post meshing boundary layer and load mesh file
  // Input:    command line arguments
  // Output:   none
  // ---------------------------------------------------------------------------
  {
    // start the timer 
    CClock Timer;
    clock_t sTime = clock();
    std::string szDateTime;
    Timer.GetDateTime (szDateTime);
    std::cout << "\nStarting out at : " << szDateTime << "\n";
    
    std::cout << "\n Loading meshfile: " << m_MeshFile << ".." << std::endl;
    //load specified mesh file
    IBERRCHK(imesh->load(0, m_MeshFile.c_str(),0), *imesh);
    m_GD = imesh->getGeometricDimension();
    std::cout  <<"Geometric dimension of meshfile = "<< m_GD <<std::endl;
    if (debug) {
      std::cout << "\nIn execute this : " <<  "\n";
    }
    // obtain the boundary layer surface faces
    moab::Tag GDTag, GIDTag, NTag, MTag;
    MBERRCHK(mb->tag_get_handle("GEOM_DIMENSION", 1, moab::MB_TYPE_INTEGER, GDTag),mb);
    MBERRCHK(mb->tag_get_handle("NEUMANN_SET", 1, moab::MB_TYPE_INTEGER, NTag),mb);
    MBERRCHK(mb->tag_get_handle("MATERIAL_SET", 1, moab::MB_TYPE_INTEGER, MTag),mb);
    MBERRCHK(mb->tag_get_handle("GLOBAL_ID", 1, moab::MB_TYPE_INTEGER, GIDTag),mb);

    moab::Range sets, verts, n_sets;
    int geom_dim = 2;
    const void* gdim[] = {&geom_dim};
    
    // get all the entity sets with Geom Dimension=2
    MBERRCHK(mb->get_entities_by_type_and_tag(0, moab::MBENTITYSET, &GDTag, gdim, 1 , sets, moab::Interface::INTERSECT, false), mb);

    MBERRCHK(mb->get_entities_by_type_and_tag(0, moab::MBENTITYSET, &NTag, 0, 1 , n_sets),mb);

    moab::Range::iterator set_it;
    moab::EntityHandle this_set;
    for (set_it = n_sets.begin(); set_it != n_sets.end(); set_it++)  {
      
      this_set = *set_it;
      
      // get the id for this set
      int set_id;
      MBERRCHK(mb->tag_get_data(NTag, &this_set, 1, &set_id), mb);
      if(set_id == m_NeumannSet)
	break;
      this_set = NULL;
    }
    if (debug)
      std::cout << "Looking for NS with id " << m_NeumannSet << ". Total NS found are: "<< n_sets.size() << std::endl;

    
    // get the quads for the surface with input global ids
    moab::EntityHandle s1;
    MBRange quads, nodes;
    int dims; // variable to store global id of all sets with GD=2


    if(m_NeumannSet != -1){
      MBERRCHK(mb->get_entities_by_type(this_set, moab::MBQUAD, quads),mb);
      MBERRCHK(mb->get_adjacencies(quads, 0, false, nodes, MBInterface::UNION),mb);
      if (debug) {
	std::cout << "Found NeumannSet with id : " << m_NeumannSet <<  std::endl;
	std::cout << "#Quads in this surface: " << quads.size() << std::endl;
	std::cout << "#Nodes in this surface: " << nodes.size() << std::endl;
      }
    }
    else{
      // from the input surface get the quads and nodes in a range
      for(MBRange::iterator rit=sets.begin(); rit != sets.end(); ++rit){
	s1 = *rit;
	MBERRCHK(mb->tag_get_data(GIDTag, &s1, 1, &dims),mb);

	if(dims == m_SurfId && m_SurfId != -1){
	  MBERRCHK(mb->get_entities_by_type(s1, moab::MBQUAD, quads),mb);
	  MBERRCHK(mb->get_adjacencies(quads, 0, false, nodes, MBInterface::UNION),mb);
	  if (debug) {
	    std::cout << "Found surface with id : " << m_SurfId <<  std::endl;
	    std::cout << "#Quads in this surface: " << quads.size() << std::endl;
	    std::cout << "#Nodes in this surface: " << nodes.size() << std::endl;
	  }
	}
      }
    }
 

    //call extrude and extrude these quads 
    std::vector <bool> node_status(false); // size of verts of bl surface
    node_status.resize(nodes.size());
    MBRange edges, hexes, hex_edge, quad_verts, adj_quads;  
    std::vector<EntityHandle> conn(m_Intervals*8), qconn(4), new_vert(m_Intervals*nodes.size()), old_hex, old_hex_conn(8);
    double coords_bl_quad[3], coords_new_quad[3], xdisp = 0.0, ydisp = 0.0, zdisp = 0.0;
    EntityHandle hex;
    int qcount = 0;
    int ncount = 0;
    for (Range::iterator kter = quads.begin(); kter != quads.end(); ++kter){
      qcount++;
      MBERRCHK(mb->get_connectivity(&(*kter), 1, qconn),mb);
      
      // MBERRCHK(mb->get_adjacencies(&quads[qcount], 1, false, adj_quads, MBInterface::UNION),mb);
    
      get_normal_quad (qconn, xdisp, ydisp, zdisp); 
      if (debug) 
	std::cout << "\n\n*** QUAD: " << qcount << "\n\ndirection vector: " << xdisp << " " << ydisp << " " << zdisp << " " << std::endl;
      //adj_quads.clear();
      MBERRCHK(mb->get_adjacencies(&(*kter), 1, 3, false, old_hex),mb);
      MBERRCHK(mb->get_connectivity(&old_hex[0], 1, old_hex_conn),mb);
      
      for (int i=0; i<4; i++){
	MBERRCHK(mb->get_coords(&qconn[i], 1, coords_bl_quad),mb);
	
	// check to see if this node is already dealt with
	int tmp = 0;
	for(int n=0; n< (int) nodes.size(); n++){
	  if(nodes[n] == qconn[i]){
	    tmp = n;
	  }
	}
	if (debug) 
	  std::cout << "Working on Node : " << tmp << std::endl;
	if(node_status[tmp] == false){

	  // another loop to number of boundary layers
	  for(int j=0; j< m_Intervals; j++){
	    
	    // now compute the coords of the new vertex
	    coords_new_quad[0] = coords_bl_quad [0]-(j+1)*m_Thickness*xdisp;
	    coords_new_quad[1] = coords_bl_quad [1]-(j+1)*m_Thickness*ydisp;
	    coords_new_quad[2] = coords_bl_quad [2]-(j+1)*m_Thickness*zdisp;
	    
	    // TODO: Check to see if this vertex is possible?
	    int nid = tmp*m_Intervals+j;
	    mb->create_vertex(coords_new_quad, new_vert[nid]);
	    if (debug) 
	      std::cout << tmp << ": created vert:" << nid << " of " << new_vert.size() << std::endl;
	    if(j==0) // set connectivity of boundary layer hex
	      conn[8*j + i+4] = qconn[i];
	    else
	      conn[8*j + i+4] = new_vert[nid-1];
	    conn[8*j + i] = new_vert[nid];
	  }
	  node_status[tmp] = true;
	  ncount++;
	}
	else{
	  for(int j=0; j< m_Intervals; j++){
	    
	    int nid = tmp*m_Intervals+j;
	    if(j==0) // set connectivity of boundary layer hex
	      conn[8*j + i+4] = qconn[i];
	    else
	      conn[8*j + i+4] = new_vert[nid-1];
	    conn[8*j +i] = new_vert[nid];
	  }
	}
      }
      if (debug) 
	std::cout << "set connectivity of the old_hex " << std::endl;
      // set connectivity of the old_hex
      for(int p=0; p<8; p++){
	for(int q=0; q<4; q++){ 
	  if (old_hex_conn[p] == qconn[q]){
	    old_hex_conn[p] = conn[8*(m_Intervals-1) + q];
	  }
	}
      }
      
      MBERRCHK(mb->set_connectivity(old_hex[0], &old_hex_conn[0], 8), mb); 
      old_hex.clear();
      if (debug) 
	std::cout << "create boundary layer hexes" << std::endl;
      // create boundary layer hexes
      for(int j=0; j< m_Intervals; j++){
	MBERRCHK(mb->create_element(MBHEX, &conn[j*8], 8, hex),mb);
      }
    }
    
    //save the final boundary layer mesh	       
    MBERRCHK(mb->write_mesh(m_OutFile.c_str()),mb);
    std::cout << "Wrote Mesh File: " << m_OutFile << std::endl;  
    // get the current date and time
    Timer.GetDateTime (szDateTime);
    std::cout << "Ending at : " << szDateTime;
 
    // compute the elapsed time
    std::cout << "Elapsed wall clock time: " << Timer.DiffTime ()
	      << " seconds or " << (Timer.DiffTime ())/60.0 << " mins\n";

    std::cout << "Total CPU time used: " << (double) (clock() - sTime)/CLOCKS_PER_SEC \
	      << " seconds" << std::endl; 
  }

  void Pbl::PrepareIO (int argc, char *argv[], std::string  TestDir)
  // ---------------------------------------------------------------------------
  // Function: Obtains parameters for post meshing boundary layer and load mesh file
  // Input:    command line arguments
  // Output:   none
  // ---------------------------------------------------------------------------
  {
    std::cout << '\n';
    std::cout << "\t\t---------------------------------------------------------" << '\n';
    std::cout << "\t         Tool to generate postmesh boundary layers      " << '\n';
    std::cout << "\t\t\t\tArgonne National Laboratory" << '\n';
    std::cout << "\t\t\t\t        2012         " << '\n';
    std::cout << "\t\t---------------------------------------------------------" << '\n';
    std::cout << "\nsee README file for using the program and details on various cards.\n"<< std::endl;

    // set and open input output files
    bool bDone = false;
    do{
      if (2 == argc) {
	m_InputFile = argv[1];
      }
      else if (1 == argc){
	std::cout << "\nRunning default case:\n" << std::endl;
	m_InputFile = TestDir + "/" + (char *)DEFAULT_TEST_PBL;
      }

      // open the file
      m_FileInput.open (m_InputFile.c_str(), std::ios::in); 
      if (!m_FileInput){
	std::cout << "Unable to open file: " << m_InputFile << std::endl;
	m_FileInput.clear ();
	exit(1);
      }
      else
	bDone = true; // file opened successfully
    }while (!bDone);

    // Get the meshfile name, surface(s), thickness, intervals and bias
    CParser Parse;
 
    // count the total number of cylinder commands in each pincellh
    for(;;){
      if (!Parse.ReadNextLine (m_FileInput, m_nLineNumber, szInputString,
			       MAXCHARS, szComment))
	IOErrorHandler (INVALIDINPUT);

      // Get MeshFile name
      if (szInputString.substr(0,8) == "meshfile"){
	std::istringstream szFormatString (szInputString);
	szFormatString >> m_Card >> m_MeshFile;
	if(szFormatString.fail())
	  IOErrorHandler(INVALIDINPUT);
	std::cout << m_Card << " name read: "<< m_MeshFile << std::endl;
	if (argc == 1){
	  m_MeshFile = TestDir + "/" + m_MeshFile;
	}
      }
      // Get BL surface 
      if (szInputString.substr(0,8) == "surfaces"){
	std::istringstream szFormatString (szInputString);
	szFormatString >> m_Card >> m_SurfId;
	if(m_SurfId < 0 || szFormatString.fail())
	  IOErrorHandler(INVALIDINPUT);
	std::cout  << m_Card << " read: "<< m_SurfId <<std::endl;
      }
      // Get BL surface via neumann set or sideset
      if (szInputString.substr(0,10) == "neumannset"){
	std::istringstream szFormatString (szInputString);
	szFormatString >> m_Card >> m_NeumannSet;
	if(m_NeumannSet < 0 || szFormatString.fail())
	  IOErrorHandler(INVALIDINPUT);
	std::cout  << m_Card << " read: "<< m_NeumannSet <<std::endl;
      }
 
      // Get thickness
      if (szInputString.substr(0,9) == "thickness"){
	std::istringstream szFormatString (szInputString);
	szFormatString >> m_Card >> m_Thickness;
	if(m_Thickness < 0 || szFormatString.fail())
	  IOErrorHandler(INVALIDINPUT);
	std::cout  << m_Card << " read: "<< m_Thickness <<std::endl;
      }
      // Get intervals
      if (szInputString.substr(0,9) == "intervals"){
	std::istringstream szFormatString (szInputString);
	szFormatString >> m_Card >> m_Intervals;
	if(m_Intervals < 0 || szFormatString.fail())
	  IOErrorHandler(INVALIDINPUT);
	std::cout  << m_Card << " read: "<< m_Intervals <<std::endl;
      }
      // Get bias
      if (szInputString.substr(0,4) == "bias"){
	std::istringstream szFormatString (szInputString);
	szFormatString >> m_Card >> m_Bias;
	if(m_Bias < 0 || szFormatString.fail())
	  IOErrorHandler(INVALIDINPUT);
	std::cout  << m_Card << " read: "<< m_Bias <<std::endl;
      }
      // Output file name
      if (szInputString.substr(0,7) == "outfile"){
	std::istringstream szFormatString (szInputString);
	szFormatString >> m_Card >> m_OutFile;
	if(m_Bias < 0 || szFormatString.fail())
	  IOErrorHandler(INVALIDINPUT);
	std::cout  << m_Card << " name read: "<< m_OutFile <<std::endl;
      }
      if (szInputString.substr(0,3) == "end"){
	break;
      }
    }
  }
  void Pbl::IOErrorHandler (ErrorStates ECode) const
  // ---------------------------------------------------------------------------
  // Function: displays error messages related to input data
  // Input:    error code
  // Output:   none
  // ---------------------------------------------------------------------------
  {
    std::cerr << '\n';
    if (ECode == INVALIDINPUT) // invalid input
      std::cerr << "Invalid input.";
    else
      std::cerr << "Unknown error ...?";

    std::cerr << '\n' << "Error in input file line : " << m_nLineNumber;
    std::cerr << std::endl;
    exit (1);
  }
  void Pbl::get_normal_quad (std::vector<EntityHandle>conn, double &x, double &y, double &z)
  // ---------------------------------------------------------------------------
  // Function: get normal of a quad
  // Input:    conn
  // Output:   vector x, y and z
  // ---------------------------------------------------------------------------
  {
    CartVect coords[3];
    MBERRCHK(mb->get_coords(&conn[0], 3, (double*) &coords[0]), mb);
    CartVect AB(coords[1] - coords[0]);
    CartVect BC(coords[2] - coords[1]);
    CartVect normal = AB*BC;
    normal = normal/normal.length();
    x = normal[0];
    y = normal[1];
    z = normal[2];
  }
} // namespace MeshKit


