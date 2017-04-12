#include "meshkit/PostBL.hpp"

namespace MeshKit
{
  
  // static registration of this mesh scheme
  moab::EntityType PostBL_tps[] = {moab::MBHEX,
                                   moab::MBMAXTYPE};
  const moab::EntityType* PostBL::output_types()
  { return PostBL_tps; }
  
  PostBL::PostBL( MKCore *mk, const MEntVector &me_vec)
    : MeshScheme( mk, me_vec),
      igeom(mk->igeom_instance()), imesh(mk->imesh_instance()),
      mb (mk->moab_instance())
    // ---------------------------------------------------------------------------
    //! Function: Constructor \n
    //! Input:    Initialize mesh and geometry instances and parameters \n
    //! Output:   none
    // ---------------------------------------------------------------------------
  {
    tri_sch = 2;
    m_Conn = 0;
    m_BElemNodes = 0;
    m_SurfId = -1;
    check_bl_edge_length = false;
    debug = false;
    hybrid = false;
    m_NeumannSet = -1;
    m_Material = 999999;
    m_nLineNumber = 0;
    szComment = "!";
    MAXCHARS = 300;
    m_JacCalls = 0;
    m_JLo = 0.0;
    m_JHi = 0.0;
    err = 0;
    fixmat = -1;
    hex27 = 0;
//    hex = NULL;
//    hex1 = NULL;
//    hex2 = NULL;
  }
  
  PostBL::~PostBL()
  // ---------------------------------------------------------------------------
  //! Function: Destructor, does nothing..\n
  //! Input:    none \n
  //! Output:   none \n
  // ---------------------------------------------------------------------------
  {}
  
  bool PostBL::add_modelent(ModelEnt *model_ent)
  // ---------------------------------------------------------------------------
  //! Function: Adds entities for PosBL graph node.\n
  //! Input:    ModelEnt \n
  //! Output:   none \n
  // ---------------------------------------------------------------------------
  {
    return MeshOp::add_modelent(model_ent);
  }
  
  void PostBL::setup_this()
  // ---------------------------------------------------------------------------
  //! Function: Setup the graph node for PostBL \n
  //! Input:    none \n
  //! Output:   none \n
  // ---------------------------------------------------------------------------
  {
    if (debug) {
        m_LogFile <<  "\nIn setup this : " <<  std::endl;
      }
    
  }
  
  void PostBL::execute_this()
  // ---------------------------------------------------------------------------
  //! Function: Read user input from file and run the PostBL algorithm \n
  //! Input:     Uses the file name (.inp) with keywords predefined by PosBL algorithm. \n
  //! Output:    Resulting mesh file is saved. \n
  // ---------------------------------------------------------------------------
  {
    m_LogFile << "\nIn execute this : creating boundary layer elements.." <<  std::endl;
    // start the timer
    CClock Timer;
    clock_t sTime = clock();
    std::string szDateTime;
    Timer.GetDateTime (szDateTime);
    VerdictWrapper vw(mb);
    
    m_LogFile <<  "\nStarting out at : " << szDateTime << std::endl;
    m_LogFile <<  "\n Loading meshfile: " << m_MeshFile << std::endl;
    
    // load specified mesh file
    IBERRCHK(imesh->load(0, m_MeshFile.c_str(),0), *imesh);
    
    //check if intervals is read from input file
    if(m_Intervals == 0){
        m_LogFile << "Please specify desired number of intervals using 'Intervals' keyword" << std::endl;
        exit(0);
      }
    // if no NS specified, take all d-1 elements as BL input
    if (m_NeumannSet == -1 && m_SurfId == -1) {
        // create a neumann set with id 99999 in the model and set it for BL generation
        moab::Range neuEnts;
        moab::EntityHandle neuEntSet;
        moab::Tag neuTag;
        mb->create_meshset(moab::MESHSET_SET, neuEntSet);
        mb->get_entities_by_dimension(mb->get_root_set(), 2, neuEnts, true);
        mb->add_entities(neuEntSet, neuEnts);
        mb->tag_get_handle("NEUMANN_SET", neuTag);
        m_NeumannSet = 999999;
        mb->tag_set_data(neuTag, &neuEntSet, 1, (void*) &m_NeumannSet);
      }
    
    moab::Range all_elems, all_verts;
    MBERRCHK(mb->get_entities_by_dimension(0, 3, all_elems,true),mb);
    if (all_elems.size() == 0)
      m_GD = 2;
    else if (all_elems.size() > 0)
      m_GD = 3;
    else
      exit(0);
    all_elems.clear();
    m_LogFile << "Geometric dimension of meshfile = "<< m_GD <<std::endl;
    
    // obtain existing tag handles
    MBERRCHK(mb->tag_get_handle("GEOM_DIMENSION", 1, moab::MB_TYPE_INTEGER, GDTag),mb);
    MBERRCHK(mb->tag_get_handle("NEUMANN_SET", 1, moab::MB_TYPE_INTEGER, NTag),mb);
    MBERRCHK(mb->tag_get_handle("MATERIAL_SET", 1, moab::MB_TYPE_INTEGER, MTag),mb);
    MBERRCHK(mb->tag_get_handle("GLOBAL_ID", 1, moab::MB_TYPE_INTEGER, GIDTag),mb);
    // create smoothset and fixed tag for mesquite
    MBERRCHK(mb->tag_get_handle("SMOOTHSET", 1, moab::MB_TYPE_INTEGER, STag,
                                moab::MB_TAG_SPARSE|moab::MB_TAG_CREAT),mb);
    MBERRCHK(mb->tag_get_handle("fixed", 1, moab::MB_TYPE_INTEGER, FTag,
                                moab::MB_TAG_SPARSE|moab::MB_TAG_CREAT),mb);
    MBERRCHK(mb->tag_get_handle("mnode", 1, moab::MB_TYPE_INTEGER, MNTag,
                                moab::MB_TAG_SPARSE|moab::MB_TAG_CREAT),mb);
    MBERRCHK(mb->tag_get_handle("matid", 1, moab::MB_TYPE_INTEGER, MatIDTag,
                                moab::MB_TAG_SPARSE|moab::MB_TAG_CREAT),mb);
    MBERRCHK(mb->tag_get_handle("BLNODEID", 1, moab::MB_TYPE_INTEGER, BLNodeIDTag,
                                moab::MB_TAG_SPARSE|moab::MB_TAG_CREAT),mb);
    
    // get all the entity sets with boundary layer geom dimension, neumann sets and material sets
    m_BLDim = m_GD - 1;
    
    const void* gdim[] = {&m_BLDim};
    
    MBERRCHK(mb->get_entities_by_type_and_tag(0, moab::MBENTITYSET, &GDTag,
                                              gdim, 1 , sets, moab::Interface::INTERSECT, false), mb);
    MBERRCHK(mb->get_entities_by_type_and_tag(0, moab::MBENTITYSET, &NTag, 0, 1 , n_sets),mb);
    MBERRCHK(mb->get_entities_by_type_and_tag(0, moab::MBENTITYSET, &MTag, 0, 1 , m_sets),mb);
    
    // Handling NeumannSets (if BL surf in input via NS)
    moab::EntityHandle this_set = 0;
    for (set_it = n_sets.begin(); set_it != n_sets.end(); set_it++)  {
        
        this_set = *set_it;
        
        // get entity handle of NS specified in the input file
        int set_id;
        MBERRCHK(mb->tag_get_data(NTag, &this_set, 1, &set_id), mb);
        if(set_id == m_NeumannSet)
          break;
        this_set = 0;
      }
    if (debug && m_NeumannSet != -1 && this_set != 0){
        m_LogFile <<  "Looking for NS with id " << m_NeumannSet <<
                      ". Total NS found are: "<< n_sets.size() << std::endl;
      }
    
    // For specified surface: get the  all the quads and nodes in a range
    moab::EntityHandle s1;
    // variable to store global id of boundary layer specified in the input file
    int dims; 
    
    // Method 1: INPUT by NeumannSet
    if(m_NeumannSet != -1 && this_set != 0){
        MBERRCHK(mb->get_entities_by_dimension(this_set, m_BLDim, quads,true),mb);
        if (quads.size() <=0){
            m_LogFile <<  " No quads found, aborting.. " << std::endl;
            exit(0);
          }
        
        MBERRCHK(mb->get_adjacencies(quads, 0, false, nodes, moab::Interface::UNION),mb);
        if(m_GD == 3)
          MBERRCHK(mb->get_adjacencies(quads, 1, true, edges, moab::Interface::UNION),mb);
        
        if (debug) {
            m_LogFile <<  "Found NeumannSet with id : " << m_NeumannSet <<  std::endl;
            m_LogFile <<  "#Quads in this surface: " << quads.size() << std::endl;
            m_LogFile <<  "#Nodes in this surface: " << nodes.size() << std::endl;
            m_LogFile << "#New nodes to be created:" << m_Intervals*nodes.size() << std::endl;
          }
      }
    // Method 2: INPUT by surface id (geom dimension)
    else if (m_SurfId !=-1){
        for(moab::Range::iterator rit=sets.begin(); rit != sets.end(); ++rit){
            s1 = *rit;
            MBERRCHK(mb->tag_get_data(GIDTag, &s1, 1, &dims),mb);
            
            if(dims == m_SurfId && m_SurfId != -1){
                MBERRCHK(mb->get_entities_by_dimension(s1, m_BLDim, quads,true),mb);
                if (quads.size() <=0){
                    m_LogFile <<  " No quads found, aborting.. " << std::endl;
                    exit(0);
                  }
                
                MBERRCHK(mb->get_adjacencies(quads, 0, false, nodes, moab::Interface::UNION),mb);
                if(m_GD == 3)
                  MBERRCHK(mb->get_adjacencies(edges, 0, false, nodes, moab::Interface::UNION),mb);
                if (debug) {
                    m_LogFile <<  "Found surface with id : " << m_SurfId <<  std::endl;
                    m_LogFile <<  "#Quads in this surface: " << quads.size() << std::endl;
                    m_LogFile <<  "#Nodes in this surface: " << nodes.size() << std::endl;
                    m_LogFile << "#New nodes to be created:" << m_Intervals*nodes.size() << std::endl;
                  }
              }
          }
      }
    
    if (quads.size() == 0 || nodes.size() == 0) {
        m_LogFile <<  "Invalid boundary layer specification, aborting.." <<  std::endl;
        exit(0);
      }
    
    // placeholder for storing smoothing entities
    moab::EntityHandle smooth_set;
    int s_id = 100;
    MBERRCHK(mb->create_meshset(moab::MESHSET_SET, smooth_set, 1), mb);
    MBERRCHK(mb->tag_set_data(STag, &smooth_set, 1, &s_id), mb);
    
    // placeholder for storing gd on new entities
    MBERRCHK(mb->create_meshset(moab::MESHSET_SET, geom_set, 1), mb);
    MBERRCHK(mb->tag_set_data(GDTag, &geom_set, 1, &m_GD), mb);
    
    // declare variables before starting BL creation
    std::vector <bool> node_status(false); // size of verts of bl surface
    node_status.resize(nodes.size());
    blmaterial_id.resize(quads.size());
    
    //size of the following is based on element type  
    new_vert.resize(m_Intervals*nodes.size());
    // element on the boundary
    Range::iterator kter = quads.begin();
    
    MBERRCHK(mb->get_adjacencies(&(*kter), 1, m_GD, false, old_hex),mb);
    std::cout << "old_hex size is - before filtering fixmat - " << old_hex.size() << std::endl;
    if((int) old_hex.size() == 0){
        m_LogFile << "unable to find adjacent hex for BL quad, aborting...";
        exit(0);
      }
    
    // allocate space for connectivity/adjacency during the first pass of this loop
    if(mb->type_from_handle(old_hex[0]) == moab::MBHEX){
        m_Conn = 8;
        m_BElemNodes = 4;
        m_HConn = 8;
        //allocating based on element type
        conn.resize(m_Intervals*m_Conn), qconn.resize(m_BElemNodes), adj_qconn.resize(m_BElemNodes),
            old_hex_conn.resize(m_Conn), adj_hex_nodes1.resize(m_Conn);
      }
    else if(mb->type_from_handle(old_hex[0]) == MBTET){
        m_Conn = 4;
        m_BElemNodes = 3;
        //allocating based on element type - thrice the number of elements
        if(hybrid){
            m_HConn = 6;
            conn.resize(m_Intervals*6);
          }
        else{
            m_HConn = 6; // we use the prism for making tets
            tet_conn.resize(3*m_Intervals*m_Conn);
            conn.resize(m_Intervals*6);
          }
        qconn.resize(m_BElemNodes), adj_qconn.resize(m_BElemNodes),
            old_hex_conn.resize(m_Conn), adj_hex_nodes1.resize(m_Conn);
      }
    else if(mb->type_from_handle(old_hex[0]) == MBQUAD){
        m_Conn = 4;
        m_HConn = 4;
        m_BElemNodes = 2;
        //allocating based on element type
        conn.resize(m_Intervals*m_Conn), qconn.resize(m_BElemNodes), adj_qconn.resize(m_BElemNodes),
            old_hex_conn.resize(m_Conn), adj_hex_nodes1.resize(m_Conn);
      }
    else if(mb->type_from_handle(old_hex[0]) == MBTRI){
        m_Conn = 3;
        m_HConn = 4;
        m_BElemNodes = 2;
        //allocating based on element type - twice the number of elements
        if(hybrid){
            m_HConn = 4;
            conn.resize(m_Intervals*m_HConn);
          }
        else{
            tri_conn.resize(2*m_Intervals*m_Conn);
            conn.resize(2*m_Intervals*m_Conn);
          }
        qconn.resize(m_BElemNodes), adj_qconn.resize(m_BElemNodes),
            old_hex_conn.resize(m_Conn), adj_hex_nodes1.resize(m_Conn);
      }
    else if(m_Conn == 0 || m_BElemNodes == 0){
        m_LogFile << "This mesh type is not supported by this tool" << std::endl;
        exit(0);
      }
    
    // Tag all nodes on outer boundary with a unique number
    int node_id = 0;
    std::vector<int> NId(nodes.size());
    for(moab::Range::iterator nodes_iter = nodes.begin(); nodes_iter != nodes.end(); nodes_iter++){
        NId[node_id] = node_id;
        MBERRCHK(mb->tag_set_data(BLNodeIDTag, &(*nodes_iter),1, &NId[node_id]), mb);
        ++node_id;
      }
    
    // Handling MaterialSet
    int mset_id = 0, found = 0;
    for (mset_it = m_sets.begin(); mset_it != m_sets.end(); mset_it++)  {
        
        mthis_set = *mset_it;
        
        // get entity handle of MS specified in the input file
        MBERRCHK(mb->tag_get_data(MTag, &mthis_set, 1, &mset_id), mb);
        
        // if no material set is specified, we'll have to resolve else just set all the MNTag to 0
        if(mset_id == m_Material){
            found = 1;
            break;
          }
        else if(mset_id == fixmat){
            MBERRCHK(mb->get_entities_by_dimension(mthis_set, m_GD, fixmat_ents ,true),mb);
          }
        
        // get all the nodes in the material and tag bl nodes
        moab::Range mat_nodes, mat_hexes;
        if(m_GD == 3){
            if(m_Conn == 8)
              MBERRCHK(mb->get_entities_by_type(mthis_set, moab::MBHEX, mat_hexes, true), mb);
            else if(m_Conn ==4)
              MBERRCHK(mb->get_entities_by_type(mthis_set, moab::MBTET, mat_hexes, true), mb);
          }
        else if(m_GD == 2){
            if(m_Conn == 4)
              MBERRCHK(mb->get_entities_by_type(mthis_set, moab::MBQUAD, mat_hexes, true), mb);
            else if(m_Conn == 3)
              MBERRCHK(mb->get_entities_by_type(mthis_set, moab::MBTRI, mat_hexes, true), mb);
          }
        
        // tag all the mat_hexes with matid
        std::vector<int> matID(mat_hexes.size(), mset_id);
        MBERRCHK(mb->tag_set_data(MatIDTag, mat_hexes, &matID[0]), mb);
        //
        MBERRCHK(mb->get_adjacencies(mat_hexes, 0, false, mat_nodes, Interface::UNION), mb);
        moab::Range mat_b_nodes = intersect(nodes, mat_nodes);
        
        std::vector<int> bl_node_data(mat_b_nodes.size(), 0);
        std::vector<int> node_tag_data(mat_b_nodes.size(),-1);
        
        // don't do error check, as it is supposed to give error when multiple material case is encountered
        mb->tag_get_data(MNTag, mat_b_nodes, &node_tag_data[0]);
        for(int i=0; i< (int)mat_b_nodes.size(); i++){
            // already a part of some material
            if(node_tag_data[i] != -1){
                bl_node_data[i] = node_tag_data[i]+1;
              }
          }
        
        MBERRCHK(mb->tag_set_data(MNTag, mat_b_nodes, &bl_node_data[0]), mb);
        mat_hexes.clear();
        mat_b_nodes.clear();
        mat_nodes.clear();
        mthis_set = 0;
      } // end handling material set
    
    // if fixmat specified, filter old hex, we don't have to correct both sides of the boundary
    if (fixmat !=-1 && (int) old_hex.size() > 1){
        moab::EntityHandle old_hex_set;
        MBERRCHK(mb->create_meshset(moab::MESHSET_SET, old_hex_set, 1), mb);
        MBERRCHK(mb->add_entities(old_hex_set,&old_hex[0], (int) old_hex.size()), mb);
        // TODO: Find a faster way of doing this
        MBERRCHK(mb->remove_entities(old_hex_set, fixmat_ents), mb);
        old_hex.clear();
        old_hex.empty();
        // the the old hex to be modified
        MBERRCHK(mb->get_entities_by_dimension(old_hex_set, m_GD, old_hex), mb);
      }
    else if(fixmat ==-1 && (int) old_hex.size()>1){
        m_LogFile << "FIXMAT not defined, elements found on either side of specified BL surface, aborting...";
        m_LogFile << "\n\n Define FIXMAT keyword with material id that remains fixed." << std::endl;
        exit(0);
      }
    
    MBERRCHK(mb->get_connectivity(&(*kter), 1, qconn),mb);
    MBERRCHK(mb->get_connectivity(&old_hex[0], 1, old_hex_conn),mb);
    
    // get the normal to the plane of the 2D mesh
    if(m_GD==2)
      get_normal_quad (old_hex_conn, surf_normal);
    if(found == 1 && m_Material !=999999){
        m_LogFile << "Found material set with id " << m_Material << std::endl;
      }
    else if (m_Material !=999999 && found == 0){
        // material set found, but non-existant. Now create this material set
        m_LogFile << "Creating material set with id " << m_Material << std::endl;
        MBERRCHK(mb->create_meshset(moab::MESHSET_SET, mthis_set, 1), mb);
        MBERRCHK(mb->tag_set_data(MTag, &mthis_set, 1, &m_Material), mb);
      }
    else if (m_Material == 999999 && found == 0){
        m_LogFile << "Use old materials for new boundary layer elements. No material specified." << std::endl;
        
        // If no material tag exists in the model create one now
        if((int)m_sets.size() == 0){
            m_LogFile << "Creating material set with id " << m_Material << std::endl;
            MBERRCHK(mb->create_meshset(moab::MESHSET_SET, mthis_set, 1), mb);
            MBERRCHK(mb->tag_set_data(MTag, &mthis_set, 1, &m_Material), mb);
          }
      }
    else{
        m_LogFile << "Unhandled case" << std::endl;
        exit(0);
      }
    
    err = compute_normals();
    if (err!=0){
        m_LogFile << "Failed to compute normals" << std::endl;
        exit(0);
      }
      
    err = push_bulk_mesh(vw);
    if (err!=0){
        m_LogFile << "Failed to push bulk mesh" << std::endl;
        exit(0);
      }
    err = create_bl_elements(vw);
    if (err!=0){
        m_LogFile << "Failed to create boundary layer elements after pushing bulk mesh" << std::endl;
        exit(0);
      }
      
    m_LogFile << "\nTotal Jacobian calls/Min/Max of penultimate hex elements:" << m_JacCalls << ", " << m_JLo << ", " << m_JHi << std::endl;
    
    // convert the final mesh to hex27 for Nek5000 
    if(hex27 == 1){
        moab::Range entities;
        moab::EntityHandle meshset;
        mb->get_entities_by_type(0, MBHEX, entities);
        mb->create_meshset(MESHSET_SET, meshset);
        mb->add_entities(meshset, entities);
        // Add nodes along mid- edge, face and region
        mb->convert_entities(meshset, true, true, true);
        m_LogFile << "Mesh converted from hex 8 to hex 27 elements" << std::endl;
      }
    
    // save the final mesh with boundary layer elements
    MBERRCHK(mb->write_mesh(m_OutFile.c_str()),mb);
    m_LogFile <<  "\n\nWrote Mesh File: " << m_OutFile << std::endl;
    // get the current date and time
    Timer.GetDateTime (szDateTime);
    m_LogFile << "Ending at : " << szDateTime;
    // report/compute the elapsed time
    m_LogFile <<  "Elapsed wall clock time: " << Timer.DiffTime ()
               << " seconds or " << (Timer.DiffTime ())/60.0 << " mins\n";
    m_LogFile <<  "Total CPU time used: " << (double) (clock() - sTime)/CLOCKS_PER_SEC \
               << " seconds" << std::endl;
  }
  
  void PostBL::PrepareIO (int argc, char *argv[], std::string  TestDir)
  // ---------------------------------------------------------------------------
  //! Function: Parser for reading the PostBL specification (.inp) file. \n
  //! Input:    Command line arguments. \n
  //! Output:   none \n
  // ---------------------------------------------------------------------------
  {
    // set and open input output files
    bool bDone = false;
    do{
        if (2 == argc) {
            m_InputFile = argv[1];
            m_LogName = m_InputFile + ".log";
          }
        else if (1 == argc){            
            m_InputFile = TestDir + "/" + (char *)DEFAULT_TEST_POSTBL;
            m_LogName = (std::string)DEFAULT_TEST_POSTBL + ".log";
          }
        
        // open input file for reading
        m_FileInput.open (m_InputFile.c_str(), std::ios::in);
        if (!m_FileInput){
            m_LogFile << "Usage: postbl <filename.inp> " << std::endl;
            m_LogFile << "Default test file can be found here <Meshkit/data>" << std::endl;
            m_LogFile << " Examples input and mesh files are located here <MeshKit/test/algs/postbl_examples>" << std::endl;
            m_FileInput.clear ();
            exit(1);
          }
        else
          bDone = true; // file opened successfully
        
        // open the log file for dumping debug/output statements
        m_LogFile.coss.open (m_LogName.c_str(), std::ios::out);
        if (!m_LogFile.coss){
            m_LogFile <<  "Unable to open file: " << m_LogName << std::endl;
            m_LogFile.coss.clear ();
            exit(1);
          }
        else
          bDone = true; // file opened successfully
        m_LogFile <<  '\n';
        m_LogFile <<  "\t\t---------------------------------------------------------" << '\n';
        m_LogFile <<  "\t\t         Tool to generate Post-mesh Boundary Layers      " << '\n';
        m_LogFile <<  "\t\t\t\tArgonne National Laboratory" << '\n';
        m_LogFile <<  "\t\t\t\t        2012         " << '\n';
        m_LogFile <<  "\t\t---------------------------------------------------------" << '\n';
        m_LogFile <<  "\nsee README file for using the program and details on various cards.\n"<< std::endl;
        
      }while (!bDone);
    
    // Get the meshfile name, surface(s), thickness, intervals and bias
    CParser Parse;
    int maxloop = 1000;
    // count the total number of cylinder commands in each pincellh
    for(;;){
        if(m_nLineNumber > maxloop){
            m_LogFile << "Warning: Didn't find End keyword, stopping after reading 1000 lines from input file" << std::endl;
            break;
          }
        if (!Parse.ReadNextLine (m_FileInput, m_nLineNumber, szInputString,
                                 MAXCHARS, szComment))
          IOErrorHandler (INVALIDINPUT);
        
        // Get tri scheme
        if (szInputString.substr(0,9) == "trischeme"){
            std::istringstream szFormatString (szInputString);
            szFormatString >> m_Card >> tri_sch;
            if(szFormatString.fail())
              IOErrorHandler(INVALIDINPUT);
            m_LogFile <<  m_Card << " name read: "<< tri_sch << std::endl;
          }
        // Get hybrid
        if (szInputString.substr(0,6) == "hybrid"){
            std::istringstream szFormatString (szInputString);
            szFormatString >> m_Card >> hybrid;
            if(szFormatString.fail())
              IOErrorHandler(INVALIDINPUT);
            m_LogFile <<  m_Card << " name read: "<< hybrid << std::endl;
          }
        // Get hybrid
        if (szInputString.substr(0,6) == "fixmat"){
            std::istringstream szFormatString (szInputString);
            szFormatString >> m_Card >> fixmat;
            if(szFormatString.fail())
              IOErrorHandler(INVALIDINPUT);
            m_LogFile <<  m_Card << " name read: "<< fixmat << std::endl;
          }
        // Get MeshFile name
        if (szInputString.substr(0,8) == "meshfile"){
            std::istringstream szFormatString (szInputString);
            szFormatString >> m_Card >> m_MeshFile;
            if(szFormatString.fail())
              IOErrorHandler(INVALIDINPUT);
            m_LogFile <<  m_Card << " name read: "<< m_MeshFile << std::endl;
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
            m_LogFile <<  m_Card << " read: "<< m_SurfId <<std::endl;
          }
        // Get BL surface via neumann set or sideset
        if (szInputString.substr(0,10) == "neumannset"){
            std::istringstream szFormatString (szInputString);
            szFormatString >> m_Card >> m_NeumannSet;
            if(m_NeumannSet < 0 || szFormatString.fail())
              IOErrorHandler(INVALIDINPUT);
            m_LogFile <<  m_Card << " read: "<< m_NeumannSet <<std::endl;
          }
        // Get BL material (block) number
        if (szInputString.substr(0,8) == "material"){
            std::istringstream szFormatString (szInputString);
            szFormatString >> m_Card >> m_Material;
            if(m_Material < 0 || szFormatString.fail())
              IOErrorHandler(INVALIDINPUT);
            m_LogFile <<  m_Card << " read: "<< m_Material <<std::endl;
          }
        
        // Get thickness
        if (szInputString.substr(0,9) == "thickness"){
            std::istringstream szFormatString (szInputString);
            szFormatString >> m_Card >> m_Thickness;
            if(m_Thickness < 0 || szFormatString.fail())
              IOErrorHandler(INVALIDINPUT);
            m_LogFile <<  m_Card << " read: "<< m_Thickness <<std::endl;
          }
        // Get intervals
        if (szInputString.substr(0,9) == "intervals"){
            std::istringstream szFormatString (szInputString);
            szFormatString >> m_Card >> m_Intervals;
            if(m_Intervals < 0 || szFormatString.fail())
              IOErrorHandler(INVALIDINPUT);
            m_LogFile <<  m_Card << " read: "<< m_Intervals <<std::endl;
          }
        // Get bias
        if (szInputString.substr(0,4) == "bias"){
            std::istringstream szFormatString (szInputString);
            szFormatString >> m_Card >> m_Bias;
            if(m_Bias < 0 || szFormatString.fail())
              IOErrorHandler(INVALIDINPUT);
            m_LogFile <<  m_Card << " read: "<< m_Bias <<std::endl;
          }
        // Output file name
        if (szInputString.substr(0,7) == "outfile"){
            std::istringstream szFormatString (szInputString);
            szFormatString >> m_Card >> m_OutFile;
            if(szFormatString.fail())
              IOErrorHandler(INVALIDINPUT);
            m_LogFile <<  m_Card << " name read: "<< m_OutFile <<std::endl;
          }
        // Output file name
        if (szInputString.substr(0,7) == "outfile"){
            std::istringstream szFormatString (szInputString);
            szFormatString >> m_Card >> m_OutFile;
            if(szFormatString.fail())
              IOErrorHandler(INVALIDINPUT);
            m_LogFile <<  m_Card << " name read: "<< m_OutFile <<std::endl;
          }
        // save hex27 mesh
        if (szInputString.substr(0,5) == "hex27"){
            std::istringstream szFormatString (szInputString);
            szFormatString >> m_Card >> hex27;
            if(szFormatString.fail())
              IOErrorHandler(INVALIDINPUT);
            m_LogFile <<  m_Card << " hex27 creation: "<< hex27 <<std::endl;
          }
        if (szInputString.substr(0,3) == "end"){
            break;
          }
      }
  }
  
  void PostBL::IOErrorHandler (ErrorStates ECode) const
  // ---------------------------------------------------------------------------
  //! Function: Displays error messages related to input data \n
  //! Input:    Error code \n
  //! Output:   none \n
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
  
  
} // namespace MeshKit


