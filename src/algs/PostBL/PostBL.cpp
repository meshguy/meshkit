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
  m_LogFile <<  "\n Loading meshfile: " << m_MeshFile << ".." << std::endl;

  // load specified mesh file
  IBERRCHK(imesh->load(0, m_MeshFile.c_str(),0), *imesh);


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
  moab::Tag GDTag, GIDTag, NTag, MTag, STag, FTag, MNTag, MatIDTag, BLNodeIDTag;
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
  moab::Range sets, n_sets, m_sets;
  m_BLDim = m_GD - 1;

  const void* gdim[] = {&m_BLDim};

  MBERRCHK(mb->get_entities_by_type_and_tag(0, moab::MBENTITYSET, &GDTag,
                                            gdim, 1 , sets, moab::Interface::INTERSECT, false), mb);
  MBERRCHK(mb->get_entities_by_type_and_tag(0, moab::MBENTITYSET, &NTag, 0, 1 , n_sets),mb);
  MBERRCHK(mb->get_entities_by_type_and_tag(0, moab::MBENTITYSET, &MTag, 0, 1 , m_sets),mb);

  // Handling NeumannSets (if BL surf in input via NS)
  moab::Range::iterator set_it;
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
  moab::Range quads, nodes,edges, fixmat_ents;
  int dims; // variable to store global id of boundary layer specified in the input file


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
  moab::EntityHandle geom_set;
  MBERRCHK(mb->create_meshset(moab::MESHSET_SET, geom_set, 1), mb);
  MBERRCHK(mb->tag_set_data(GDTag, &geom_set, 1, &m_GD), mb);

  // declare variables before starting BL creation
  std::vector <bool> node_status(false); // size of verts of bl surface
  node_status.resize(nodes.size());
  moab::Range hexes, hex_edge, quad_verts;
  double coords_new_quad[3];
  moab::EntityHandle hex, hex1, hex2;
  int qcount = 0;

  //size of the following is based on element type
  std::vector<moab::EntityHandle> conn, qconn, adj_qconn, tri_conn, tet_conn,
          new_vert(m_Intervals*nodes.size()), old_hex_conn, adj_hexes, adj_quads, adj_hex_nodes1;
  moab::CartVect surf_normal;

  // element on the boundary
  Range::iterator kter = quads.begin();

  std::vector<moab::EntityHandle> old_hex;
  MBERRCHK(mb->get_adjacencies(&(*kter), 1, m_GD, false, old_hex),mb);
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
  moab::Range::iterator mset_it;
  moab::EntityHandle mthis_set;
  int mset_id = 0, found = 0;
  for (mset_it = m_sets.begin(); mset_it != m_sets.end(); mset_it++)  {

      mthis_set = *mset_it;

      // get entity handle of MS specified in the input file
      MBERRCHK(mb->tag_get_data(MTag, &mthis_set, 1, &mset_id), mb);

      // if no material set is specified, we'll have to resolve

      // else just set all the MNTag to 0

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

      //
      std::vector<int> bl_node_data(mat_b_nodes.size(), 0);
      std::vector<int> node_tag_data(mat_b_nodes.size(),-1);
      // don't error check, as it is supposed to give error when multiple material case is encountered
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
      //
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

  // COMPUTE NORMALS
  // get tag data and print
  std::vector<int> all_bl(nodes.size(), 0);
  mb->tag_get_data(MNTag, nodes, &all_bl[0]);
  int count = -1;
  for (Range::iterator kter = nodes.begin(); kter != nodes.end(); ++kter){
      ++count;
      // only one normal in case of single material
      double xdisp = 0.0, ydisp = 0.0, zdisp = 0.0;
      // check if node belongs to one or more materials
      if(all_bl[count] == 0){
          moab::Range adj_for_normal;
          MBERRCHK(mb->get_adjacencies(&(*kter), 1, m_BLDim, false, adj_for_normal, Interface::UNION), mb);
          //create the coordinates of the innermost node corresponding to this node
          moab::Range adj_for_norm = intersect(quads, adj_for_normal);
          // now compute the average normal direction for this vertex
          moab::CartVect rt(0.0, 0.0, 0.0), v(0.0, 0.0, 0.0);

          for(Range::iterator qter = adj_for_norm.begin(); qter != adj_for_norm.end(); ++qter){
              MBERRCHK(mb->get_connectivity(&(*qter), 1, adj_qconn),mb);

              int side_number = 0, sense = 1, offset = 0;
              mb->side_number(old_hex[0], (*qter), side_number, sense, offset);

              if(m_GD==3){
                  get_normal_quad (adj_qconn, v);
                  if(sense == 1){
                      // do nothing
                  }
                  else{
                      v=-v;
                  }
              }
              else if(m_GD==2){
                  if(sense == 1){
                      get_normal_edge(adj_qconn, surf_normal, v);
                  }
                  else{
                      get_normal_edge(adj_qconn, -surf_normal, v);
                  }
              }
              rt = rt + v;
          }
          if(rt.length() !=0){
              xdisp=rt[0]/rt.length();
              ydisp=rt[1]/rt.length();
              zdisp=rt[2]/rt.length();
          }
          else{
              xdisp=0.0;
              ydisp=0.0;
              zdisp=0.0;
          }
      }
      else if(all_bl[count] > 0 && fixmat != -1){ // node belongs to more than one material and fixmat specified
          // NODE B/W MATERIALS
          moab::Range adj_for_normal;
          MBERRCHK(mb->get_adjacencies(&(*kter), 1, m_BLDim, false, adj_for_normal, Interface::UNION), mb);
          //create the coordinates of the innermost node corresponding to this node
          moab::Range adj_for_norm = intersect(quads, adj_for_normal);
          // now compute the average normal direction for this vertex
          moab::CartVect rt(0.0, 0.0, 0.0), v(0.0, 0.0, 0.0);

          for(Range::iterator qter = adj_for_norm.begin(); qter != adj_for_norm.end(); ++qter){
              MBERRCHK(mb->get_connectivity(&(*qter), 1, adj_qconn),mb);

              moab::Range this_quad_hex;
              MBERRCHK(mb->get_adjacencies(&(*qter), 1, m_GD, false, this_quad_hex, moab::Interface::UNION),mb);
              moab::Range quad_hex = intersect(fixmat_ents, this_quad_hex);

              int side_number = 0, sense = 1, offset = 0;
              Range::iterator hexter = quad_hex.begin();
              mb->side_number(*hexter, (*qter), side_number, sense, offset);

              if(m_GD==3){
                  get_normal_quad (adj_qconn, v);
                  if(sense == 1 ){
                      v=-v;
                  }
                  else{
                      // do nothing
                  }
              }
              else if(m_GD==2){
                  if(sense == 1 ){
                      get_normal_edge(adj_qconn, -surf_normal, v);
                  }
                  else{
                      get_normal_edge(adj_qconn, surf_normal, v);
                  }
              }
              rt = rt + v;
          }
          if(rt.length() !=0){
              xdisp=rt[0]/rt.length();
              ydisp=rt[1]/rt.length();
              zdisp=rt[2]/rt.length();
          }
          else{
              xdisp=0.0;
              ydisp=0.0;
              zdisp=0.0;
          }
      }
      else if(all_bl[count] > 0 && fixmat == -1){ // node belongs to more than one material and fixmat not specified
          // NODE ON BOUNDARY
          // get the edge that is not on the boundary and count how many such edges we have
          moab::Range adj_for_normal;
          int nEdgeDim = 1;
          MBERRCHK(mb->get_adjacencies(&(*kter), 1, nEdgeDim, true, adj_for_normal, Interface::UNION), mb);
          moab::Range edge_normal;
          if(m_GD == 2)
              edge_normal = subtract(adj_for_normal, quads);
          else if(m_GD == 3)
              edge_normal = subtract(adj_for_normal, edges);


          if(edge_normal.size() > 1){
              double ncoord[3];
              MBERRCHK(mb->get_coords(&(*kter), 1, ncoord),mb);

              m_LogFile << "Work in progress !! :- MULTIPLE NORMALS ARE NEEDED AT" << ncoord[0]
                        << ", " << ncoord[1] << ", " << ncoord[2] << " #normals " << edge_normal.size() << std::endl;
              exit(0);
          }
          else{
              m_LogFile << "We've one edge seperating materials 1 NORMAL IS NEEDED" << edge_normal.size() << std::endl;
              moab::Range edge_conn;
              MBERRCHK(mb->get_connectivity(&(*edge_normal.begin()), 1, edge_conn),mb);
              // now get the normal direction for this edge
              moab::CartVect coords[2], bl_coord[1];
              moab::CartVect AB;
              MBERRCHK(mb->get_coords(&(*kter), 1, (double*) &bl_coord[0]), mb);
              MBERRCHK(mb->get_coords(edge_conn, (double*) &coords[0]), mb);
              for(int d = 0; d<2; d++){
                  if(bl_coord[0][0] == coords[d][0] &&
                          bl_coord[0][1] == coords[d][1] &&
                          bl_coord[0][2] == coords[d][2]){
                      // do nothing
                  }
                  else{
                      AB = (bl_coord[0] - coords[d]);
                  }
              }

              xdisp=AB[0]/AB.length();
              ydisp=AB[1]/AB.length();
              zdisp=AB[2]/AB.length();
          }
      }
      else if(all_bl[count] < 0 ){
          m_LogFile << "Material must have associated with BLNode: Error, shouldn't have gotten here: " << count << std::endl;
      }


      // after the normal computation is done create new BL nodes
      double coords_bl_quad[3], move = 0.0;
      MBERRCHK(mb->get_coords(&(*kter), 1, coords_bl_quad),mb);

      double temp = 0;
      double num = m_Thickness*(m_Bias-1)*(pow(m_Bias, m_Intervals -1));
      double deno = pow(m_Bias, m_Intervals) - 1;
      if (deno !=0)
          temp = num/deno;
      else
          temp = m_Thickness/m_Intervals;

      // loop thru intervals to create BL nodes
      for(int j=0; j< m_Intervals; j++){

          move+= temp/pow(m_Bias,j);
          if (debug){
              m_LogFile <<  " move:" << move;
          }
          // now compute the coords of the new vertex
          coords_new_quad[0] = coords_bl_quad[0]-move*xdisp;
          coords_new_quad[1] = coords_bl_quad[1]-move*ydisp;
          coords_new_quad[2] = coords_bl_quad[2]-move*zdisp;

          int nid = count*m_Intervals+j;
          // Possible TODO's
          //TODO: Check if this vertex is possible (detect possible collision with geometry)
          // TODO: See posibility of using ray tracing
          // TODO: Parallize: Profile T-junction model and try to device an algorithm
          // TODO: Modularize node creation part and use doxygen for all code and design of code, python design and test cases - current functions in code:
          // Setup this, Execute this -- break info sub functions and classes,
          // prepareIO --make this optional when using python,
          // get normal (2d and 3d) -- can be combined to one function
          // get det jacobian (hex elements) --needs check for other elements
          //
          MBERRCHK(mb->create_vertex(coords_new_quad, new_vert[nid]), mb);
          if (debug){
              m_LogFile << std::setprecision (3) << std::scientific << " : created node:" << (nid + 1)
                        << " of " << new_vert.size() << " new nodes:- coords: " << coords_new_quad[0]
                        << ", " << coords_new_quad[1] << ", " << coords_new_quad[2]  << std::endl;
          }
      }
  }

  // swap nodal coordinates of input nodes with innermost BL nodes to push the bulk mesh
  count = -1;
  for (Range::iterator kter = nodes.begin(); kter != nodes.end(); ++kter){
      ++count;
      if(all_bl[count] == 0 ){
          double coords_new_quad[3];
          double coords_old_quad[3];

          int nid = (count+1)*m_Intervals - 1;
          MBERRCHK(mb->get_coords(&new_vert[nid], 1, coords_new_quad),mb);

          MBERRCHK(mb->get_coords(&(*kter), 1, coords_old_quad),mb);

          MBERRCHK(mb->set_coords(&(*kter), 1, coords_new_quad),mb);

          MBERRCHK(mb->set_coords(&new_vert[nid], 1, coords_old_quad),mb);
          if (debug) {
              m_LogFile << std::setprecision (3) << std::scientific << " : NID:" << (nid)
                        << coords_old_quad[0]
                        << ", " << coords_old_quad[1] << ", " << coords_old_quad[2] << " OLD:- coords: NEW" << coords_new_quad[0]
                        << ", " << coords_new_quad[1] << ", " << coords_new_quad[2]  << std::endl;
          }
      }
      else if(all_bl[count] > 0 && fixmat != -1){ // node belongs to more than one material and fixmat specified
          // NODE B/W MATERIALS
          double coords_new_quad[3];
          double coords_old_quad[3];

          int nid = (count+1)*m_Intervals - 1;
          MBERRCHK(mb->get_coords(&new_vert[nid], 1, coords_new_quad),mb);

          MBERRCHK(mb->get_coords(&(*kter), 1, coords_old_quad),mb);

          MBERRCHK(mb->set_coords(&(*kter), 1, coords_new_quad),mb);

          MBERRCHK(mb->set_coords(&new_vert[nid], 1, coords_old_quad),mb);

          m_LogFile << std::setprecision (3) << std::scientific << "FM : NID:" << (nid)
                    << coords_old_quad[0]
                    << ", " << coords_old_quad[1] << ", " << coords_old_quad[2] << " OLD:- coords: NEW" << coords_new_quad[0]
                    << ", " << coords_new_quad[1] << ", " << coords_new_quad[2]  << std::endl;

          // now find the hex in fixmat_ents that was affected by this swap and UNswap the coords
          moab::Range fhex;
          MBERRCHK(mb->get_adjacencies(&(*kter), 1, m_GD, false, fhex, Interface::UNION), mb);
          moab::Range fmhex= intersect(fhex,fixmat_ents);
          std::vector<moab::EntityHandle> fmconn;
          for(Range::iterator fmter = fmhex.begin(); fmter != fmhex.end(); ++fmter){
              MBERRCHK(mb->get_connectivity(&(*fmter), 1, fmconn),mb);
              for(int i=0; i < (int)fmconn.size(); i++){
                  if((*kter) == fmconn[i]){
                      //we have a node to unswap or reset connectivity for fixmat
                      fmconn[i] = new_vert[nid];
                  }
              }
              MBERRCHK(mb->set_connectivity(*fmter, &fmconn[0], fmconn.size()), mb);
              double jac = 0;
              vw.quality_measure(*fmter, MB_JACOBIAN, jac);
              ++m_JacCalls;
              if (jac < 0){
                  m_LogFile << "Negative Jacobian, check BL thickness/intervals. Stopping." << std::endl;
                  exit(0);
                }
              //                MBERRCHK(mb->add_adjacencies((*fmter), hexes, true), mb);
              //                MBERRCHK(mb->add_adjacencies((*fmter), fixmat_ents, true), mb);

          }

          m_LogFile << " We're here in MM case, now go along the edge --- have fun in the process !!" << std::endl;
          // material boundary - get edge direction and length

          //  find_min_edge_length(adj_qconn_r, qconn[i], nodes, m_MinEdgeLength);

          // check to see if this is a fixmat case
      }
      else if(all_bl[count] > 0 && fixmat == -1){ // node belongs to more than one material and fixmat not specified
          m_LogFile << " I'm here" << std::endl;
          // NODE ON BOUNDARY
          double coords_new_quad[3];
          double coords_old_quad[3];

          int nid = (count+1)*m_Intervals - 1;
          MBERRCHK(mb->get_coords(&new_vert[nid], 1, coords_new_quad),mb);

          MBERRCHK(mb->get_coords(&(*kter), 1, coords_old_quad),mb);

          MBERRCHK(mb->set_coords(&(*kter), 1, coords_new_quad),mb);

          MBERRCHK(mb->set_coords(&new_vert[nid], 1, coords_old_quad),mb);

          m_LogFile << std::setprecision (3) << std::scientific << " : NID:" << (nid)
                    << coords_old_quad[0]
                    << ", " << coords_old_quad[1] << ", " << coords_old_quad[2] << " OLD:- coords: NEW" << coords_new_quad[0]
                    << ", " << coords_new_quad[1] << ", " << coords_new_quad[2]  << std::endl;
      }
  }

  // check for the volume of penultimate elements if -ve volume encountered. Report.
  count = -1;
  // Try to move another layer attached to it.
  for (Range::iterator kter = nodes.begin(); kter != nodes.end(); ++kter){
      ++count;
      for(int j=0; j< m_Intervals; j++){
          int nid = count*m_Intervals+j;
          double coords_new_quad[3];
          MBERRCHK(mb->get_coords(&new_vert[nid], 1, coords_new_quad),mb);
          m_LogFile << std::setprecision (3) << std::scientific << " : NID:" << (nid)
                    << " of " << new_vert.size() << " new nodes:- coords: " << coords_new_quad[0]
                    << ", " << coords_new_quad[1] << ", " << coords_new_quad[2]  << std::endl;
      }
  }
  // shoot multiple normals for multiple materials case only.

  // This can be used of regular case also how to invoke it. mention in algorithm

  // mention local and global smoothing.

  //    qcount = -1;
  //    int flag[quads.size()];

  // Now start creating New elements
  for (Range::iterator kter = quads.begin(); kter != quads.end(); ++kter){

      qcount++;

      std::vector<moab::EntityHandle> old_hex;
      MBERRCHK(mb->get_adjacencies(&(*kter), 1, m_GD, false, old_hex),mb);
      double jac = 0;
      vw.quality_measure(old_hex[0], MB_JACOBIAN, jac);
      if (qcount == 1){
          m_JLo = jac;
        }
      ++m_JacCalls;

      if(m_JHi < jac)
        m_JHi = jac;

      if(m_JLo > jac)
        m_JLo = jac;

      if (jac < 0){
          m_LogFile << "Negative Jacobian, check BL thickness/intervals. Stopping." << std::endl;
          exit(0);
        }

      MBERRCHK(mb->get_connectivity(&(*kter), 1, qconn),mb);
      double one_node_in_quad[3];
      for (int i=0; i<m_BElemNodes; i++){

          int node_tag_id = 0;
          MBERRCHK(mb->tag_get_data(BLNodeIDTag, &qconn[i], 1, &node_tag_id) ,mb);
          MBERRCHK(mb->get_coords(&qconn[i], 1, one_node_in_quad),mb);
          m_LogFile << std::setprecision (3) << std::scientific << " new nodes:- coords: " << one_node_in_quad[0]
                    << ", " << one_node_in_quad[1] << ", " << one_node_in_quad[2]  << std::endl;

          //populate the connectivity after creating nodes for this BL node
          for(int j=0; j< m_Intervals; j++){
              if(m_Conn == 8 && m_BElemNodes == 4){ // hex
                  int nid = node_tag_id*m_Intervals + j;
                  if(m_Intervals == 1){
                      conn[m_Conn*j +i] = qconn[i];
                      conn[m_Conn*j + i+m_BElemNodes] = new_vert[nid];
                  }
                  else if(j==0){
                      conn[m_Conn*j +i] = qconn[i];
                      conn[m_Conn*j + i+m_BElemNodes] = new_vert[nid + m_Intervals - 2];
                  }
                  else if(j==(m_Intervals-1)){
                      conn[m_Conn*j +i] = new_vert[nid - m_Intervals + 1];
                      conn[m_Conn*j + i+m_BElemNodes] = new_vert[nid];
                  }
                  else {
                      conn[m_Conn*j +i] = new_vert[nid + m_Intervals - 2*j -1];
                      conn[m_Conn*j + i+m_BElemNodes] = new_vert[nid + m_Intervals - 2*j -2];
                  }
              }
              else if(m_Conn == 4 && m_BElemNodes == 2){ // Quads
                  int nid = node_tag_id*m_Intervals+j;
                  if(m_Intervals == 1){
                      conn[m_Conn*j +i] = new_vert[nid];
                      conn[m_Conn*j + i+m_BElemNodes] = qconn[m_BElemNodes-i-1];
                  }
                  else if(j==0){
                      conn[m_Conn*j +i] = qconn[m_BElemNodes-i-1];
                      conn[m_Conn*j + i+m_BElemNodes] = new_vert[nid + m_Intervals - 2];
                  }
                  else if(j==(m_Intervals-1)){
                      conn[m_Conn*j +i] = new_vert[nid];
                      conn[m_Conn*j + m_BElemNodes + 1 -i] = new_vert[nid - m_Intervals + 1];
                  }
                  else {
                      conn[m_Conn*j +i] = new_vert[nid + m_Intervals - 2*j -2];
                      conn[m_Conn*j + m_BElemNodes + 1 -i] = new_vert[nid + m_Intervals - 2*j -1];
                  }
              }
              else if(m_Conn == 4 && m_BElemNodes == 3){ // && hybrid == true make wedges aka prisms for tet mesh
                  int nid = node_tag_id*m_Intervals+j;
                  if(m_Intervals == 1){
                      conn[m_Conn*j +i] = qconn[i];
                      conn[m_Conn*j + i+m_BElemNodes] = new_vert[nid];
                  }
                  else if(j==0){
                      conn[m_HConn*j +i] = qconn[i];
                      conn[m_HConn*j + i+m_BElemNodes] = new_vert[nid + m_Intervals - 2];
                  }
                  else if(j==(m_Intervals-1)){
                      conn[m_HConn*j +i] = new_vert[nid - m_Intervals + 1];
                      conn[m_HConn*j + i+m_BElemNodes] = new_vert[nid];
                  }
                  else {
                      conn[m_HConn*j +i] = new_vert[nid + m_Intervals - 2*j -1];
                      conn[m_HConn*j + i+m_BElemNodes] = new_vert[nid + m_Intervals - 2*j -2];
                  }
              }
              else if(m_Conn == 3 && m_BElemNodes == 2){ // make quads for tri mesh
                  int nid = node_tag_id*m_Intervals+j;
                  if(m_Intervals == 1){
                      conn[m_Conn*j +i] = new_vert[nid];
                      conn[m_Conn*j + i+m_BElemNodes] = qconn[m_BElemNodes-i-1];
                  }
                  else if(j==0){
                      conn[m_HConn*j +i] = qconn[m_BElemNodes-i-1];
                      conn[m_HConn*j + i+m_BElemNodes] = new_vert[nid + m_Intervals - 2];
                  }
                  else if(j==(m_Intervals-1)){
                      conn[m_HConn*j +i] = new_vert[nid];
                      conn[m_HConn*j + m_BElemNodes + 1 - i] = new_vert[nid - m_Intervals + 1];
                  }
                  else {
                      conn[m_HConn*j +i] = new_vert[nid + m_Intervals - 2*j -2];
                      conn[m_HConn*j + m_BElemNodes + 1 - i] = new_vert[nid + m_Intervals - 2*j -1];
                  }
              }
              else {
                  std::cout << "ERROR: cannot create BL elements: element type not supported." << std::endl;
                  exit(0);
              }
          }
      }

      //TODO: Set Connectivity of tet's, break prisms into 3 tets, Another loop is required.

      if(m_Conn == 4 && m_BElemNodes == 3 && hybrid == false){
          for(int c=0; c<m_Intervals; c++){
                  // first tet
                tet_conn[m_Conn*c*3] =     conn[c*m_HConn + 0];
                tet_conn[m_Conn*c*3 + 1] = conn[c*m_HConn + 4];
                tet_conn[m_Conn*c*3 + 2] = conn[c*m_HConn + 5];
                tet_conn[m_Conn*c*3 + 3] = conn[c*m_HConn + 3];
                // middle tet
                tet_conn[m_Conn*c*3 + 4] = conn[c*m_HConn + 0];
                tet_conn[m_Conn*c*3 + 5] = conn[c*m_HConn + 1];
                tet_conn[m_Conn*c*3 + 6] = conn[c*m_HConn + 2];
                tet_conn[m_Conn*c*3 + 7] = conn[c*m_HConn + 5];
                //last tet
                tet_conn[m_Conn*c*3 + 8] = conn[c*m_HConn + 0];
                tet_conn[m_Conn*c*3 + 9] = conn[c*m_HConn + 1];
                tet_conn[m_Conn*c*3 + 10] = conn[c*m_HConn + 5];
                tet_conn[m_Conn*c*3 + 11] = conn[c*m_HConn + 4];
          }
      }



      if(m_Conn == 3 && m_BElemNodes == 2 && hybrid == false){
          for(int c=0; c<m_Intervals; c++){
              if(tri_sch == 1){
                  // lower triangle
                  tri_conn[m_Conn*c*2] =     conn[c*m_HConn + 0];
                  tri_conn[m_Conn*c*2 + 1] = conn[c*m_HConn + 1];
                  tri_conn[m_Conn*c*2 + 2] = conn[c*m_HConn + 2];
                  // upper triangl
                  tri_conn[m_Conn*c*2 + 3] = conn[c*m_HConn + 2];
                  tri_conn[m_Conn*c*2 + 4] = conn[c*m_HConn + 3];
                  tri_conn[m_Conn*c*2 + 5] = conn[c*m_HConn + 0];
              }
              else if(tri_sch == 2){
                  // lower triangle
                  tri_conn[m_Conn*c*2] =     conn[c*m_HConn + 0];
                  tri_conn[m_Conn*c*2 + 1] = conn[c*m_HConn + 1];
                  tri_conn[m_Conn*c*2 + 2] = conn[c*m_HConn + 2];
                  // upper triangle
                  tri_conn[m_Conn*c*2 + 3] = conn[c*m_HConn + 2];
                  tri_conn[m_Conn*c*2 + 4] = conn[c*m_HConn + 3];
                  tri_conn[m_Conn*c*2 + 5] = conn[c*m_HConn + 0];
              }
          }
      }

      // create boundary layer hexes
      for(int j=0; j< m_Intervals; j++){
          if(m_Conn == 8){
              MBERRCHK(mb->create_element(moab::MBHEX, &conn[j*m_Conn], m_Conn, hex),mb);
          }
          else if(m_Conn==4 && m_GD ==3 && hybrid == true){
              MBERRCHK(mb->create_element(MBPRISM, &conn[j*6], 6, hex),mb);
          }
          else if(m_Conn==4 && m_GD ==2){
              MBERRCHK(mb->create_element(MBQUAD, &conn[j*m_Conn], m_Conn, hex),mb);
          }
          else if(m_Conn==3 && m_GD ==2 && hybrid == true){
              MBERRCHK(mb->create_element(MBQUAD, &conn[j*m_HConn], m_HConn, hex),mb);
          }
          else if(m_Conn==3 && m_GD ==2 && hybrid == false){
              MBERRCHK(mb->create_element(MBTRI, &tri_conn[j*m_Conn*2], m_Conn, hex),mb);
              MBERRCHK(mb->create_element(MBTRI, &tri_conn[j*m_Conn*2+3], m_Conn, hex1),mb);
              //       MBERRCHK(mb->add_entities(mthis_set, &hex1, 1), mb);
              //        MBERRCHK(mb->add_entities(smooth_set, &hex1, 1), mb);
          }
          else if(m_Conn==4 && m_GD ==3 && hybrid == false){
              MBERRCHK(mb->create_element(MBTET, &tet_conn[j*m_Conn*3], m_Conn, hex),mb);
              MBERRCHK(mb->create_element(MBTET, &tet_conn[j*m_Conn*3+4], m_Conn, hex1),mb);
              MBERRCHK(mb->create_element(MBTET, &tet_conn[j*m_Conn*3+8], m_Conn, hex2),mb);
              //       MBERRCHK(mb->add_entities(mthis_set, &hex1, 1), mb);
              //        MBERRCHK(mb->add_entities(smooth_set, &hex1, 1), mb);
          }
          // add this hex to a block
          moab::Range adj_hex_for_mat;
          // int hmat_id = 0;

          if(mthis_set == 0){
              //               std::vector<int> hmat_id(qconn.size(), 0);
              //                MBERRCHK(mb->tag_get_data(MNTag, &qconn[0], 1, &hmat_id) ,mb);
              //        MBERRCHK(mb->get_adjacencies(&qconn[0], 1, m_GD, false, adj_hex_for_mat, moab::Interface::INTERSECT), mb);
              MBERRCHK(mb->get_adjacencies(&(*kter), 1, m_GD, false, adj_hex_for_mat, moab::Interface::INTERSECT), mb);
              MBERRCHK(mb->add_adjacencies(hex, adj_hex_for_mat, true), mb);

              std::vector<int> hmat_id(adj_hex_for_mat.size(), 0);

              // this will lead to an error, so no error checking, new adj hexes don't have matidtag
              mb->tag_get_data(MatIDTag, adj_hex_for_mat, &hmat_id[0]);//, mb);
              for(int p=0; p< (int)hmat_id.size(); p++){
                  if(hmat_id[p] !=0){
                      // this is our mat id for this hex
                      moab::EntityHandle mat_set = 0;
                      for (set_it = m_sets.begin(); set_it != m_sets.end(); set_it++)  {
                          mat_set = *set_it;
                          int set_id;
                          MBERRCHK(mb->tag_get_data(MTag, &mat_set, 1, &set_id), mb);
                          if(set_id == hmat_id[p])
                              break;
                      }
                      MBERRCHK(mb->add_entities(mat_set, &hex, 1), mb);
                      if(m_Conn==3 && m_GD ==2 && hybrid == false)
                          MBERRCHK(mb->add_entities(mthis_set, &hex1, 1), mb);
                      if(m_Conn==4 && m_GD ==3 && hybrid == false){
                          MBERRCHK(mb->add_entities(mthis_set, &hex1, 1), mb);
                          MBERRCHK(mb->add_entities(mthis_set, &hex2, 1), mb);
                      }
                  }
              }


              //    here find the material set for this hex

              //                m_LogFile << "Find out the material that this hex corresponds to?? bailing out" << std::endl;
              //                exit(0);
          }
          else {
              MBERRCHK(mb->add_entities(mthis_set, &hex, 1), mb);
              if(m_Conn==3 && m_GD ==2 && hybrid == false)
                  MBERRCHK(mb->add_entities(mthis_set, &hex1, 1), mb);
              if(m_Conn==4 && m_GD ==3 && hybrid == false){
                  MBERRCHK(mb->add_entities(mthis_set, &hex1, 1), mb);
                  MBERRCHK(mb->add_entities(mthis_set, &hex2, 1), mb);
              }

          }
          // mark entities for smoothing
          //                MBERRCHK(mb->add_entities(smooth_set, &hex, 1), mb);
          // add geom dim tag
          MBERRCHK(mb->add_entities(geom_set, &hex, 1), mb);

          //              // TODO: Add Local Smooting
      }

  }

  //    all_elems.clear();
  //    moab::Range skin_verts;
  //    MBERRCHK(mb->get_entities_by_dimension(0, 3, all_elems,true),mb);

  //    moab::Skinner skinner(mb);
  //    skinner.find_skin(0, all_elems, 0, skin_verts);

  //    m_LogFile << "setting 'fixed'' tag = 1 on verts in the skin = " <<  skin_verts.size() << std::endl;

  //    // set fixed tag = 1 on all skin verts
  //    std::vector<int> all_skin_data(skin_verts.size(), 1);
  //    MBERRCHK(mb->tag_set_data(FTag, skin_verts, &all_skin_data[0]), mb);

  m_LogFile << "\nTotal Jacobian calls/Min/Max of penultimate hex elements:" << m_JacCalls << ", " << m_JLo << ", " << m_JHi << std::endl;

  // save the final boundary layer mesh
  MBERRCHK(mb->write_mesh(m_OutFile.c_str()),mb);
  m_LogFile <<  "\n\nWrote Mesh File: " << m_OutFile << std::endl;
  // get the current date and time
  Timer.GetDateTime (szDateTime);
  m_LogFile << "Ending at : " << szDateTime;
  // report/compute the elapsed time
  m_LogFile <<  "Elapsed wall clock time: " << Timer.DiffTime ()
             << " seconds or " << (Timer.DiffTime ())/60.0 << " mins\n";
  m_LogFile <<  "AL2 Total CPU time used: " << (double) (clock() - sTime)/CLOCKS_PER_SEC \
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

    // count the total number of cylinder commands in each pincellh
    for(;;){
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
        // Debug flag
        if (szInputString.substr(0,5) == "debug"){
            std::istringstream szFormatString (szInputString);
            szFormatString >> m_Card >> debug;
            if(szFormatString.fail())
                IOErrorHandler(INVALIDINPUT);
            m_LogFile <<  m_Card << " name read: "<< debug <<std::endl;
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

void PostBL::get_normal_quad (std::vector<EntityHandle>conn, moab::CartVect &v)
// ---------------------------------------------------------------------------
//! Function: Get normal of a quad \n
//! Input:    conn \n
//! Output:   moab::CartVect v \n
// ---------------------------------------------------------------------------
{
    moab::CartVect coords[3];
    MBERRCHK(mb->get_coords(&conn[0], 3, (double*) &coords[0]), mb);
    moab::CartVect AB(coords[1] - coords[0]);
    moab::CartVect BC(coords[2] - coords[1]);
    moab::CartVect normal = AB*BC;
    normal = normal/normal.length();
    v = normal;
}

void PostBL::get_normal_edge (std::vector<EntityHandle>conn, moab::CartVect BC, moab::CartVect &v)
// ---------------------------------------------------------------------------
//! Function: Get normal of a edge along its quad \n
//! Input:    conn of edge, normal to the surf \n
//! Output:   moab::CartVect v \n
// ---------------------------------------------------------------------------
{
    moab::CartVect coords[2];
    MBERRCHK(mb->get_coords(&conn[0], 2, (double*) &coords[0]), mb);
    moab::CartVect AB(coords[1] - coords[0]);
    moab::CartVect normal = AB*BC;
    normal = normal/normal.length();
    v = normal;
}

void PostBL::find_min_edge_length(moab::Range adj_qn, moab::EntityHandle node , moab::Range bl_nodes, double &e_len)
// ---------------------------------------------------------------------------
//! Function: Get minimum edge length from several adjacent quads/edges specified \n
//! Input:    verts of quads, BL vert \n
//! Output:   distance b/w BL vert and inner vert \n
// ---------------------------------------------------------------------------
{
    // get nodes adj(a) to BL node
    moab::CartVect coords[1];
    double len = 0; e_len = 0;
    MBERRCHK(mb->get_coords(&node, 1, (double*) &coords[0]), mb);
    moab::Range non_bl;
    non_bl = subtract(adj_qn, bl_nodes);
    int n_non_bl = (int) non_bl.size();

    // if there are no bl nodes - this case has already been dealt with
    if(non_bl.size() > 0){
        moab::CartVect non_bl_coords[4];
        for(int i=0; i< n_non_bl; i++){
            MBERRCHK(mb->get_coords(non_bl,(double*) &non_bl_coords[0]), mb);
            moab::CartVect edge(coords[0] - non_bl_coords[0]);
            len = edge.length();
            if(i==0)
                e_len = len;
            if(len < e_len)
                e_len = len;
        }
    }
    // m_LogFile << " node minimum edge length" << e_len << std::endl;
}

void PostBL::get_det_jacobian(std::vector<moab::EntityHandle> conn, int offset, double &AvgJ)
// ---------------------------------------------------------------------------
//! Function: Get determinant of jacobian \n
//! Input:    conn \n
//! Output:   vector x, y and z \n
// ---------------------------------------------------------------------------
{
    //TODO: Add quality check for tri/quad and pyramids
    if(m_Conn ==8){
        ++m_JacCalls;
        moab::CartVect vertex[8], xi;
        mstream m_LogFile;
        MBERRCHK(mb->get_coords(&conn[offset], 8, (double*) &vertex[0]), mb);

        double corner[8][3] = { { -1, -1, -1 },
                                {  1, -1, -1 },
                                {  1,  1, -1 },
                                { -1,  1, -1 },
                                { -1, -1,  1 },
                                {  1, -1,  1 },
                                {  1,  1,  1 },
                                { -1,  1,  1 } };

        for (unsigned j = 0; j < 8; ++j) {
            xi[0] = corner[j][0];
            xi[1] = corner[j][1];
            xi[2] = corner[j][2];
            Matrix3 J(0.0);
            double detJ = 0;
            for (unsigned i = 0; i < 8; ++i) {
                const double   xi_p = 1 + xi[0]*corner[i][0];
                const double  eta_p = 1 + xi[1]*corner[i][1];
                const double zeta_p = 1 + xi[2]*corner[i][2];
                const double dNi_dxi   = corner[i][0] * eta_p * zeta_p;
                const double dNi_deta  = corner[i][1] *  xi_p * zeta_p;
                const double dNi_dzeta = corner[i][2] *  xi_p *  eta_p;
                J(0,0) += dNi_dxi   * vertex[i][0];
                J(1,0) += dNi_dxi   * vertex[i][1];
                J(2,0) += dNi_dxi   * vertex[i][2];
                J(0,1) += dNi_deta  * vertex[i][0];
                J(1,1) += dNi_deta  * vertex[i][1];
                J(2,1) += dNi_deta  * vertex[i][2];
                J(0,2) += dNi_dzeta * vertex[i][0];
                J(1,2) += dNi_dzeta * vertex[i][1];
                J(2,2) += dNi_dzeta * vertex[i][2];
            }
            J *= 0.125;
            detJ = J.determinant();
            if(detJ <= 0.0){
                m_LogFile << "We've negative jacobian at the hex corner: "<< j+1 << std::endl;
                exit(0);
            }
            AvgJ+=detJ;
        }
        AvgJ/=8;
        if(m_JacCalls == 1){
            m_JLo = AvgJ;
            m_JHi = AvgJ;
        }
        else if(AvgJ < m_JLo){
            m_JLo = AvgJ;
        }
        else if(AvgJ > m_JHi){
            m_JHi = AvgJ;
        }
    }
}
} // namespace MeshKit


