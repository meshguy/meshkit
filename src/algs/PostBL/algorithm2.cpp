/*********************************************
 July 31,2013
PostBL: Post-Mesh BOundary Layer Improved Algorithms 2
Argonne National Laboratory
*********************************************/

#include "../meshkit/PostBL.hpp"
using namespace MeshKit;
#include <string.h>

void PostBL:: Algo2(){

    m_LogFile << "\nIn execute this : creating boundary layer elements.." <<  std::endl;
    // start the timer
    CClock Timer;
    clock_t sTime = clock();
    std::string szDateTime;
    Timer.GetDateTime (szDateTime);

    m_LogFile <<  "\nStarting out at : " << szDateTime << std::endl;
    m_LogFile <<  "\n Loading meshfile: " << m_MeshFile << ".." << std::endl;

    // load specified mesh file
    IBERRCHK(imesh->load(0, m_MeshFile.c_str(),0), *imesh);
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
    moab::Range quads, nodes, fixmat_ents;
    int dims; // variable to store global id of boundary layer specified in the input file

    // Method 1: INPUT by NeumannSet
    if(m_NeumannSet != -1 && this_set != 0){
        MBERRCHK(mb->get_entities_by_dimension(this_set, m_BLDim, quads,true),mb);
        if (quads.size() <=0){
            m_LogFile <<  " No quads found, aborting.. " << std::endl;
            exit(0);
        }

        MBERRCHK(mb->get_adjacencies(quads, 0, false, nodes, moab::Interface::UNION),mb);
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
    moab::Range edges, hexes, hex_edge, quad_verts;
    double coords_new_quad[3];
    moab::EntityHandle hex, hex1;
    int qcount = 0;

    //size of the following is based on element type
    std::vector<moab::EntityHandle> conn, qconn, adj_qconn, tri_conn,
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
    if(mb->type_from_handle(old_hex[0]) == MBHEX){
        m_Conn = 8;
        m_BElemNodes = 4;
        m_HConn = 8;
        //allocating based on element type
        conn.resize(m_Intervals*m_Conn), qconn.resize(m_BElemNodes), adj_qconn.resize(m_BElemNodes),
                old_hex_conn.resize(m_Conn), adj_hex_nodes1.resize(m_Conn);
    }
    else if(mb->type_from_handle(old_hex[0]) == MBTET){
        m_Conn = 4;
        m_HConn = 6;
        m_BElemNodes = 3;
        //allocating based on element type - thrice the number of elements
        if(hybrid)
            conn.resize(m_Intervals*6);
        else
            conn.resize(2*m_Intervals*m_Conn);
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
        for(int i=0; i<mat_b_nodes.size(); i++){
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
    if (fixmat !=0 && (int) old_hex.size() > 1){
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
    else if(fixmat ==0 && old_hex.size()>1){
        m_LogFile << "FIXMAT not defined, elements found on either side of specified BL surface, aborting...";
        exit(0);
    }

    MBERRCHK(mb->get_connectivity(&(*kter), 1, qconn),mb);
    MBERRCHK(mb->get_connectivity(&old_hex[0], 1, old_hex_conn),mb);
    // get the normal to the plane of the 2D mesh
    if(m_GD==2)
        get_normal_quad (old_hex_conn, surf_normal);
    if(found == 1 && m_Material !=999){
        m_LogFile << "Found material set with id " << m_Material << std::endl;
    }
    else if (m_Material !=999 && found == 0){
        // material set found, but non-existant. Now create this material set
        m_LogFile << "Creating material set with id " << m_Material << std::endl;
        MBERRCHK(mb->create_meshset(moab::MESHSET_SET, mthis_set, 1), mb);
        MBERRCHK(mb->tag_set_data(MTag, &mthis_set, 1, &m_Material), mb);
    }
    else if (m_Material == 999 && found == 0){
        //
        m_LogFile << "Use old materials for new boundary layer elements. No material specified." << std::endl;
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
                MBERRCHK(mb->side_number(old_hex[0], (*qter), side_number, sense, offset), mb);

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

                int side_number = 0, sense = 1, offset = 0;
                MBERRCHK(mb->side_number(old_hex[0], (*qter), side_number, sense, offset), mb);

                if(m_GD==3){
                    get_normal_quad (adj_qconn, v);
                    if(sense == 1 && side_number >= 0){
                        // do nothing
                    }
                    else{
                        v=-v;
                    }
                }
                else if(m_GD==2){
                    if(sense == 1 && side_number >= 0){
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
        else if(all_bl[count] > 0 && fixmat == -1){ // node belongs to more than one material and fixmat not specified
            // NODE ON BOUNDARY
            // get the edge that is not on the boundary and count how many such edges we have
            moab::Range adj_for_normal;
            int nEdgeDim = 1;
            MBERRCHK(mb->get_adjacencies(&(*kter), 1, nEdgeDim, true, adj_for_normal, Interface::UNION), mb);
            moab::Range edge_normal = subtract(adj_for_normal, quads);

            if(edge_normal.size() > 1){
                m_LogFile << "MULTIPLE NORMALS ARE NEEDED" << std::endl;
            }
            else{
                m_LogFile << "We've one edge seperating materials 1 NORMAL IS NEEDED" << std::endl;
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
            m_LogFile << " Error, shouldn't have gotten here: " << count << std::endl;
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

            m_LogFile << std::setprecision (3) << std::scientific << " : NID:" << (nid)
                      << coords_old_quad[0]
                      << ", " << coords_old_quad[1] << ", " << coords_old_quad[2] << " OLD:- coords: NEW" << coords_new_quad[0]
                      << ", " << coords_new_quad[1] << ", " << coords_new_quad[2]  << std::endl;
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

            m_LogFile << std::setprecision (3) << std::scientific << " : NID:" << (nid)
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
                for(int i=0; i < fmconn.size(); i++){
                    if((*kter) == fmconn[i]){
                        //we have a node to unswap or reset connectivity for fixmat
                        fmconn[i] = new_vert[nid];
                    }
                }
                MBERRCHK(mb->set_connectivity(*fmter, &fmconn[0], fmconn.size()), mb);
            }

            //         m_LogFile << " We're here in MM case, now go along the edge --- have fun in the process !!" << std::endl;
            // material boundary - get edge direction and length

            //  find_min_edge_length(adj_qconn_r, qconn[i], nodes, m_MinEdgeLength);

            // check to see if this is a fixmat case
        }
        else if(all_bl[count] > 0 && fixmat == -1){ // node belongs to more than one material and fixmat not specified
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
        //        if(flag[qcount] != 1){

        MBERRCHK(mb->get_connectivity(&(*kter), 1, qconn),mb);
        double one_node_in_quad[3];
        for (int i=0; i<m_BElemNodes; i++){

            int node_tag_id = 0;
            MBERRCHK(mb->tag_get_data(BLNodeIDTag, &qconn[i], 1, &node_tag_id) ,mb);
            std::cout << "this is node id : " << node_tag_id << std::endl;

            MBERRCHK(mb->get_coords(&qconn[i], 1, one_node_in_quad),mb);
            m_LogFile << std::setprecision (3) << std::scientific << " new nodes:- coords: " << one_node_in_quad[0]
                      << ", " << one_node_in_quad[1] << ", " << one_node_in_quad[2]  << std::endl;



            //populate the connectivity after creating nodes for this BL node
            for(int j=0; j< m_Intervals; j++){
                if(m_Conn == 8 && m_BElemNodes == 4){ // hex
                    int nid = node_tag_id*m_Intervals + j;
                    if(j==0){
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
                    if(j==0){
                        conn[m_Conn*j +i] = qconn[m_BElemNodes-i-1];
                        conn[m_Conn*j + i+m_BElemNodes] = new_vert[nid + m_Intervals - 2];
                    }
                    else if(j==(m_Intervals-1)){
                        conn[m_Conn*j +i] = new_vert[nid - m_Intervals + 1];
                        conn[m_Conn*j + m_BElemNodes + 1 -i] = new_vert[nid];
                    }
                    else {
                        conn[m_Conn*j +i] = new_vert[nid + m_Intervals - 2*j -1];
                        conn[m_Conn*j + m_BElemNodes + 1 -i] = new_vert[nid + m_Intervals - 2*j -2];
                    }
                }
                else if(m_Conn == 4 && m_BElemNodes == 3 && hybrid == true){ // make wedges aka prisms for tet mesh
                    int nid = node_tag_id*m_Intervals+j;
                    if(j==0){
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
                    if(j==0){
                        conn[m_HConn*j +i] = qconn[m_BElemNodes-i-1];
                        conn[m_HConn*j + i+m_BElemNodes] = new_vert[nid + m_Intervals - 2];
                    }
                    else if(j==(m_Intervals-1)){
                        conn[m_HConn*j +i] = new_vert[nid - m_Intervals + 1];
                        conn[m_HConn*j + m_BElemNodes + 1 - i] = new_vert[nid];
                    }
                    else {
                        conn[m_HConn*j +i] = new_vert[nid + m_Intervals - 2*j -1];
                        conn[m_HConn*j + m_BElemNodes + 1 - i] = new_vert[nid + m_Intervals - 2*j -2];
                    }
                }
            }
        }

        //TODO: Set Connectivity of tet's, break prisms into 3 tets, Another loop is required.
        if(m_Conn == 3 && m_BElemNodes == 2 && hybrid == false){
            for(int c=0; c<m_Intervals; c++){
                if(tri_sch == 1){
                    // lower triangle
                    tri_conn[m_Conn*c*2] =     conn[c*m_HConn];
                    tri_conn[m_Conn*c*2 + 1] = conn[c*m_HConn + 1];
                    tri_conn[m_Conn*c*2 + 2] = conn[c*m_HConn + 3];
                    // upper triangle
                    tri_conn[m_Conn*c*2 + 3] = conn[c*m_HConn + 1];
                    tri_conn[m_Conn*c*2 + 4] = conn[c*m_HConn + 2];
                    tri_conn[m_Conn*c*2 + 5] = conn[c*m_HConn + 3];
                }
                else if(tri_sch == 2){
                    // lower triangle
                    tri_conn[m_Conn*c*2] =     conn[c*m_HConn];
                    tri_conn[m_Conn*c*2 + 1] = conn[c*m_HConn + 1];
                    tri_conn[m_Conn*c*2 + 2] = conn[c*m_HConn + 2];
                    // upper triangle
                    tri_conn[m_Conn*c*2 + 3] = conn[c*m_HConn + 2];
                    tri_conn[m_Conn*c*2 + 4] = conn[c*m_HConn + 3];
                    tri_conn[m_Conn*c*2 + 5] = conn[c*m_HConn];
                }
            }
        }

        // create boundary layer hexes
        for(int j=0; j< m_Intervals; j++){
            //double j_hex = 0.0;
            //  get_det_jacobian(conn, j*m_Conn, j_hex);
            if(m_Conn == 8){
                MBERRCHK(mb->create_element(MBHEX, &conn[j*m_Conn], m_Conn, hex),mb);
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
            // add this hex to a block
            if(mthis_set == 0){
                // here find the material set for this hex
                m_LogFile << "Find out the material that this hex corresponds to?? bailing out" << std::endl;
                exit(0);
            }
            else {
                MBERRCHK(mb->add_entities(mthis_set, &hex, 1), mb);
            }
            // mark entities for smoothing
            //                MBERRCHK(mb->add_entities(smooth_set, &hex, 1), mb);
            // add geom dim tag
            MBERRCHK(mb->add_entities(geom_set, &hex, 1), mb);
            //              // TODO: Add Local Smooting
        }

    }

    all_elems.clear();
    moab::Range skin_verts;
    MBERRCHK(mb->get_entities_by_dimension(0, 3, all_elems,true),mb);

    moab::Skinner skinner(mb);
    skinner.find_skin(0, all_elems, 0, skin_verts);

    m_LogFile << "setting 'fixed'' tag = 1 on verts in the skin = " <<  skin_verts.size() << std::endl;

    // set fixed tag = 1 on all skin verts
    std::vector<int> all_skin_data(skin_verts.size(), 1);
    MBERRCHK(mb->tag_set_data(FTag, skin_verts, &all_skin_data[0]), mb);

    m_LogFile << "\nTotal Jacobian calls/Min/Max: " << m_JacCalls << ", " << m_JLo << ", " << m_JHi << std::endl;

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
