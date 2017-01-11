#include "meshkit/PostBL.hpp"
namespace MeshKit
{
  int PostBL::create_bl_elements(VerdictWrapper vw)
  // ---------------------------------------------------------------------------
  //! Function: Parser for reading the PostBL specification (.inp) file. \n
  //! Input:    Command line arguments. \n
  //! Output:   none \n
  // ---------------------------------------------------------------------------
  {
    int qcount = 0;
    
    // Now start creating New elements
    for (Range::iterator kter = quads.begin(); kter != quads.end(); ++kter){
        qcount++;        
        std::vector<moab::EntityHandle> old_hex;
        MBERRCHK(mb->get_adjacencies(&(*kter), 1, m_GD, false, old_hex),mb);
        //      int this_quad_mat = -1;
        //      MBERRCHK(mb->tag_get_data(MatIDTag, &old_hex[0], 1, &this_quad_mat), mb);
        double jac = 0;
        
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
                    if(j == 0){
                        conn[m_Conn*j +i] = qconn[i];
                        conn[m_Conn*j + i+m_BElemNodes] = new_vert[nid];
                      }
                    else if(j==(m_Intervals-1)){
                        conn[m_Conn*j +i] = new_vert[nid];
                        conn[m_Conn*j + i+m_BElemNodes] = new_vert[nid -1];
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
                    if(j==0){
                        conn[m_HConn*j +i] = qconn[i];
                        conn[m_HConn*j + i+m_BElemNodes] = new_vert[nid];
                      }
                    else if(j==(m_Intervals-1)){
                        conn[m_HConn*j +i] = new_vert[nid];
                        conn[m_HConn*j + i+m_BElemNodes] = new_vert[nid-1];
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
        // hex material will be found from the innermost hex adjacency. Let's create that first.       
        for(int j=0; j < m_Intervals; j++){
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
            
            vw.quality_measure(hex, MB_JACOBIAN, jac);
            if(jac < 0){
                m_LogFile <<  "New BL elements have negative jacobian, trying to fix.. " << jac << std::endl;
                // swap 0-4 1-5 2-6 3-7 for positive jacobian -> posibly a multi material case
                moab::EntityHandle temp;
                for(int ii = j*m_Conn; ii < j*m_Conn + m_BElemNodes; ii++){
                    temp = conn[ii];
                    conn[ii] = conn[ii+m_BElemNodes];
                    conn[ii+m_BElemNodes] = temp;
                  }
                
                //reverse the tet connectivity based on the prism reverse
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
                
                // recreate this element 
                if(m_Conn == 8){
                    MBERRCHK(mb->create_element(moab::MBHEX, &conn[j*m_Conn], m_Conn, hex),mb);
                  }
                else if(m_Conn==4 && m_GD ==3 && hybrid == false){
                    MBERRCHK(mb->create_element(MBTET, &tet_conn[j*m_Conn*3], m_Conn, hex),mb);
                    MBERRCHK(mb->create_element(MBTET, &tet_conn[j*m_Conn*3+4], m_Conn, hex1),mb);
                    MBERRCHK(mb->create_element(MBTET, &tet_conn[j*m_Conn*3+8], m_Conn, hex2),mb);
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
                // check jacobian again
                vw.quality_measure(hex, MB_JACOBIAN, jac);              
                if(jac < 0){
                    m_LogFile << "Negative jacobian still found for the Bl elements, invalid mesh!"  << jac << std::endl;
                  }      
                // end recreation of element -> negative hex was encontered while creating Bl elements
              }
            
            if(mthis_set == 0){ // No material specified for BL hexes find the material to append the hex             
                moab::EntityHandle mat_set = 0;
                for (set_it = m_sets.begin(); set_it != m_sets.end(); set_it++)  {
                    int set_id = -1;
                    mat_set = *set_it;
                    MBERRCHK(mb->tag_get_data(MTag, &mat_set, 1, &set_id), mb);
                    if(set_id == blmaterial_id[qcount -1])
                      break;
                  }
                if(mat_set != 0){
                    if(hex!=0)
                      MBERRCHK(mb->add_entities(mat_set,&hex,1), mb);
                    if(hex1!=0)
                      MBERRCHK(mb->add_entities(mat_set,&hex1,1), mb);
                    if(hex2!=0)
                      MBERRCHK(mb->add_entities(mat_set,&hex2,1), mb);
                  }
                else{
                    exit(0);
                  }
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
            // mark entities for smoothing - careful: smoothing BL elements will cause loss of biased elements
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
    return 0;
  }
}
