#include "meshkit/PostBL.hpp"
namespace MeshKit
{
  int PostBL:: compute_normals()
  // ---------------------------------------------------------------------------
  //! Function: Compute the normal direction for boundary layer compution. \n
  //! Input:    From included .hpp file  \n
  //! Output:   none \n
  // ---------------------------------------------------------------------------
  {
    // COMPUTE NORMALS
    // get tag data and print
    all_bl.resize(nodes.size(), 0);
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
                
                m_LogFile << "Multiple normals are needed at: " << ncoord[0]
                          << ", " << ncoord[1] << ", " << ncoord[2] << " #normals " << edge_normal.size() << std::endl;
                exit(0);
              }
            else{
                m_LogFile << "We've one edge seperating materials 1, edge normal is needed" << edge_normal.size() << std::endl;
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
            m_LogFile << "Material must have an associated BLNode: Error! - shouldn't have gotten here: " << count << std::endl;
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
    
    return 0;
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
}
