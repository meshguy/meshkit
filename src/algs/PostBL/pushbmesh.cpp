#include "meshkit/PostBL.hpp"
namespace MeshKit
{
  int PostBL::push_bulk_mesh(VerdictWrapper vw)
  // ---------------------------------------------------------------------------
  //! Function: After normal creation push the bulk mesh and make room for creation of new BL elements  \n
  //! Input:    mesh quality verdict handle and included variables from PostBL.hpp file \n
  //! Output:   save a mesh file in debug mode \n
  // ---------------------------------------------------------------------------
  {
    // swap nodal coordinates of input nodes with innermost BL nodes to push the bulk mesh
    int count = -1;
    
    int matindx = -1;
    for (Range::iterator kter = nodes.begin(); kter != nodes.end(); ++kter){
        ++count;
        if(all_bl[count] == 0 ){
            
            int nid = (count+1)*m_Intervals - 1;
            MBERRCHK(mb->get_coords(&new_vert[nid], 1, coords_new_quad),mb);
            MBERRCHK(mb->get_coords(&(*kter), 1, coords_old_quad),mb);
            
            //TODO: Set connectivity for pushed hexes
            if (debug) {
                m_LogFile << std::setprecision (3) << std::scientific << " : NID: " << (nid)
                          << coords_old_quad[0]
                          << ", " << coords_old_quad[1] << ", " << coords_old_quad[2] << " OLD:- coords: NEW" << coords_new_quad[0]
                          << ", " << coords_new_quad[1] << ", " << coords_new_quad[2]  << std::endl;
              }
            moab::Range deformed_hex;
            MBERRCHK(mb->get_adjacencies(&(*kter), 1, m_GD, false, deformed_hex, Interface::UNION), mb);
            std::vector<moab::EntityHandle> dhex_conn;
            for(Range::iterator fmter = deformed_hex.begin(); fmter != deformed_hex.end(); ++fmter){
                mb->get_connectivity(&(*fmter), 1, dhex_conn);
                if(dhex_conn.size () > 0){
                    for(int i=0; i < (int)dhex_conn.size(); i++){
                        if((*kter) == dhex_conn[i]){
                            // push the bulk mesh
                            dhex_conn[i] = new_vert[nid];
                          }
                      }
                    MBERRCHK(mb->set_connectivity(*fmter, &dhex_conn[0], dhex_conn.size()), mb);
                  }
                double jac = 0;
                vw.quality_measure(*fmter, MB_JACOBIAN, jac);
                ++m_JacCalls;
                if (jac < 0){
                    m_LogFile << "ck BL thickness/intervals. Stopping." << std::endl;
                    exit(0);
                  }
              }
          }
        else if(all_bl[count] > 0 && fixmat != -1){ 
            // node belongs to more than one material and fixmat specified
            int nid = (count+1)*m_Intervals - 1;
            MBERRCHK(mb->get_coords(&new_vert[nid], 1, coords_new_quad),mb);
            MBERRCHK(mb->get_coords(&(*kter), 1, coords_old_quad),mb);

            if (debug) {
                m_LogFile << std::setprecision (3) << std::scientific << " : NID: " << (nid)
                          << coords_old_quad[0]
                          << ", " << coords_old_quad[1] << ", " << coords_old_quad[2] << " OLD:- coords: NEW" << coords_new_quad[0]
                          << ", " << coords_new_quad[1] << ", " << coords_new_quad[2]  << std::endl;
              }
            
            // now find the hex in fixmat_ents that was affected by this swap and UNswap the coords
            moab::Range fhex;
            MBERRCHK(mb->get_adjacencies(&(*kter), 1, m_GD, false, fhex, Interface::UNION), mb);
            moab::Range fmhex= intersect(fhex,fixmat_ents);
            moab::Range non_fixhex = subtract(fhex, fmhex);
            
            std::vector<int> tag_non_fixhex(non_fixhex.size(),0);
            MBERRCHK(mb->tag_get_data(MatIDTag, non_fixhex, &tag_non_fixhex[0]), mb);
            m_LogFile << quads.size() << std::endl;
         
            if( matindx < (int) quads.size() ){
                ++matindx;
                // handling more than 2 materials in MM case is not supported
                blmaterial_id[matindx] = tag_non_fixhex[0];
              }
            std::vector<moab::EntityHandle> fmconn;
            for(Range::iterator fmter = non_fixhex.begin(); fmter != non_fixhex.end(); ++fmter){
                mb->get_connectivity(&(*fmter), 1, fmconn);
                if(fmconn.size () > 0){
                    for(int i=0; i < (int)fmconn.size(); i++){
                        if((*kter) == fmconn[i]){
                            // push the bulk mesh
                            fmconn[i] = new_vert[nid];
                          }
                      }
                    MBERRCHK(mb->set_connectivity(*fmter, &fmconn[0], fmconn.size()), mb);
                  }
                double jac = 0;
                vw.quality_measure(*fmter, MB_JACOBIAN, jac);
                ++m_JacCalls;
                if (jac < 0){
                    m_LogFile << "ck BL thickness/intervals. Stopping." << std::endl;
                    exit(0);
                  }
              }   
            m_LogFile << "Multiple material case: working along the edge" << std::endl;
          }
        else if(all_bl[count] > 0 && fixmat == -1){ // node belongs to more than one material and fixmat not specified
            // NODE ON BOUNDARY
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
      } // for loop ends
    
    if(debug == true){
        MBERRCHK(mb->write_mesh("bulkpushed.exo"),mb);
        m_LogFile <<  "\n\nWrote Mesh File: bulkpushed.exo" << std::endl;
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
    return 0;
  }
}
