#include <math.h>
#include <iomanip>
#include "../meshkit/Pbl.hpp"
#ifdef HAVE_MOAB
#include "MBSkinner.hpp"
#include "MBAdaptiveKDTree.hpp"
#include "MBRange.hpp"
#include "MBCartVect.hpp"
#include "moab/Matrix3.hpp"
#endif

namespace MeshKit
{

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
	debug = false;
        m_NeumannSet = -1;
	m_Material = 999;
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
            m_LogFile <<  "\nIn setup this : " <<  "\n";
        }
        
    }
    
    void Pbl::execute_this()
    // ---------------------------------------------------------------------------
    // Function: Obtains parameters for post meshing boundary layer and load mesh file
    // Input:    command line arguments
    // Output:   none
    // ---------------------------------------------------------------------------
    {
      std::cout << "\nIn execute this : creating boundary layer elements.." <<  "\n";
      m_LogFile << "\nIn execute this : creating boundary layer elements.." <<  "\n";
        // start the timer
        CClock Timer;
        clock_t sTime = clock();
        std::string szDateTime;
        Timer.GetDateTime (szDateTime);

        std::cout << "\nStarting out at : " << szDateTime << "\n";   
        std::cout << "\n Loading meshfile: " << m_MeshFile << ".." << std::endl;
        m_LogFile <<  "\nStarting out at : " << szDateTime << "\n";
        m_LogFile <<  "\n Loading meshfile: " << m_MeshFile << ".." << std::endl;

        //load specified mesh file
        IBERRCHK(imesh->load(0, m_MeshFile.c_str(),0), *imesh);
        m_GD = imesh->getGeometricDimension();
        m_LogFile << "Geometric dimension of meshfile = "<< m_GD <<std::endl;
        std::cout  <<"Geometric dimension of meshfile = "<< m_GD <<std::endl;
 
        // obtain the boundary layer surface faces
        moab::Tag GDTag, GIDTag, NTag, MTag;
        MBERRCHK(mb->tag_get_handle("GEOM_DIMENSION", 1, moab::MB_TYPE_INTEGER, GDTag),mb);
        MBERRCHK(mb->tag_get_handle("NEUMANN_SET", 1, moab::MB_TYPE_INTEGER, NTag),mb);
        MBERRCHK(mb->tag_get_handle("MATERIAL_SET", 1, moab::MB_TYPE_INTEGER, MTag),mb);
        MBERRCHK(mb->tag_get_handle("GLOBAL_ID", 1, moab::MB_TYPE_INTEGER, GIDTag),mb);
        
        moab::Range sets, verts, n_sets, m_sets;
        int geom_dim = 2;
        const void* gdim[] = {&geom_dim};
        
        // get all the entity sets with Geom Dimension=2
        MBERRCHK(mb->get_entities_by_type_and_tag(0, moab::MBENTITYSET, &GDTag, 
						  gdim, 1 , sets, moab::Interface::INTERSECT, false), mb);
        
        MBERRCHK(mb->get_entities_by_type_and_tag(0, moab::MBENTITYSET, &NTag, 0, 1 , n_sets),mb);

        MBERRCHK(mb->get_entities_by_type_and_tag(0, moab::MBENTITYSET, &MTag, 0, 1 , m_sets),mb);
        
	// Handling NeumannSets (if BL surf in input via NS)
        moab::Range::iterator set_it;
        moab::EntityHandle this_set;
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
            std::cout << "Looking for NS with id " << m_NeumannSet << 
	      ". Total NS found are: "<< n_sets.size() << std::endl;
            m_LogFile <<  "Looking for NS with id " << m_NeumannSet << 
	      ". Total NS found are: "<< n_sets.size() << std::endl;
	}

        // get the  all the quads and nodes in a range
        moab::EntityHandle s1;
        MBRange quads, nodes;
        int dims; // variable to store global id of all sets with GD=2

	// INPUT by NeumannSet
        if(m_NeumannSet != -1 && this_set != 0){
	  MBERRCHK(mb->get_entities_by_type(this_set, moab::MBQUAD, quads,true),mb);
	    if (quads.size() <=0){
	      std::cout << " No quads found, aborting.. " << std::endl;
	      m_LogFile <<  " No quads found, aborting.. " << std::endl;
	      exit(0);
	    }
            MBERRCHK(mb->get_adjacencies(quads, 0, false, nodes, MBInterface::UNION),mb);
            if (debug) {
                std::cout << "Found NeumannSet with id : " << m_NeumannSet <<  std::endl;
                std::cout << "#Quads in this surface: " << quads.size() << std::endl;
                std::cout << "#Nodes in this surface: " << nodes.size() << std::endl;
		std::cout << "#New Nodes to be created:" << m_Intervals*nodes.size() << std::endl;
                m_LogFile <<  "Found NeumannSet with id : " << m_NeumannSet <<  std::endl;
                m_LogFile <<  "#Quads in this surface: " << quads.size() << std::endl;
                m_LogFile <<  "#Nodes in this surface: " << nodes.size() << std::endl;
		m_LogFile << "#New nodes to be created:" << m_Intervals*nodes.size() << std::endl;
            }
        }
	// INPUT by surface id 
        else if (m_SurfId !=-1){
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
			std::cout << "#New Nodes to be created:" << m_Intervals*nodes.size() << std::endl;
                        m_LogFile <<  "Found surface with id : " << m_SurfId <<  std::endl;
                        m_LogFile <<  "#Quads in this surface: " << quads.size() << std::endl;
                        m_LogFile <<  "#Nodes in this surface: " << nodes.size() << std::endl;
			m_LogFile << "#New nodes to be created:" << m_Intervals*nodes.size() << std::endl;
                    }
                }
            }
        }
	else {
	  std::cout << "Invalid boundary layer specification, aborting.." << std::endl;
	  m_LogFile <<  "Invalid boundary layer specification, aborting.." <<  std::endl; 
	  exit(0);
        }
	
	// Handling MaterialSet
	moab::Range::iterator mset_it;
        moab::EntityHandle mthis_set;
	int mset_id = 0, found = 0;
        for (mset_it = m_sets.begin(); mset_it != m_sets.end(); mset_it++)  {
            
            mthis_set = *mset_it;
            
            // get entity handle of MS specified in the input file
            MBERRCHK(mb->tag_get_data(MTag, &mthis_set, 1, &mset_id), mb);
            if(mset_id == m_Material){
	      found = 1;
	      break;
	    }
            mthis_set = 0;
        }
	if(found == 1 && m_Material !=999){
	  std::cout << "Found material set with id " << m_Material << std::endl;
	  m_LogFile << "Found material set with id " << m_Material << std::endl;
	}
	else{
	  // No material set found, creating material set 999 for BL elements
	  std::cout << "Creating material set with id " << m_Material << std::endl;
	  m_LogFile << "Creating material set with id " << m_Material << std::endl; 
	  MBERRCHK(mb->create_meshset(moab::MESHSET_SET, mthis_set, 1), mb);
	  MBERRCHK(mb->tag_set_data(MTag, &mthis_set, 1, &m_Material), mb);
	}

        std::vector <bool> node_status(false); // size of verts of bl surface
        node_status.resize(nodes.size());
        MBRange edges, hexes, hex_edge, quad_verts, adj_quads, adj_hex_nodes;
        std::vector<EntityHandle> conn(m_Intervals*8), qconn(4), 
	  new_vert(m_Intervals*nodes.size()), old_hex, old_hex_conn(8), adj_hexes, adj_hex_nodes1;
        double coords_bl_quad[3], coords_new_quad[3], xdisp = 0.0, ydisp = 0.0, zdisp = 0.0;
        EntityHandle hex;
        int qcount = 0;
        int ncount = 0;

	// Now start creating New elements
        for (Range::iterator kter = quads.begin(); kter != quads.end(); ++kter){
            qcount++;
            MBERRCHK(mb->get_connectivity(&(*kter), 1, qconn),mb);
            
            get_normal_quad (qconn, xdisp, ydisp, zdisp);
            if (debug){
                std::cout << "\n\n*** QUAD: " << qcount << "\n\ndirection vector: " << xdisp 
			  << " " << ydisp << " " << zdisp << " " << std::endl;
                m_LogFile <<  "\n\n*** QUAD: " << qcount << "\n\ndirection vector: " << xdisp 
			  << " " << ydisp << " " << zdisp << " " << std::endl;
	    }
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
			break;
                    }
                }
                if (debug){
                    std::cout << "\n*Working on Node: " << tmp << " coords: " << coords_bl_quad[0] 
				    << ", " << coords_bl_quad[1] << ", " << coords_bl_quad[2]  << std::endl;
                    m_LogFile <<  "\n*Working on Node: " << tmp << " coords: " << coords_bl_quad[0] 
				    << ", " << coords_bl_quad[1] << ", " << coords_bl_quad[2]  << std::endl;
		}
		double temp;

                if(node_status[tmp] == false){
		  adj_hexes.clear();
		  MBERRCHK(mb->get_adjacencies(&qconn[i], 1, 3, true, adj_hexes, MBInterface::UNION), mb);
		  std::cout << adj_hexes.size() << " - adjacent hexes found " << std::endl;

		  double num = m_Thickness*(m_Bias-1)*(pow(m_Bias, m_Intervals -1));
		  double deno = pow(m_Bias, m_Intervals) - 1;
		  if (deno !=0)
		    temp = num/deno;
		  else
		    temp = m_Thickness/m_Intervals;
		  // another loop to number of boundary layers
		  double move = 0.0;
                    for(int j=0; j< m_Intervals; j++){
		      
		      move+= temp/pow(m_Bias,j);
		      if (debug){
			std::cout <<  " move :" << move;
			  m_LogFile <<  " move:" << move;
		      }
		      // now compute the coords of the new vertex
		      coords_new_quad[0] = coords_bl_quad[0]-move*xdisp;
		      coords_new_quad[1] = coords_bl_quad[1]-move*ydisp;
		      coords_new_quad[2] = coords_bl_quad[2]-move*zdisp;
		      
                        int nid = tmp*m_Intervals+j;
                        mb->create_vertex(coords_new_quad, new_vert[nid]);
                        if (debug){
			  std::cout << std::setprecision (3) << std::scientific << " : created node:" << nid 
				    << " of " << new_vert.size() << " new nodes:- coords: " << coords_new_quad[0] 
				    << ", " << coords_new_quad[1] << ", " << coords_new_quad[2]  << std::endl;
			  m_LogFile << std::setprecision (3) << std::scientific << " : created node:" << nid 
				    << " of " << new_vert.size() << " new nodes:- coords: " << coords_new_quad[0] 
				    << ", " << coords_new_quad[1] << ", " << coords_new_quad[2]  << std::endl;
			}
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
		// if a hex does have a quad on BL, but, has a node or edge on the boundary layer
		for(int k=0; k < (int) adj_hexes.size(); k++){
		  adj_hex_nodes1.clear();		    
		  MBERRCHK(mb->get_connectivity(&adj_hexes[k], 1, adj_hex_nodes1), mb);
		  std::vector<EntityHandle> inodes;
		  int var = 0;
		  for(int n = 0; n < nodes.size(); n++){
		    for(int m = 0; m < 8; m++){
		      if(nodes[n] == adj_hex_nodes1[m]){
			inodes.push_back(adj_hex_nodes1[m]);
			var++;
		      }
		    }
		  }
		  // doesn't have a quad or when inodes is 0 it mean this is a newly created hex
		  if(inodes.size() != 4 && inodes.size() != 0){
		    if(debug){
		      std::cout << "Hex on BL surface is connected by a node or edge" <<
			"\n nodes on BL - " << inodes.size() << " aborting .." << std::endl;
		      m_LogFile << "Hex on BL surface is connected by a node or edge" <<  
			"\n nodes on BL - " << inodes.size() << " aborting .." <<std::endl;
		    }
		    // TODO: mark this hex and set it's connectivity later
		    int t = 0;
		    for(int p = 0; p < 8; p++){
			if(qconn[i] == adj_hex_nodes1[p]){
			    adj_hex_nodes1[p] = conn[8*(m_Intervals-1)+i];
		      }
		    }	  
		    MBERRCHK(mb->set_connectivity(adj_hexes[k], &adj_hex_nodes1[0], 8), mb);
		  }
		    inodes.clear();
		}
		
            } // Loop thru 4 nodes of a quad ends

            if (debug)
	      {
                m_LogFile <<  "\nsetting connectivity of the old BL hex " << std::endl;
                std::cout << "\nsetting connectivity of the old BL hex " << std::endl;
	      }
            // set connectivity of the old_hex
            for(int p=0; p<8; p++){
                for(int q=0; q<4; q++){
                    if (old_hex_conn[p] == qconn[q]){
                        old_hex_conn[p] = conn[8*(m_Intervals-1) + q];
                    }
                }
            }
            double j_old_hex = 0.0;
            get_det_jacobian(old_hex_conn, 0, j_old_hex);
            if (j_old_hex <= 0){
                std::cout << "\n ERRROR: Negative Jacobian of old hex\n Bailing out.. ";
		m_LogFile <<  "\n ERRROR: Negative Jacobian of old hex\n Bailing out.. ";
                exit(0);
            }
            
            MBERRCHK(mb->set_connectivity(old_hex[0], &old_hex_conn[0], 8), mb);
            old_hex.clear();
	    
	    // handle 3 or 4 hex adjacency cases, set connectivity of non BL hexes.
	    
            if (debug){
                std::cout << "creating new boundary layer hexes" << std::endl;
                m_LogFile <<  "creating new boundary layer hexes" << std::endl;
	    }
            // create boundary layer hexes
            for(int j=0; j< m_Intervals; j++){
                double j_hex = 0.0;
                get_det_jacobian(conn, j*8, j_hex);
                if (j_old_hex <= 0){
                    std::cout << "\n ERRROR: Negative Jacobian of old hex\n Bailing out.. ";
                    m_LogFile <<  "\n ERRROR: Negative Jacobian of old hex\n Bailing out.. ";
                    exit(0);
                }
                MBERRCHK(mb->create_element(MBHEX, &conn[j*8], 8, hex),mb);
		// add this hex to a block
		MBERRCHK(mb->add_entities(mthis_set, &hex, 1), mb);
            }
        } // Loop thru quads ends
        
        //save the final boundary layer mesh
        MBERRCHK(mb->write_mesh(m_OutFile.c_str()),mb);
        std::cout << "\n\nWrote Mesh File: " << m_OutFile << std::endl;
        m_LogFile <<  "\n\nWrote Mesh File: " << m_OutFile << std::endl;
        // get the current date and time
        Timer.GetDateTime (szDateTime);
        std::cout << "Ending at : " << szDateTime;
	m_LogFile << "Ending at : " << szDateTime;
        // report/compute the elapsed time
        std::cout << "Elapsed wall clock time: " << Timer.DiffTime ()
        << " seconds or " << (Timer.DiffTime ())/60.0 << " mins\n";
        m_LogFile <<  "Elapsed wall clock time: " << Timer.DiffTime ()
        << " seconds or " << (Timer.DiffTime ())/60.0 << " mins\n";
	std::cout << "Total CPU time used: " << (double) (clock() - sTime)/CLOCKS_PER_SEC \
        << " seconds" << std::endl;        
        m_LogFile <<  "Total CPU time used: " << (double) (clock() - sTime)/CLOCKS_PER_SEC \
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
        std::cout << "\t\t         Tool to generate Post-mesh Boundary Layers      " << '\n';
        std::cout << "\t\t\t\tArgonne National Laboratory" << '\n';
        std::cout << "\t\t\t\t        2012         " << '\n';
        std::cout << "\t\t---------------------------------------------------------" << '\n';
        std::cout << "\nsee README file for using the program and details on various cards.\n"<< std::endl;
        // set and open input output files
        bool bDone = false;
        do{
            if (2 == argc) {
                m_InputFile = argv[1];
		m_LogName = m_InputFile + ".log";
            }
            else if (1 == argc){
                std::cout << "\nRunning default case:\n" << std::endl;
                m_LogFile << "\nRunning default case:\n" << std::endl;

                m_InputFile = TestDir + "/" + (char *)DEFAULT_TEST_PBL;
		m_LogName = m_InputFile + ".log";
            }
            
            // open input file for reading
            m_FileInput.open (m_InputFile.c_str(), std::ios::in);
            if (!m_FileInput){
                std::cout << "Unable to open file: " << m_InputFile << std::endl;
                m_LogFile << "Unable to open file: " << m_InputFile << std::endl;

                m_FileInput.clear ();
                exit(1);
            }
            else
                bDone = true; // file opened successfully

	    // open the log file for dumping debug/output statements
            m_LogFile.open (m_LogName.c_str(), std::ios::out);
            if (!m_LogFile){
                std::cout << "Unable to open file: " << m_LogName << std::endl;
                m_LogFile <<  "Unable to open file: " << m_LogName << std::endl;
                m_LogFile.clear ();
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
            
            // Get MeshFile name
            if (szInputString.substr(0,8) == "meshfile"){
                std::istringstream szFormatString (szInputString);
                szFormatString >> m_Card >> m_MeshFile;
                if(szFormatString.fail())
                    IOErrorHandler(INVALIDINPUT);
                std::cout << m_Card << " name read: "<< m_MeshFile << std::endl;
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
                std::cout  << m_Card << " read: "<< m_SurfId <<std::endl;
                m_LogFile <<  m_Card << " read: "<< m_SurfId <<std::endl;
            }
            // Get BL surface via neumann set or sideset
            if (szInputString.substr(0,10) == "neumannset"){
                std::istringstream szFormatString (szInputString);
                szFormatString >> m_Card >> m_NeumannSet;
                if(m_NeumannSet < 0 || szFormatString.fail())
                    IOErrorHandler(INVALIDINPUT);
                std::cout  << m_Card << " read: "<< m_NeumannSet <<std::endl;
                m_LogFile <<  m_Card << " read: "<< m_NeumannSet <<std::endl;
            }
              // Get BL material (block) number
            if (szInputString.substr(0,8) == "material"){
                std::istringstream szFormatString (szInputString);
                szFormatString >> m_Card >> m_Material;
                if(m_Material < 0 || szFormatString.fail())
                    IOErrorHandler(INVALIDINPUT);
                std::cout  << m_Card << " read: "<< m_Material <<std::endl;
                m_LogFile <<  m_Card << " read: "<< m_Material <<std::endl;
            }
                      
            // Get thickness
            if (szInputString.substr(0,9) == "thickness"){
                std::istringstream szFormatString (szInputString);
                szFormatString >> m_Card >> m_Thickness;
                if(m_Thickness < 0 || szFormatString.fail())
                    IOErrorHandler(INVALIDINPUT);
                std::cout  << m_Card << " read: "<< m_Thickness <<std::endl;
                m_LogFile <<  m_Card << " read: "<< m_Thickness <<std::endl;
            }
            // Get intervals
            if (szInputString.substr(0,9) == "intervals"){
                std::istringstream szFormatString (szInputString);
                szFormatString >> m_Card >> m_Intervals;
                if(m_Intervals < 0 || szFormatString.fail())
                    IOErrorHandler(INVALIDINPUT);
                std::cout  << m_Card << " read: "<< m_Intervals <<std::endl;
                m_LogFile <<  m_Card << " read: "<< m_Intervals <<std::endl;
            }
            // Get bias
            if (szInputString.substr(0,4) == "bias"){
                std::istringstream szFormatString (szInputString);
                szFormatString >> m_Card >> m_Bias;
                if(m_Bias < 0 || szFormatString.fail())
                    IOErrorHandler(INVALIDINPUT);
                std::cout  << m_Card << " read: "<< m_Bias <<std::endl;
                m_LogFile <<  m_Card << " read: "<< m_Bias <<std::endl;
            }
            // Output file name
            if (szInputString.substr(0,7) == "outfile"){
                std::istringstream szFormatString (szInputString);
                szFormatString >> m_Card >> m_OutFile;
                if(szFormatString.fail())
                    IOErrorHandler(INVALIDINPUT);
                std::cout  << m_Card << " name read: "<< m_OutFile <<std::endl;
                m_LogFile <<  m_Card << " name read: "<< m_OutFile <<std::endl;
            }
            // Debug flag
            if (szInputString.substr(0,5) == "debug"){
                std::istringstream szFormatString (szInputString);
                szFormatString >> m_Card >> debug;
                if(szFormatString.fail())
                    IOErrorHandler(INVALIDINPUT);
                std::cout  << m_Card << " name read: "<< debug <<std::endl;
                m_LogFile <<  m_Card << " name read: "<< debug <<std::endl;
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
    
    void Pbl::get_det_jacobian(std::vector<EntityHandle>conn, int offset, double &detJ)
    // ---------------------------------------------------------------------------
    // Function: get normal of a quad
    // Input:    conn
    // Output:   vector x, y and z
    // ---------------------------------------------------------------------------
    {
        CartVect vertex[8], xi(0.0);
        MBERRCHK(mb->get_coords(&conn[offset], 8, (double*) &vertex[0]), mb);
        
        double corner[8][3] = { { -1, -1, -1 },
            {  1, -1, -1 },
            {  1,  1, -1 },
            { -1,  1, -1 },
            { -1, -1,  1 },
            {  1, -1,  1 },
            {  1,  1,  1 },
            { -1,  1,  1 } };
        Matrix3 J(0.0);
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
    }
    
} // namespace MeshKit


