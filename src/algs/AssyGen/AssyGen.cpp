/*********************************************
AssyGen Tool: Reactor Geometry Generator
Argonne National Laboratory
*********************************************/

#include "meshkit/AssyGen.hpp"

namespace MeshKit
{
  // static registration of this  mesh scheme
  moab::EntityType AssyGen_tps[] = { moab::MBVERTEX,
                                     moab::MBEDGE,
                                     moab::MBTRI,
                                     moab::MBHEX,
                                     moab::MBMAXTYPE};
  const moab::EntityType* AssyGen::output_types()
  { return AssyGen_tps; }

  AssyGen::AssyGen( MKCore *mk, const MEntVector &me_vec)
    : MeshScheme( mk, me_vec),
      igeomImpl(mk->igeom_instance())
  {
    err = 0;
    tmpSB = 1;
    m_nPlanar = 0; //default is 3D
    m_nLineNumber = 0;
    root_set= NULL;
    szComment = "!";
    MAXCHARS = 10000;
    MAXLINES = 1000;
    pi = M_PI;
    m_nTotalPincells = 0;
    m_dRadialSize = -1.0;
    m_dTetMeshSize = -1.0;
    m_nDimensions = 0;
    m_nMaterialSetId = 1;
    m_nNeumannSetId = 1;
    m_szEngine = "acis";
    m_szMeshType = "hex";
    m_nDuct = 0;
    m_nDuctNum = 0;
    m_nJouFlag = 0;
    m_szSideset = "yes";
    m_nAssyGenInputFiles = 0;
    m_dMergeTol = 1e-4;
    m_edgeInterval = 99;
    m_nStartpinid = 1;
    m_szInfo = "off";
    m_szMeshScheme = "pave";
    pin_name = "";
    m_nHblock = -1;
    m_nPincells = 0;
    m_bCreateMatFiles = false;
    m_nSuperBlocks = 0;
    m_bmerge = false;
    m_bimprint = false;
    save_exodus = false;
    have_common = true;
    com_run_count = 0;
    m_nBLAssemblyMat = 0;
    m_szInnerDuct = "";
    m_szSmooth  = "off";
  }

  AssyGen::~AssyGen()
  {
    iGeom_dtor(igeomImpl->instance(), &err);
    //CHECK( "Interface destruction didn't work properly." );
    // close the input and output files
    m_FileInput.close ();
    m_FileOutput.close ();
    m_SchemesFile.close ();
    if(strcmp(m_szInfo.c_str(),"on") == 0)
      m_AssmInfo.close ();
  }


  bool AssyGen::add_modelent(ModelEnt *model_ent)
  {
    return MeshOp::add_modelent(model_ent);
  }

  void AssyGen::setup_this()
  {

    // start the timer
    CClock Timer;
    clock_t sTime = clock();
    std::string szDateTime;
    Timer.GetDateTime (szDateTime);
    std::cout << "\nStarting out at : " << szDateTime << "\n";

    if (have_common == true)
      ReadCommonInp();

    //count pin cylinders and cell material, needed for setting array size before actual read
    ReadInputPhase1 ();

    if (have_common == true)
      ReadCommonInp();

    // read the problem size and create pincell
    ReadAndCreate ();

    // create the .jou file
    CreateCubitJournal();

    CreateAssyGenInputFiles();

    // get the current date and time
    Timer.GetDateTime (szDateTime);
    std::cout << "Ending at : " << szDateTime;

    // compute the elapsed time
    std::cout << "Elapsed wall clock time: " << Timer.DiffTime ()
              << " seconds or " << (Timer.DiffTime ())/60.0 << " mins\n";

    std::cout << "## Total CPU time used := " << (double) (clock() - sTime)/CLOCKS_PER_SEC \
              << " seconds" << std::endl;
  }

  void AssyGen::execute_this()
  {
  }


  void AssyGen::PrepareIO (int argc, char *argv[],  std::string TestDir)
  // ---------------------------------------------------------------------------
  // Function: Obtains file names and opens input/output files
  // Input:    command line arguments
  // Output:   none
  // ---------------------------------------------------------------------------
  {
    std::cout << '\n';
    std::cout << "\t\t---------------------------------------------------------" << '\n';
    std::cout << "\t\tProgram to Generate Nuclear Reactor Assembly Geometries      " << '\n';
    std::cout << "\t\t\t\tArgonne National Laboratory" << '\n';
    std::cout << "\t\t---------------------------------------------------------" << '\n';
    std::cout << "\nsee http://press3.mcs.anl.gov/sigma/meshkit-library/rgg/ for details.\n"<< std::endl;
    // set and open input output files
    bool bDone = false;
    do{
        if (2 == argc) {
            m_szFile = argv[1];
            m_szInFile=m_szFile+".inp";
            m_szJouFile = m_szFile+".jou";
            m_szSchFile = m_szFile+".template.jou";
            m_szAssmInfo = m_szFile + "_info.csv";
            m_szLogFile = m_szFile + ".scrŒeenlog";
            m_szPyCubGeom = m_szFile + ".py";
            m_szCommonFile = "common.inp";
          }
        else if (3 == argc) {
            int i=1;// will loop through arguments, and process them
            for (i=1; i<argc-1 ; i++) {
                if (argv[i][0]=='-') {
                    switch (argv[i][1])
                      {
                      case 'j':
                        {
                          m_nJouFlag = 1;
                          std::cout << "Creating journal file only.\n Geometry file must exist in the same directory." << std::endl;
                          m_szFile = argv[2];
                          m_szInFile=m_szFile+".inp";
                          m_szJouFile = m_szFile+".jou";
                          m_szSchFile = m_szFile+".template.jou";
                          m_szAssmInfo = m_szFile + "_info.csv";
                          m_szLogFile = m_szFile + ".screenlog";
                          m_szPyCubGeom = m_szFile + ".py";
                          m_szCommonFile = "common.inp";
                          break;
                        }
                      case 'h':
                        {
                          std::cout << "\nInstruction on writing assygen input file can also be found at: " << std::endl;
                          std::cout << "        http://press3.mcs.anl.gov/sigma/meshkit/rgg/assygen-input-file-keyword-definitions/" << std::endl;
                          std::cout << "Usage: assygen [-j -h] <input file name without extension>"<< std::endl;
                          std::cout << "        -j create journal file only" << std::endl;
                          std::cout << "        -h print help" << std::endl;

                          exit(0);
                          break;
                        }
                      }
                  }
              }
          }
        else if (1 == argc){
            std::cout << "\nInstruction on writing assygen input file can also be found at: " << std::endl;
            std::cout << "        http://press3.mcs.anl.gov/sigma/meshkit/rgg/assygen-input-file-keyword-definitions/" << std::endl;
            std::cout << "Usage: assygen [-t -j -h] <input file name without extension>"<< std::endl;
            std::cout << "        -t print timing and memory usage info in each step" << std::endl;
            std::cout << "        -j create journal file only" << std::endl;
            std::cout << "        -h print help" << std::endl;

            m_szInFile = TestDir + "/" + (char *)DEFAULT_TEST_FILE;
            m_szGeomFile = (char *)TEST_FILE_NAME;
            m_szJouFile = (char *)TEST_FILE_NAME;
            m_szFile =  (char *)TEST_FILE_NAME;
            m_szInFile+=".inp";
            m_szJouFile+=".jou";
            m_szSchFile = m_szFile+".template.jou";
            m_szAssmInfo = m_szFile + "_info.csv";
            m_szLogFile = m_szFile + ".screenlog";
            m_szPyCubGeom = m_szFile + ".py";
            m_szCommonFile = TestDir + "/" + "common.inp";

            std::cout <<"Default case input file is located here <MeshKit/data> "<< std::endl;
          }
        // open the file
        m_FileInput.open (m_szInFile.c_str(), std::ios::in);
        if (!m_FileInput){
            std::cout << "Usage: assygen <input filename WITHOUT EXTENSION>"<< std::endl;
            m_FileInput.clear ();
            exit(1);
          }
        else
          bDone = true; // file opened successfully

        // open common.inp file, if not found do nothing.
        m_FileCommon.open (m_szCommonFile.c_str(), std::ios::in);
        if (!m_FileCommon){
            have_common = false;
            std::cout << "common.inp file not specified." << std::endl;
            m_FileCommon.clear ();
          }
        else {
            have_common = true;
          }
        std::cout << " opened file " << m_szCommonFile << " have common is "
                  << have_common << std::endl;
      } while (!bDone);
    std::cout << "\nEntered input file name: " <<  m_szInFile <<std::endl;

    // open the file
    do{
        m_FileOutput.open (m_szJouFile.c_str(), std::ios::out);
        if (!m_FileOutput){
            std::cout << "Unable to open o/p file: " << m_szJouFile << std::endl;
            m_FileOutput.clear ();
            exit(1);
          }
        else
          bDone = true; // file opened successfully
      } while (!bDone);

    // open the template journal file for writing
    do{
        m_SchemesFile.open (m_szSchFile.c_str(), std::ios::out);
        if (!m_SchemesFile){
            std::cout << "Unable to open o/p file: " << m_szSchFile << std::endl;
            m_SchemesFile.clear ();
            exit(1);
          }
        else
          bDone = true; // file opened successfully
      } while (!bDone);

    do{
        m_PyCubGeomFile.open (m_szPyCubGeom.c_str(), std::ios::out);
        if (!m_PyCubGeomFile){
            std::cout << "Unable to open o/p file: " << m_szPyCubGeom << std::endl;
            m_PyCubGeomFile.clear ();
            exit(1);
          }
        else
          bDone = true; // file opened successfully
      } while (!bDone);

    std::cout<<"\no/p Cubit journal file name: "<< m_szJouFile
            << std::endl;


    //ACIS ENGINE
#ifdef HAVE_ACIS
    //  if(m_szEngine == "acis"){
    m_szGeomFile = m_szFile+".sat";
    //  }
#elif defined(HAVE_OCC)
    //  OCC ENGINE
    //  if (m_szEngine == "occ"){
    m_szGeomFile = m_szFile+".brep";
    m_szGeomFile1 = m_szFile+".sat";
    //  }o
#endif
    std::cout << "\no/p geometry file name: " <<  m_szGeomFile <<std::endl;

    // writing schemes .jou file ends, now write the main journal file.
    // stuff common to both surface and volume
    m_FileOutput << "## This file is created by rgg program in MeshKit ##\n";
    m_FileOutput << "#User needs to specify mesh interval and schemes in this file\n#" << std::endl;

    m_PyCubGeomFile << "## This python script is created by the RGG AssyGen program in MeshKit ##\n";
    m_PyCubGeomFile << "# Here the RGG AssyGen program creates the assembly geometry and mesh\n#" << std::endl;
    m_PyCubGeomFile << "\nimport cubit" << std::endl;


    // write the name faces python function here
    m_PyCubGeomFile << "def name_faces(name, body):\n"
                       "    vector_locs = cubit.get_bounding_box(\"volume\", body.id())\n"
                       "    topno = vector_locs[7] - 1e-2\n"
                       "    botno = vector_locs[6] + 1e-2\n"
                       "    cubit.cmd('group \"g1\" equals surf in vol {0} '.format(body.id()))\n"
                       "    cubit.cmd('group \"g2\" equals surf  in g1 with z_coord  < {0} and z_coord > {1}'.format(topno,botno))\n"
                       "    cubit.cmd('group  \"g3\" subtract g2 from g1')\n"
                       "    cubit.cmd('group \"gtop\" equals surf in g3 with z_coord > {0}'.format(topno) )\n"
                       "    cubit.cmd('group \"gbot\" equals surf in g3 with z_coord < {0}'.format(botno) )\n"
                       "    g2id = cubit.get_id_from_name(\"g2\")\n"
                       "    ssurfs = cubit.get_group_surfaces(g2id)\n"
                       "    side_surfs = len(ssurfs)\n"
                       "    for i in range(0,side_surfs):\n"
                       "      sname = name + \"_side\" + str(i+1)\n"
                       "      cubit.cmd('surf {0} name \"{1}\"'.format( ssurfs[i] , sname )  )\n"
                       "    top_surf = name + \"_top\"\n"
                       "    bot_surf = name + \"_bot\"\n"
                       "    cubit.cmd('surf in gtop name \"{0}\"'.format(top_surf) )\n"
                       "    cubit.cmd('surf in gbot name \"{0}\"'.format(bot_surf) )\n"
                       "    cubit.cmd('delete group g1 g2 g3 gtop gbot')\n" << std::endl;


    m_PyCubGeomFile << "\ncubit.cmd('reset')" << std::endl;


    m_FileOutput << "{include(\"" << m_szSchFile << "\")}" <<std::endl;
    m_FileOutput << "#" << std::endl;
    m_FileOutput << "set logging on file '" << m_szLogFile << "'" <<std::endl;
    m_FileOutput << "Timer Start" << std::endl;
    // import the geometry file
    m_FileOutput << "# Import geometry file " << std::endl;
    //ACIS ENGINE
#ifdef HAVE_ACIS
    m_FileOutput << "import '" << m_szGeomFile <<"'" << std::endl;

#elif defined(HAVE_OCC)
    //  OCC ENGINE
    // Use sat file always as step isn't supported
    m_FileOutput << "import '" << m_szGeomFile1 <<"'" << std::endl;
#endif

    m_FileOutput << "#" << std::endl;

  }


  void AssyGen::ReadCommonInp ()
  // -------------------------------------------------------------------------------------------
  // Function: reads the input file to count the no. of cyl in a pincell, before the actual read
  // Input:    none
  // Output:   none
  // -------------------------------------------------------------------------------------------
  {
    ++com_run_count;
    if(com_run_count > 1){
        //Rewind the reader for common.inp file
        m_FileCommon.clear (std::ios_base::goodbit);
        m_FileCommon.seekg (0L, std::ios::beg);
      }
    CParser Parse1;
    bool found = false;
    std::string card;
    m_nLineNumber = 0;
    std::cout << "Reading from common.inp file." << std::endl;
    for(;;){
        if (!Parse1.ReadNextLine (m_FileCommon, m_nLineNumber, szInputString,
                                  MAXCHARS, szComment))
          IOErrorHandler (INVALIDINPUT);

        if (szInputString.substr(0,10) == "geomengine"){
            found = true;
            std::istringstream szFormatString (szInputString);
            szFormatString >> card >> m_szEngine;
            if( ((strcmp (m_szEngine.c_str(), "acis") != 0) &&
                 (strcmp (m_szEngine.c_str(), "occ") != 0)) || szFormatString.fail())
              IOErrorHandler(EGEOMENGINE);
          }
        // start id for pin number
        if (szInputString.substr(0, 10) == "startpinid") {
            found = true;
            std::istringstream szFormatString(szInputString);
            szFormatString >> card >> m_nStartpinid;
          }
        if (szInputString.substr(0,8) == "meshtype"){
            found = true;
            std::istringstream szFormatString (szInputString);
            szFormatString >> card >> m_szMeshType;
            if( ((strcmp (m_szMeshType.c_str(), "hex") != 0) &&
                 (strcmp (m_szMeshType.c_str(), "tet") != 0)) || szFormatString.fail())
              IOErrorHandler(INVALIDINPUT);
          }
        // info flag
        if (szInputString.substr(0,4) == "info"){
            found = true;
            std::istringstream szFormatString (szInputString);
            szFormatString >> card >> m_szInfo;
            std::cout <<"--------------------------------------------------"<<std::endl;
          }
        // hex block along z
        if (szInputString.substr(0,6) == "hblock"){
            found = true;
            std::istringstream szFormatString (szInputString);
            szFormatString >> card >> m_nHblock >> m_dZstart >> m_dZend;
            std::cout <<"--------------------------------------------------"<<std::endl;
          }
        // Hex or Rect geometry type
        if (szInputString.substr(0,12) == "geometrytype"){
            found = true;
            std::istringstream szFormatString (szInputString);
            szFormatString >> card >> m_szGeomType;
            if( ((strcmp (m_szGeomType.c_str(), "hexagonal") != 0) &&
                 (strcmp (m_szGeomType.c_str(), "rectangular") != 0)) || szFormatString.fail())
              IOErrorHandler(EGEOMTYPE);

            // set the number of sides in the geometry
            if(m_szGeomType == "hexagonal")
              m_nSides = 6;
            else  if(m_szGeomType == "rectangular")
              m_nSides = 4;
          }
        // Default if volume, set geometry type to surface for 2D assemblies
        if (szInputString.substr(0,8) == "geometry"){
            found = true;
            std::string outfile;
            std::istringstream szFormatString (szInputString);
            szFormatString >> card >> outfile;
            if(strcmp (outfile.c_str(), "surface") == 0 || szFormatString.fail())
              m_nPlanar=1;
          }
        // 'yes' or 'no' for creating sidesets
        if (szInputString.substr(0,13) == "createsideset"){
            found = true;
            std::istringstream szFormatString (szInputString);
            szFormatString >> card >> m_szSideset;
            std::cout <<"--------------------------------------------------"<<std::endl;
          }
        // Create specified number of files with varying material ids
        if (szInputString.substr(0,11) == "createfiles"){
            found = true;
            std::istringstream szFormatString (szInputString);
            szFormatString >> card >> m_nAssyGenInputFiles;
            std::cout <<"--------------------------------------------------"<<std::endl;
          }
        // Create specified number of files with varying material ids
        if (szInputString.substr(0,11) == "save_exodus"){
            found = true;
            save_exodus = true;
            std::cout <<"--------------------------------------------------"<<std::endl;
          }
        // specify a merge tolerance value for cubit journal file
        if (szInputString.substr(0,14) == "mergetolerance"){
            found = true;
            std::istringstream szFormatString (szInputString);
            szFormatString >> card >> m_dMergeTol;
            std::cout <<"--------------------------------------------------"<<std::endl;
          }
        // Handle mesh size inputs
        if (szInputString.substr(0,14) == "radialmeshsize"){
            found = true;
            std::istringstream szFormatString (szInputString);
            szFormatString >> card >> m_dRadialSize;
            if(m_dRadialSize < 0 || szFormatString.fail())
              IOErrorHandler(ENEGATIVE);
            std::cout <<"--------------------------------------------------"<<std::endl;

          }
        // Handle mesh size inputs
        if (szInputString.substr(0,11) == "tetmeshsize"){
            found = true;
            std::istringstream szFormatString (szInputString);
            szFormatString >> card >> m_dTetMeshSize;
            if(m_dTetMeshSize < 0 || szFormatString.fail())
              IOErrorHandler(ENEGATIVE);
            std::cout <<"--------------------------------------------------"<<std::endl;

          }
        // Handle mesh size inputs
        if (szInputString.substr(0,13) == "axialmeshsize"){
            found = true;
            if(com_run_count > 1){
                std::istringstream szFormatString (szInputString);
                szFormatString >> card;
                m_dAxialSize.SetSize(m_nDuct);
                int num_ams_specified = std::distance(std::istream_iterator<std::string>(szFormatString),
                                                      std::istream_iterator<std::string>());
                std::istringstream szFormatStringAgain (szInputString);
                szFormatStringAgain >> card;
                for (int p = 1; p <= m_nDuct; p++){
                    if(p <= num_ams_specified)
                      szFormatStringAgain >> m_dAxialSize(p);
                    else
                      m_dAxialSize(p) = m_dAxialSize(num_ams_specified);
                    if(m_dAxialSize(p) < 0)
                      IOErrorHandler(ENEGATIVE);
                  }
                std::cout <<"--------------------------------------------------"<<std::endl;
              }
          }
        // edge interval
        if (szInputString.substr(0, 12) == "edgeinterval") {
            found = true;
            std::istringstream szFormatString(szInputString);
            szFormatString >> card >> m_edgeInterval;
          }
        // mesh scheme - hole or pave
        if (szInputString.substr(0, 10) == "meshscheme") {
            std::istringstream szFormatString(szInputString);
            szFormatString >> card >> m_szMeshScheme;
          }
        // breaking condition
        if(szInputString.substr(0,3) == "end" || m_nLineNumber == MAXLINES){
            found = true;
            break;
          }
        if (found == false){
            std::cout << "Cannot specify: " << szInputString << " in common.inp files" << std::endl;
          }
      }

  }

  void AssyGen::ReadInputPhase1 ()
  // -------------------------------------------------------------------------------------------
  // Function: reads the input file to count the no. of cyl in a pincell, before the actual read
  // Input:    none
  // Output:   none
  // -------------------------------------------------------------------------------------------
  {
    std::cout << "Reading from AssyGen input file." << std::endl;
    CParser Parse;  bool bDone = false;
    int nCyl =0, nCellMat=0, nInputLines=0;
    std::string card, szVolId, szVolAlias;
    m_nLineNumber = 0;
    // count the total number of cylinder commands in each pincell
    for(;;){
        if (!Parse.ReadNextLine (m_FileInput, m_nLineNumber, szInputString,
                                 MAXCHARS, szComment))
          IOErrorHandler (INVALIDINPUT);

        if (szInputString.substr(0,10) == "geomengine"){
            std::istringstream szFormatString (szInputString);
            szFormatString >> card >> m_szEngine;
            if( ((strcmp (m_szEngine.c_str(), "acis") != 0) &&
                 (strcmp (m_szEngine.c_str(), "occ") != 0)) || szFormatString.fail())
              IOErrorHandler(EGEOMENGINE);
          }
        // Read material data
        if ((szInputString.substr(0,9) == "materials") && (szInputString.substr(0,19) != "materialset_startid")){

            std::istringstream szFormatString (szInputString);
            szFormatString >> card >> m_nAssemblyMat;
            if(szFormatString.fail())
              IOErrorHandler(INVALIDINPUT);
            m_szAssmMat.SetSize(m_nAssemblyMat); m_szAssmMatAlias.SetSize(m_nAssemblyMat);
            for (int j=1; j<=m_nAssemblyMat; j++){
                szFormatString >> m_szAssmMat(j) >> m_szAssmMatAlias(j);
                if( (strcmp (m_szAssmMat(j).c_str(), "") == 0) ||
                    (strcmp (m_szAssmMatAlias(j).c_str(), "") == 0)){
                    IOErrorHandler(EMAT);
                  }
                // checking if & inserted at the end of the material by mistake
                if (j == m_nAssemblyMat){
                    std::string dummy = "";
                    szFormatString >> dummy;
                    if (strcmp (dummy.c_str(), "") != 0)
                      IOErrorHandler(EMAT);
                  }
              }
          }
        // Read material data
        if (szInputString.substr(0,11) == "blmaterials"){

            std::istringstream szFormatString (szInputString);
            szFormatString >> card >> m_nBLAssemblyMat;
            if(szFormatString.fail())
              IOErrorHandler(INVALIDINPUT);
            m_szBLAssmMat.SetSize(m_nBLAssemblyMat); m_dBLMatBias.SetSize(m_nBLAssemblyMat);
            m_nBLMatIntervals.SetSize(m_nBLAssemblyMat);
            for (int j=1; j<=m_nBLAssemblyMat; j++){
                szFormatString >> m_szBLAssmMat(j) >> m_dBLMatBias(j) >> m_nBLMatIntervals(j);
                if( (strcmp (m_szBLAssmMat(j).c_str(), "") == 0) ||
                    (m_nBLMatIntervals(j) < 0) ){
                    IOErrorHandler(EMAT);
                  }
                // checking if & inserted at the end of the material by mistake
                if (j == m_nBLAssemblyMat){
                    std::string dummy = "";
                    szFormatString >> dummy;
                    if (strcmp (dummy.c_str(), "") != 0)
                      IOErrorHandler(EMAT);
                  }
              }
          }

        // start id for pin number
        if (szInputString.substr(0, 10) == "startpinid") {
            std::istringstream szFormatString(szInputString);
            szFormatString >> card >> m_nStartpinid;
          }
        if (szInputString.substr(0,8) == "meshtype"){
            std::istringstream szFormatString (szInputString);
            szFormatString >> card >> m_szMeshType;
            if( ((strcmp (m_szMeshType.c_str(), "hex") != 0) &&
                 (strcmp (m_szMeshType.c_str(), "tet") != 0)) || szFormatString.fail())
              IOErrorHandler(INVALIDINPUT);
          }
        if (szInputString.substr(0,4) == "duct" || szInputString.substr(0,10) == "dimensions"){
            ++m_nDuct;
          }
        if (szInputString.substr(0,8) == "pincells"){
            std::istringstream szFormatString (szInputString);
            szFormatString >> card >> m_nPincells;
            if(m_nPincells>0)
              m_Pincell.SetSize(m_nPincells);
            else if(m_nPincells ==0)
              m_Pincell.SetSize(1); // assume for using dummy pincell

            // count the number of cylinder lines for each pincell
            for (int i=1; i<=m_nPincells; i++){
                // read the no. of input lines first pincell
                if (!Parse.ReadNextLine (m_FileInput, m_nLineNumber, szInputString,
                                         MAXCHARS, szComment))
                  IOErrorHandler (INVALIDINPUT);
                std::istringstream szFormatString1 (szInputString);
                szFormatString1 >> szVolId >> szVolAlias >> nInputLines;
                if(szFormatString1.fail())
                  IOErrorHandler(INVALIDINPUT);
                // loop thru the input lines of each pincell
                for(int l=1; l<=nInputLines; l++){
                    if (!Parse.ReadNextLine (m_FileInput, m_nLineNumber, szInputString,
                                             MAXCHARS, szComment))
                      IOErrorHandler (INVALIDINPUT);
                    if (szInputString.substr(0,8) == "cylinder" || szInputString.substr(0,7) == "frustum"){
                        ++nCyl;
                      }
                    if (szInputString.substr(0,12) == "cellmaterial"){
                        ++nCellMat;
                      }
                  }

                // set the sizes
                if(nCyl>0){
                    if  (nCellMat!=0){
                        m_Pincell(i).SetCellMatSize(nCyl);
                      }
                    m_Pincell(i).SetNumCyl(nCyl);
                  }
                else if(nCyl ==0){
                    if(nInputLines >0)
                      m_Pincell(i).SetCellMatSize(nCellMat);
                  }
                nCyl = 0;
                nCellMat = 0;
              }
          }
        // info flag
        if (szInputString.substr(0,4) == "info"){
            std::istringstream szFormatString (szInputString);
            szFormatString >> card >> m_szInfo;
            std::cout <<"--------------------------------------------------"<<std::endl;
          }
        // merge
        if (szInputString.substr(0,7) == "imprint"){
            m_bimprint = true;
            std::cout <<"--------------------------------------------------"<<std::endl;
          }
        // imprint
        if (szInputString.substr(0,5) == "merge"){
            m_bmerge = true;
            std::cout <<"--------------------------------------------------"<<std::endl;
          }
        // info flag
        if (szInputString.substr(0,6) == "smooth"){
            std::istringstream szFormatString (szInputString);
            szFormatString >> card >> m_szSmooth;
            std::cout <<"--------------------------------------------------"<<std::endl;
          }
        // mesh scheme - hole or pave
        if (szInputString.substr(0, 10) == "meshscheme") {
            std::istringstream szFormatString(szInputString);
            szFormatString >> card >> m_szMeshScheme;
          }
        // hex block along z
        if (szInputString.substr(0,6) == "hblock"){
            std::istringstream szFormatString (szInputString);
            szFormatString >> card >> m_nHblock >> m_dZstart >> m_dZend;
            std::cout <<"--------------------------------------------------"<<std::endl;
          }
        // breaking condition
        if(szInputString.substr(0,3) == "end" || m_nLineNumber == MAXLINES){
            break;
          }
      }

    if(strcmp(m_szInfo.c_str(),"on") == 0){
        std::cout  << m_szAssmInfo <<std::endl;
        do{
            m_AssmInfo.open (m_szAssmInfo.c_str(), std::ios::out);
            if (!m_AssmInfo){
                std::cout << "Unable to open o/p file: " << m_szAssmInfo << std::endl;
                m_AssmInfo.clear ();
                exit(1);
              }
            else
              bDone = true; // file opened successfully
          } while (!bDone);

        // write header for info file
        m_AssmInfo <<"pincell"<<  " \t" <<
                     "m" << " \t" << "n" << " \t" << "dX" << " \t" <<
                     "dY" << " \t" << "dZ"  << std::endl;
      }
    // set the size of cp_inpins matrix
    //  cp_inpins.resize(m_nDuct);
    for (int j=0; j<m_nDuct ; j++)
      cp_inpins.push_back(std::vector<iBase_EntityHandle>());

  }

  void AssyGen::ReadPinCellData (int i)
  //---------------------------------------------------------------------------
  //Function: reading pincell i from file and storing the data
  //Input:    none
  //Output:   none
  //---------------------------------------------------------------------------
  {
    CParser Parse;
    std::string card, szVolId, szVolAlias, szIFlag;
    int nInputLines, nMaterials, nCyl = 0, nRadii=0, nCellMat=0;
    double dLZ=0.0, dFlatF=0.0, dPX=0.0, dPY=0.0, dPZ=0.0;
    CVector <std::string> szVMatName, szVMatAlias, szVCylMat, szVCellMat;
    CVector<double> dVCoor(2), dVCylRadii, dVCylZPos, dZVStart, dZVEnd;

    //loop over input lines
    if (m_szGeomType == "rectangular"){

        std::cout << "\ngetting volume id";
        if (!Parse.ReadNextLine (m_FileInput, m_nLineNumber, szInputString,
                                 MAXCHARS, szComment))
          IOErrorHandler (INVALIDINPUT);
        std::istringstream szFormatString (szInputString);
        szFormatString >> szVolId >> szVolAlias >> nInputLines >> szIFlag;

        // error checking
        if( (strcmp (szVolAlias.c_str(), "") == 0) ||
            (strcmp (szVolId.c_str(), "") == 0))
          IOErrorHandler(EPIN);
        if( nInputLines < 0 )
          IOErrorHandler(ENEGATIVE);

        m_Pincell(i).SetLineOne (szVolId, szVolAlias, nInputLines);
        if(szIFlag == "intersect"){
            m_Pincell(i).SetIntersectFlag(1);
          }
        else{
            m_Pincell(i).SetIntersectFlag(0);
          }
        for(int l=1; l<=nInputLines; l++){
            if (!Parse.ReadNextLine (m_FileInput, m_nLineNumber, szInputString,
                                     MAXCHARS, szComment))
              IOErrorHandler (INVALIDINPUT);
            if (szInputString.substr(0,5) == "pitch"){

                std::istringstream szFormatString (szInputString);
                std::cout << "\ngetting pitch data";
                szFormatString >> card >> dPX >> dPY >> dPZ;

                if( dPX < 0 || dPY < 0 || dPZ < 0 || szFormatString.fail())
                  IOErrorHandler(ENEGATIVE);
                m_Pincell(i).SetPitch (dPX, dPY, dPZ);
              }
            if (szInputString.substr(0,9) == "materials"){

                std::istringstream szFormatString (szInputString);
                szFormatString >> card >> nMaterials;
                if(szFormatString.fail())
                  IOErrorHandler(INVALIDINPUT);

                //setting local arrays
                szVMatName.SetSize(nMaterials);
                szVMatAlias.SetSize(nMaterials);

                //set class variable sizes
                m_Pincell(i).SetMatArray(nMaterials);
                std::cout << "\ngetting material data";
                for(int j=1; j<= nMaterials; j++){
                    szFormatString >> szVMatName(j) >> szVMatAlias(j);
                    if(szFormatString.fail())
                      IOErrorHandler(INVALIDINPUT);
                  }
                m_Pincell(i).SetMat(szVMatName, szVMatAlias);
              }
            if (szInputString.substr(0,8) == "cylinder"){

                ++nCyl;
                std::cout << "\ngetting cylinder data";
                std::istringstream szFormatString (szInputString);
                szFormatString >> card >> nRadii >> dVCoor(1) >> dVCoor(2);
                if(szFormatString.fail())
                  IOErrorHandler(INVALIDINPUT);
                m_Pincell(i).SetCylSizes(nCyl, nRadii);
                m_Pincell(i).SetCylPos(nCyl, dVCoor);
                m_Pincell(i).SetCellType(nCyl, 0);

                //set local array
                dVCylRadii.SetSize(2*nRadii);
                szVCylMat.SetSize(nRadii);
                dVCylZPos.SetSize(2);
                m_Pincell(i).SetCylSizes(nCyl, nRadii);

                // reading ZCoords
                for(int k=1; k<=2; k++){
                    szFormatString >> dVCylZPos(k);
                    if(szFormatString.fail())
                      IOErrorHandler(INVALIDINPUT);
                  }
                m_Pincell(i).SetCylZPos(nCyl, dVCylZPos);

                // reading Radii
                for(int l=1; l<= nRadii; l++){
                    szFormatString >> dVCylRadii(l);
                    if( dVCylRadii(l) < 0  || szFormatString.fail())
                      IOErrorHandler(ENEGATIVE);
                  }
                m_Pincell(i).SetCylRadii(nCyl, dVCylRadii);

                // reading Material alias
                for(int m=1; m<= nRadii; m++){
                    szFormatString >> szVCylMat(m);
                    if(strcmp (szVCylMat(m).c_str(), "") == 0 || szFormatString.fail())
                      IOErrorHandler(EALIAS);
                    //                  // setting stuff for hole scheme determination for meshing
                    //                  if (m > 1 && m_szMeshScheme == "hole" && m_nBLAssemblyMat == 0){
                    //                      // find material name for this alias
                    //                      for (int ll=1; ll<= m_nAssemblyMat; ll++){
                    //                          if(szVCylMat(m) == m_szAssmMatAlias(ll))
                    //                            m_FileOutput << "group 'hole_surfaces' add surface name '"<< m_szAssmMat(ll)  << "_top'" << std::endl;
                    //                        }
                    //                    }
                    m_Pincell(i).SetCylMat(nCyl, szVCylMat);
                  }
              }
            if (szInputString.substr(0,7) == "frustum"){

                ++nCyl;
                std::cout << "\ngetting frustum data";
                std::istringstream szFormatString (szInputString);
                szFormatString >> card >> nRadii >> dVCoor(1) >> dVCoor(2);
                if(szFormatString.fail())
                  IOErrorHandler(INVALIDINPUT);
                m_Pincell(i).SetCylSizes(nCyl, nRadii);
                m_Pincell(i).SetCylPos(nCyl, dVCoor);
                m_Pincell(i).SetCellType(nCyl, 1);

                //set local array
                dVCylRadii.SetSize(2*nRadii);
                szVCylMat.SetSize(nRadii);
                dVCylZPos.SetSize(2);
                m_Pincell(i).SetCylSizes(nCyl, nRadii);

                // reading ZCoords
                for(int k=1; k<=2; k++){
                    szFormatString >> dVCylZPos(k);
                    if(szFormatString.fail())
                      IOErrorHandler(INVALIDINPUT);
                  }
                m_Pincell(i).SetCylZPos(nCyl, dVCylZPos);

                // reading Radii
                for(int l=1; l<= 2*nRadii; l++){
                    szFormatString >> dVCylRadii(l);
                    if( dVCylRadii(l) < 0  || szFormatString.fail())
                      IOErrorHandler(ENEGATIVE);
                  }
                m_Pincell(i).SetCylRadii(nCyl, dVCylRadii);

                // reading Material alias
                for(int m=1; m<= nRadii; m++){
                    szFormatString >> szVCylMat(m);
                    if(strcmp (szVCylMat(m).c_str(), "") == 0 || szFormatString.fail())
                      IOErrorHandler(EALIAS);
                    // setting stuff for hole scheme determination for meshing
                    if (m > 2 && m_szMeshScheme == "hole"){
                        // find material name for this alias
                        for (int ll=1; ll<= m_nAssemblyMat; ll++){
                            if(szVCylMat(m) == m_szAssmMatAlias(ll))
                              m_FileOutput << "group 'hole_surfaces' add surface name '"<< m_szAssmMat(ll)  << "_top'" << std::endl;
                          }
                      }
                  }
                m_Pincell(i).SetCylMat(nCyl, szVCylMat);
              }
            if (szInputString.substr(0,12) == "cellmaterial"){

                std::cout << "\ngetting cell material data\n";
                std::istringstream szFormatString (szInputString);
                szFormatString >> card;

                //set local arrays
                m_Pincell(i).GetCellMatSize(nCellMat); // since size of cell material already set equal to number of cylinders
                dZVStart.SetSize(nCellMat);
                dZVEnd.SetSize(nCellMat);
                szVCellMat.SetSize(nCellMat);

                for(int k=1; k<=nCellMat; k++){
                    szFormatString >> dZVStart(k)>> dZVEnd(k) >> szVCellMat(k);
                    if(strcmp (szVCellMat(k).c_str(), "") == 0 || szFormatString.fail())
                      IOErrorHandler(EALIAS);
                  }
                m_Pincell(i).SetCellMat(dZVStart, dZVEnd, szVCellMat);
              }
          }
      }//if rectangular ends

    if (m_szGeomType == "hexagonal"){

        std::cout << "\ngetting volume id";
        if (!Parse.ReadNextLine (m_FileInput, m_nLineNumber, szInputString,
                                 MAXCHARS, szComment))
          IOErrorHandler (INVALIDINPUT);
        std::istringstream szFormatString (szInputString);
        szFormatString >> szVolId >> szVolAlias >> nInputLines >> szIFlag;

        // error checking
        if( (strcmp (szVolAlias.c_str(), "") == 0) ||
            (strcmp (szVolId.c_str(), "") == 0))
          IOErrorHandler(EPIN);
        if( nInputLines < 0 )
          IOErrorHandler(ENEGATIVE);

        m_Pincell(i).SetLineOne (szVolId, szVolAlias, nInputLines);
        if(szIFlag == "intersect"){
            m_Pincell(i).SetIntersectFlag(1);
          }
        else{
            m_Pincell(i).SetIntersectFlag(0);
          }
        for(int l=1; l<=nInputLines; l++){
            if (!Parse.ReadNextLine (m_FileInput, m_nLineNumber, szInputString,
                                     MAXCHARS, szComment))
              IOErrorHandler (INVALIDINPUT);
            if (szInputString.substr(0,5) == "pitch"){

                std::istringstream szFormatString (szInputString);
                std::cout << "\ngetting pitch data";
                szFormatString >> card >> dFlatF >> dLZ;
                if( dFlatF < 0 || dLZ < 0  || szFormatString.fail())
                  IOErrorHandler(ENEGATIVE);
                m_Pincell(i).SetPitch (dFlatF, dLZ);
              }
            if (szInputString.substr(0,9) == "materials"){

                std::istringstream szFormatString (szInputString);
                szFormatString >> card >> nMaterials;
                if(szFormatString.fail())
                  IOErrorHandler(INVALIDINPUT);
                //setting local arrays
                szVMatName.SetSize(nMaterials);
                szVMatAlias.SetSize(nMaterials);

                //set class variable sizes
                m_Pincell(i).SetMatArray(nMaterials);
                std::cout << "\ngetting material data";
                for(int j=1; j<= nMaterials; j++){
                    szFormatString >> szVMatName(j) >> szVMatAlias(j);
                    if(szFormatString.fail())
                      IOErrorHandler(INVALIDINPUT);
                  }
                m_Pincell(i).SetMat(szVMatName, szVMatAlias);
              }
            if (szInputString.substr(0,8) == "cylinder"){

                ++nCyl;
                std::cout << "\ngetting cylinder data";
                std::istringstream szFormatString (szInputString);
                szFormatString >> card >> nRadii >> dVCoor(1) >> dVCoor(2);
                if(szFormatString.fail())
                  IOErrorHandler(INVALIDINPUT);
                m_Pincell(i).SetCylSizes(nCyl, nRadii);
                m_Pincell(i).SetCylPos(nCyl, dVCoor);
                m_Pincell(i).SetCellType(nCyl, 0);

                //set local array
                dVCylRadii.SetSize(2*nRadii);
                szVCylMat.SetSize(nRadii);
                dVCylZPos.SetSize(2);
                //
                m_Pincell(i).SetCylSizes(nCyl, nRadii);

                // reading ZCoords - max and min 2 always
                for(int k=1; k<=2; k++)
                  szFormatString >> dVCylZPos(k);
                m_Pincell(i).SetCylZPos(nCyl, dVCylZPos);

                // reading Radii
                for(int l=1; l<= nRadii; l++){
                    szFormatString >> dVCylRadii(l);
                    if( dVCylRadii(l) < 0 || szFormatString.fail())
                      IOErrorHandler(ENEGATIVE);
                  }
                m_Pincell(i).SetCylRadii(nCyl, dVCylRadii);

                // reading Material alias
                for(int m=1; m<= nRadii; m++){
                    szFormatString >> szVCylMat(m);
                    if(strcmp (szVCylMat(m).c_str(), "") == 0 || szFormatString.fail())
                      IOErrorHandler(EALIAS);
                    // setting stuff for hole scheme determination for meshing
                    if (m > 1 && m_szMeshScheme == "hole" && m_nBLAssemblyMat == 0){
                        // find material name for this alias
                        for (int ll=1; ll<= m_nAssemblyMat; ll++){
                            //   if(szVCylMat(m) == m_szAssmMatAlias(ll))
                            if(strcmp (m_szAssmMatAlias(ll).c_str(), szVCylMat(m).c_str()) == 0)
                              m_FileOutput << "group 'hole_surfaces' add surface name '"<< m_szAssmMat(ll)  << "_top'" << std::endl;
                          }
                      }
                  }

                m_Pincell(i).SetCylMat(nCyl, szVCylMat);
              }
            if (szInputString.substr(0,7) == "frustum"){

                ++nCyl;
                std::cout << "\ngetting frustum data";
                std::istringstream szFormatString (szInputString);
                szFormatString >> card >> nRadii >> dVCoor(1) >> dVCoor(2);
                if(szFormatString.fail())
                  IOErrorHandler(INVALIDINPUT);
                m_Pincell(i).SetCylSizes(nCyl, nRadii);
                m_Pincell(i).SetCylPos(nCyl, dVCoor);
                m_Pincell(i).SetCellType(nCyl, 1);

                //set local array
                dVCylRadii.SetSize(2*nRadii);
                szVCylMat.SetSize(nRadii);
                dVCylZPos.SetSize(2);
                m_Pincell(i).SetCylSizes(nCyl, nRadii);

                // reading ZCoords
                for(int k=1; k<=2; k++){
                    szFormatString >> dVCylZPos(k);
                    if(szFormatString.fail())
                      IOErrorHandler(INVALIDINPUT);
                  }
                m_Pincell(i).SetCylZPos(nCyl, dVCylZPos);

                // reading Radii
                for(int l=1; l<= 2*nRadii; l++){
                    szFormatString >> dVCylRadii(l);
                    if( dVCylRadii(l) < 0  || szFormatString.fail())
                      IOErrorHandler(ENEGATIVE);
                  }
                m_Pincell(i).SetCylRadii(nCyl, dVCylRadii);

                // reading Material alias
                for(int m=1; m<= nRadii; m++){
                    szFormatString >> szVCylMat(m);
                    if(strcmp (szVCylMat(m).c_str(), "") == 0 || szFormatString.fail())
                      IOErrorHandler(EALIAS);
                    // setting stuff for hole scheme determination for meshing
                    if (m > 2 && m_szMeshScheme == "hole"){
                        // find material name for this alias
                        for (int ll=1; ll<= m_nAssemblyMat; ll++){
                            if(szVCylMat(m) == m_szAssmMatAlias(ll))
                              m_FileOutput << "group 'hole_surfaces' add surface name '"<< m_szAssmMat(ll)  << "_top'" << std::endl;
                          }
                      }
                  }
                m_Pincell(i).SetCylMat(nCyl, szVCylMat);
              }
            if (szInputString.substr(0,12) == "cellmaterial"){

                std::cout << "\ngetting cell material data";
                std::istringstream szFormatString (szInputString);
                szFormatString >> card;

                //set local arrays
                m_Pincell(i).GetCellMatSize(nCellMat); // since size of cell material already set equal to number of cylinders
                dZVStart.SetSize(nCellMat);
                dZVEnd.SetSize(nCellMat);
                szVCellMat.SetSize(nCellMat);

                for(int k=1; k<=nCellMat; k++){
                    szFormatString >> dZVStart(k)>> dZVEnd(k) >> szVCellMat(k);
                    if(strcmp (szVCellMat(k).c_str(), "") == 0 || szFormatString.fail())
                      IOErrorHandler(EALIAS);
                  }
                m_Pincell(i).SetCellMat(dZVStart, dZVEnd, szVCellMat);
              }
          }
      }// if hexagonal end

  }


  void AssyGen::ReadAndCreate()
  //---------------------------------------------------------------------------
  //Function: reads the input file and creates assembly
  //Input:    none
  //Output:   none
  //---------------------------------------------------------------------------
  {
    //Rewind the input file
    m_FileInput.clear (std::ios_base::goodbit);
    m_FileInput.seekg (0L, std::ios::beg);
    m_nLineNumber = 0;
    CParser Parse;
    std::string card;

    // start reading the input file break when encounter end
    for(;;){
        if (!Parse.ReadNextLine (m_FileInput, m_nLineNumber, szInputString,
                                 MAXCHARS, szComment))
          IOErrorHandler (INVALIDINPUT);
        if (szInputString.substr(0,12) == "geometrytype"){
            std::istringstream szFormatString (szInputString);
            szFormatString >> card >> m_szGeomType;
            if( ((strcmp (m_szGeomType.c_str(), "hexagonal") != 0) &&
                 (strcmp (m_szGeomType.c_str(), "rectangular") != 0)) || szFormatString.fail())
              IOErrorHandler(EGEOMTYPE);

            // set the number of sides in the geometry
            if(m_szGeomType == "hexagonal")
              m_nSides = 6;
            else  if(m_szGeomType == "rectangular")
              m_nSides = 4;
          }

        if (szInputString.substr(0,8) == "geometry"){
            std::string outfile;
            std::istringstream szFormatString (szInputString);
            szFormatString >> card >> outfile;
            if(strcmp (outfile.c_str(), "surface") == 0 || szFormatString.fail())
              m_nPlanar=1;
          }
        if( (szInputString.substr(0,10) == "dimensions") ||
            (szInputString.substr(0,4) == "duct") ){

            ++m_nDuctNum;
            std::cout << "getting assembly dimensions " << m_nDuctNum << std::endl;

            if(m_szGeomType =="hexagonal"){
                std::istringstream szFormatString (szInputString);

                if(m_nDuctNum == 1){
                    m_dMXYAssm.SetSize(m_nDuct, 2); m_dMZAssm.SetSize(m_nDuct, 2);
                  }
                szFormatString >> card >> m_nDimensions
                    >> m_dMXYAssm(m_nDuctNum, 1) >> m_dMXYAssm(m_nDuctNum, 2)
                                                 >> m_dMZAssm(m_nDuctNum, 1) >> m_dMZAssm(m_nDuctNum, 2);
                if(m_nDuctNum == 1){
                    m_dMAssmPitch.SetSize(m_nDuct, m_nDimensions); m_szMMAlias.SetSize(m_nDuct, m_nDimensions);

                    assms.resize(m_nDimensions*m_nDuct); // setup while reading the problem size
                    // Declaration for python script
                    m_PyCubGeomFile << "assms = range(" << m_nDimensions*m_nDuct << ")\ncp_inpins  = []\n" << std::endl;
                  }

                for (int i=1; i<=m_nDimensions; i++){
                    szFormatString >> m_dMAssmPitch(m_nDuctNum, i);
                    if( m_dMAssmPitch(m_nDuctNum, i) < 0 )
                      IOErrorHandler(ENEGATIVE);
                  }

                for (int i=1; i<=m_nDimensions; i++){
                    szFormatString >> m_szMMAlias(m_nDuctNum, i);
                    if(strcmp (m_szMMAlias(m_nDuctNum, i).c_str(), "") == 0)
                      IOErrorHandler(EALIAS);
                    m_szDuctMats.push_back(m_szMMAlias(m_nDuctNum, i));
                    // this is the innermost duct add to group for journaling later
                    if (i==1){
                        // find material name for this alias
                        for (int ll=1; ll<= m_nAssemblyMat; ll++){
                            //   if(szVCylMat(m) == m_szAssmMatAlias(ll))
                            if(strcmp (m_szAssmMatAlias(ll).c_str(),  m_szMMAlias(m_nDuctNum, i).c_str()) == 0)
                              m_FileOutput << "group 'innerduct' add surface name '"<< m_szAssmMat(ll)  << "_top'" << std::endl;
                          }
                      }
                    // setting stuff for hole scheme determination for meshing
                    if (i > 1 && m_szMeshScheme == "hole"){
                        // find material name for this alias
                        for (int ll=1; ll<= m_nAssemblyMat; ll++){
                            //   if(szVCylMat(m) == m_szAssmMatAlias(ll))
                            if(strcmp (m_szAssmMatAlias(ll).c_str(),  m_szMMAlias(m_nDuctNum, i).c_str()) == 0)
                              m_FileOutput << "group 'hole_surfaces' add surface name '"<< m_szAssmMat(ll)  << "_top'" << std::endl;
                          }
                      }
                  }
              }
            if(m_szGeomType =="rectangular"){
                std::istringstream szFormatString (szInputString);
                if(m_nDuctNum == 1){
                    m_dMXYAssm.SetSize(m_nDuct, 2);
                    m_dMZAssm.SetSize(m_nDuct, 2);
                  }
                szFormatString >> card >> m_nDimensions
                    >> m_dMXYAssm(m_nDuctNum, 1) >> m_dMXYAssm(m_nDuctNum, 2)
                                                 >> m_dMZAssm(m_nDuctNum, 1) >> m_dMZAssm(m_nDuctNum, 2);
                if (szFormatString.fail())
                  IOErrorHandler(INVALIDINPUT);
                if(m_nDuctNum == 1){
                    m_dMAssmPitchX.SetSize(m_nDuct, m_nDimensions);
                    m_dMAssmPitchY.SetSize(m_nDuct, m_nDimensions);
                    m_szMMAlias.SetSize(m_nDuct, m_nDimensions);
                    assms.resize(m_nDimensions*m_nDuct);
                    // Declaration for python script
                    m_PyCubGeomFile << "assms = range(" << m_nDimensions*m_nDuct << ")\ncp_inpins  = []\n" << std::endl;
                  }
                for (int i=1; i<=m_nDimensions; i++){
                    szFormatString >> m_dMAssmPitchX(m_nDuctNum, i) >> m_dMAssmPitchY(m_nDuctNum, i);
                    if( m_dMAssmPitchX(m_nDuctNum, i) < 0 || m_dMAssmPitchY(m_nDuctNum, i) < 0 || szFormatString.fail())
                      IOErrorHandler(ENEGATIVE);
                  }

                for (int i=1; i<=m_nDimensions; i++){
                    szFormatString >> m_szMMAlias(m_nDuctNum, i);
                    if(strcmp (m_szMMAlias(m_nDuctNum, i).c_str(), "") == 0 || szFormatString.fail())
                      IOErrorHandler(EALIAS);
                    m_szDuctMats.push_back(m_szMMAlias(m_nDuctNum, i));
                    // this is the innermost duct
                    if (i==1){
                        // find material name for this alias
                        for (int ll=1; ll<= m_nAssemblyMat; ll++){
                            //   if(szVCylMat(m) == m_szAssmMatAlias(ll))
                            if(strcmp (m_szAssmMatAlias(ll).c_str(),  m_szMMAlias(m_nDuctNum, i).c_str()) == 0)
                              m_FileOutput << "group 'innerduct' add surface name '"<< m_szAssmMat(ll)  << "_top'" << std::endl;
                          }
                      }
                    // setting stuff for hole scheme determination for meshing
                    if (i > 1 && m_szMeshScheme == "hole"){
                        // find material name for this alias
                        for (int ll=1; ll<= m_nAssemblyMat; ll++){
                            //   if(szVCylMat(m) == m_szAssmMatAlias(ll))
                            if(strcmp (m_szAssmMatAlias(ll).c_str(),  m_szMMAlias(m_nDuctNum, i).c_str()) == 0)
                              m_FileOutput << "group 'hole_surfaces' add surface name '"<< m_szAssmMat(ll)  << "_top'" << std::endl;
                          }
                      }
                  }
              }
          }
        if (szInputString.substr(0,8) == "pincells"){
            std::istringstream szFormatString (szInputString);

            szFormatString >> card >> m_nPincells >> m_dPitch;
            if(m_nPincells < 0)
              IOErrorHandler(ENEGATIVE);

            // this is an option if a user wants to specify pitch here
            double dTotalHeight = 0.0;

            //get the number of cylinder in each pincell
            int nTemp = 1;
            if(m_nDimensions > 0){
                dTotalHeight = m_dMZAssm(nTemp, 2)-m_dMZAssm(nTemp, 1);
              }
            else{
                dTotalHeight = 0; // nothing specified only pincells in the model
              }

            // loop thro' the pincells and read/store pincell data
            for (int i=1; i<=m_nPincells; i++){

                // set pitch if specified in pincell card
                if(m_dPitch > 0.0)
                  m_Pincell(i).SetPitch(m_dPitch, dTotalHeight);

                ReadPinCellData(i);
                //ERRORR("Error in ReadPinCellData", err);
                std::cout << "\nread pincell " << i << std::endl;
              }
          }
        if (szInputString.substr(0,8) == "assembly"){
            if(m_szGeomType =="hexagonal"){
                Create_HexAssm(szInputString);
                //ERRORR("Error in Create_HexAssm", err);
              }
            if(m_szGeomType =="rectangular"){
                Create_CartAssm(szInputString);
                //ERRORR("Error in Create_CartAssm", err);
              }
            if (m_nJouFlag == 0){
                CreateOuterCovering();
                //ERRORR("Error in CreateOuterCovering", err);

                // subtract pins before save
                if(m_nDuct > 0){
                    Subtract_Pins();
                    clock_t s_subtract = clock();
                    std::cout << "## Subract Pins CPU time used := " << (double) (clock() - s_subtract)/CLOCKS_PER_SEC
                              << " seconds" << std::endl;
                  }
                if(m_nPlanar ==1){
                    Create2DSurf();
                    //ERRORR("Error in Create2DSurf", err);
                  }
              }
          }

        // section the assembly as described in section card
        if (szInputString.substr(0,7) == "section" && m_nJouFlag == 0){
            std::cout << "Sectioning geometry .." << std::endl;
            char cDir;
            double dOffset;
            std::string szReverse = "";
            std::istringstream szFormatString (szInputString);
            szFormatString >> card >> cDir >> dOffset >> szReverse;
            Section_Assm(cDir, dOffset, szReverse);
            //ERRORR("Error in Section_Assm", err);
            std::cout <<"--------------------------------------------------"<<std::endl;

          }
        if (szInputString.substr(0,4) == "move" && m_nJouFlag == 0){
            std::cout << "Moving geometry .." << std::endl;
            double dX, dY, dZ;
            std::istringstream szFormatString (szInputString);
            szFormatString >> card >> dX >> dY >> dZ;
            if(szFormatString.fail())
              IOErrorHandler(INVALIDINPUT);
            Move_Assm(dX, dY, dZ);
            //ERRORR("Error in Move_Assm", err);
            std::cout <<"--------------------------------------------------"<<std::endl;

          }
        // center the assembly
        if (szInputString.substr(0,6) == "center" && m_nJouFlag == 0){

            char rDir = 'q';
            std::istringstream szFormatString (szInputString);
            szFormatString >> card >> rDir;
            if (rDir != 'q')
              std::cout << "Positioning assembly to "<< rDir << " center" << std::endl;
            else
              std::cout << "Positioning assembly to xy center" << std::endl;
            Center_Assm(rDir);
            //ERRORR("Error in Center_Assm", err);
            std::cout <<"--------------------------------------------------"<<std::endl;
          }
        // rotate the assembly if rotate card is specified
        if (szInputString.substr(0,6) == "rotate" && m_nJouFlag == 0){
            char cDir;
            double dAngle;
            std::cout << "Rotating geometry .." << std::endl;
            std::istringstream szFormatString (szInputString);
            szFormatString >> card >> cDir >> dAngle;
            if(szFormatString.fail())
              IOErrorHandler(INVALIDINPUT);
            Rotate_Assm(cDir, dAngle);
            //ERRORR("Error in Rotate_Assm", err);
            std::cout <<"--------------------------------------------------"<<std::endl;

          }
        // 'yes' or 'no' for creating sidesets
        if (szInputString.substr(0,13) == "createsideset"){
            std::istringstream szFormatString (szInputString);
            szFormatString >> card >> m_szSideset;
            std::cout <<"--------------------------------------------------"<<std::endl;
          }
        // Create specified number of files with varying material ids
        if (szInputString.substr(0,11) == "createfiles"){
            std::istringstream szFormatString (szInputString);
            szFormatString >> card >> m_nAssyGenInputFiles;
            std::cout <<"--------------------------------------------------"<<std::endl;
          }
        // Create specified number of files with varying material ids
        if (szInputString.substr(0,14) == "creatematfiles"){
            m_bCreateMatFiles = true;
            std::istringstream szFormatString (szInputString);
            szFormatString >> card >> m_nAssyGenInputFiles;
            std::cout <<"--------------------------------------------------"<<std::endl;
          }
        // Create specified number of files with varying material ids
        if (szInputString.substr(0,11) == "save_exodus"){
            save_exodus = true;
            std::cout <<"--------------------------------------------------"<<std::endl;
          }
        // specify a merge tolerance value for cubit journal file
        if (szInputString.substr(0,14) == "mergetolerance"){
            std::istringstream szFormatString (szInputString);
            szFormatString >> card >> m_dMergeTol;
            std::cout <<"--------------------------------------------------"<<std::endl;
          }
        // Handle mesh size inputs
        if (szInputString.substr(0,14) == "radialmeshsize"){
            std::istringstream szFormatString (szInputString);
            szFormatString >> card >> m_dRadialSize;
            if(m_dRadialSize < 0 || szFormatString.fail())
              IOErrorHandler(ENEGATIVE);
            std::cout <<"--------------------------------------------------"<<std::endl;

          }
        // Handle mesh size inputs
        if (szInputString.substr(0,11) == "tetmeshsize"){
            std::istringstream szFormatString (szInputString);
            szFormatString >> card >> m_dTetMeshSize;
            if(m_dTetMeshSize < 0 || szFormatString.fail())
              IOErrorHandler(ENEGATIVE);
            std::cout <<"--------------------------------------------------"<<std::endl;

          }
        // Handle mesh size inputs
        if (szInputString.substr(0,13) == "axialmeshsize"){
            std::istringstream szFormatString (szInputString);
            szFormatString >> card;
            if(m_nDuct > 0){
                m_dAxialSize.SetSize(m_nDuct);
                int num_ams_specified = std::distance(std::istream_iterator<std::string>(szFormatString),
                                                      std::istream_iterator<std::string>());
                std::istringstream szFormatStringAgain (szInputString);
                szFormatStringAgain >> card;
                for (int p = 1; p <= m_nDuct; p++){
                    if(p <= num_ams_specified)
                      szFormatStringAgain >> m_dAxialSize(p);
                    else
                      m_dAxialSize(p) = m_dAxialSize(num_ams_specified);
                    if(m_dAxialSize(p) < 0)
                      IOErrorHandler(ENEGATIVE);
                  }
              }
            else{
                m_dAxialSize.SetSize(1);
                szFormatString >> m_dAxialSize(1);
              }
            std::cout <<"--------------------------------------------------"<<std::endl;

          }
        // edge interval
        if (szInputString.substr(0, 12) == "edgeinterval") {
            std::istringstream szFormatString(szInputString);
            szFormatString >> card >> m_edgeInterval;
          }
        // Handle mesh size inputs
        if (szInputString.substr(0,18) == "neumannset_startid"){
            std::istringstream szFormatString (szInputString);
            szFormatString >> card >> m_nNeumannSetId;
            if(m_nNeumannSetId < 0 || szFormatString.fail())
              IOErrorHandler(ENEGATIVE);
            std::cout <<"--------------------------------------------------"<<std::endl;

          }
        // Handle mesh size inputs
        if (szInputString.substr(0,19) == "materialset_startid"){
            std::istringstream szFormatString (szInputString);
            szFormatString >> card >> m_nMaterialSetId;
            if(m_nMaterialSetId < 0 || szFormatString.fail())
              IOErrorHandler(ENEGATIVE);
            std::cout <<"--------------------------------------------------"<<std::endl;

          }
        if ((szInputString.substr(0,23) == "list_neumannset_startid") ){
            std::istringstream szFormatString (szInputString);
            int num_nset_ids = 0;
            szFormatString >> card >> num_nset_ids;
            m_nListNeuSet.SetSize(num_nset_ids);
            for (int p = 1; p <= num_nset_ids; p++){
                szFormatString >> m_nListNeuSet(p);
                if(m_nListNeuSet(p) < 0 || szFormatString.fail())
                  IOErrorHandler(ENEGATIVE);
              }
          }
        if ((szInputString.substr(0,14) == "numsuperblocks") ){
            std::istringstream szFormatString (szInputString);
            szFormatString >> card >> m_nSuperBlocks;
            sb.SetSize(m_nSuperBlocks);
          }
        if ((szInputString.substr(0,10) == "superblock") ){
            std::istringstream szFormatString (szInputString);
            szFormatString >> card >> sb(tmpSB).m_nSuperBlockId  >> sb(tmpSB).m_szSuperBlockAlias >> sb(tmpSB).m_nNumSBContents;
            sb(tmpSB).m_nSBContents.SetSize(sb(tmpSB).m_nNumSBContents);
            for (int p = 1; p <= sb(tmpSB).m_nNumSBContents; p++){
                szFormatString >> sb(tmpSB).m_nSBContents(p);
                if(sb(tmpSB).m_nSBContents(p) < 0 || szFormatString.fail())
                  IOErrorHandler(ENEGATIVE);
              }
            ++tmpSB;
          }
        if ((szInputString.substr(0,24) == "list_materialset_startid") ){
            std::istringstream szFormatString (szInputString);
            int num_mset_ids = 0;
            szFormatString >> card >> num_mset_ids;
            m_nListMatSet.SetSize(num_mset_ids);
            for (int p = 1; p <= num_mset_ids; p++){
                szFormatString >> m_nListMatSet(p);
                if(m_nListMatSet(p) < 0 || szFormatString.fail())
                  IOErrorHandler(ENEGATIVE);
              }
          }
        if (szInputString.substr(0,3) == "end" || m_nLineNumber == MAXLINES){


            if ( m_nJouFlag == 0){
                // impring merge before saving
                Imprint_Merge(m_bimprint, m_bmerge);

                clock_t s_save= clock();
                // save .sat file
                iGeom_save(igeomImpl->instance(), m_szGeomFile.c_str(), NULL, &err, m_szGeomFile.length() , 0);
                std::cout << "## Saving CPU time used := " << (double) (clock() - s_save)/CLOCKS_PER_SEC
                          << " seconds" << std::endl;

                m_PyCubGeomFile << "cubit.cmd('export acis \"" <<m_szGeomFile1 << "\" over')" << std::endl;

                std::cout << "Normal Termination.\n"<< "Geometry file: " << m_szGeomFile << " saved." << std::endl;
//                // Now run the journal file from the python script
//                m_PyCubGeomFile << "infile = open(\""<< m_szJouFile << "\", \"r\")" << std::endl;
//                m_PyCubGeomFile << "for line in infile:\n  cubit.cmd(line)" << std::endl;

                // Reloading file to check load times
                bool if_loadagain = false;
                if (if_loadagain == true){
                    clock_t s_load= clock();
                    iGeom_load(igeomImpl->instance(), m_szGeomFile.c_str(), NULL, &err,
                               strlen(m_szGeomFile.c_str()), 0);
                    std::cout << "## Load again CPU time used := " << (double) (clock() - s_load)/CLOCKS_PER_SEC
                              << " seconds" << std::endl;
                  }
              }
            break;
          }
      }

  }

  void AssyGen::CreateAssyGenInputFiles()
  //---------------------------------------------------------------------------
  //Function: Create Cubit Journal File for generating mesh
  //Input:    none
  //Output:   none
  //---------------------------------------------------------------------------
  {
    // create file names
    std::ostringstream os;
    std::string file, temp, temp1;
    int counter_ms = 0;
    int counter_ns = 0;
    for(int i=1; i<=m_nAssyGenInputFiles; i++){
        // form the name of the AssyGen input file
        if(m_bCreateMatFiles)
          os << m_nListMatSet(i) << ".inp";
        else
          os << m_szFile << i << ".inp";
        std::ofstream ofs;
        // open input file

        std::cout << os.str() << std::endl;
        bool bDone = false;
        do{
            file = os.str();
            ofs.open (file.c_str(), std::ios::out);
            if(!ofs){
                ofs.clear();
                std::cout << "Unable to open AssyGen Input File(s) for writing" << std::endl;
                exit(1);
              }
            else {
                bDone = true;
                std::cout << "File Opened" << std::endl;
              }
            os.str("");

            // write the input deck
            ofs << "! ## This is an automatically created AssyGen Input File: "<< file << std::endl;
            ofs << "MeshType " << m_szMeshType << std::endl;
            ofs << "GeomEngine " << m_szEngine << std::endl;
            ofs << "GeometryType " << m_szGeomType << std::endl;
            ofs << "Materials " << m_nAssemblyMat;

            // list all materials and alias with subscripts
            for (int j =1;j<=m_nAssemblyMat; j++){
                if(!m_bCreateMatFiles){
                    os << m_szAssmMat(j) << "_" << (i-1)*m_nAssemblyMat + j;
                    temp = os.str();

                    ofs << "  " << temp << " " << m_szAssmMatAlias(j);
                    os.str("");
                  }
                else{
                    ofs << "  " << m_szAssmMat(j) << " " << m_szAssmMatAlias(j);
                  }
              }
            ofs << "\n";

            //Rewind the input file
            int nDumpAllLinesAfter = 0, nMid = 0;
            m_FileInput.clear (std::ios_base::goodbit);
            m_FileInput.seekg (0L, std::ios::beg);
            m_nLineNumber = 0;
            CParser Parse;
            std::string card;

            // start reading the input file break when encounter end
            for(;;){
                if (!Parse.ReadNextLine (m_FileInput, m_nLineNumber, szInputString,
                                         MAXCHARS, szComment))
                  IOErrorHandler (INVALIDINPUT);
                // dump all the lines after duct command, limiting the format for writing AssyGen Files
                if( (szInputString.substr(0,10) == "dimensions") ||
                    (szInputString.substr(0,4) == "duct") ){
                    ofs << szInputString << std::endl;
                    nDumpAllLinesAfter = m_nLineNumber;
                  }
                else if ((szInputString.substr(0,19) == "materialset_startid") ){
                    nMid = (i-1)*m_nAssemblyMat + 1;
                    ofs << "MaterialSet_StartId " << nMid << std::endl;
                  }
                else if ((szInputString.substr(0,18) == "neumannset_startid") ){
                    nMid = (i-1)*m_nAssemblyMat + 1;
                    ofs << "NeumannSet_StartId " << nMid << std::endl;
                  }
                else if ((szInputString.substr(0,24) == "list_materialset_startid") ){
                    ++counter_ms;
                    ofs << "MaterialSet_StartId " << m_nListMatSet(counter_ms) << std::endl;
                  }
                else if ((szInputString.substr(0,23) == "list_neumannset_startid") ){
                    ++counter_ns;
                    ofs << "NeumannSet_StartId " << m_nListNeuSet(counter_ns) << std::endl;
                  }
                else if((szInputString.substr(0,11) == "createfiles")){
                    //skip this line
                  }
                else if((szInputString.substr(0,14) == "creatematfiles")){
                    //skip this line
                  }
                else if (nDumpAllLinesAfter > 0 && m_nLineNumber > nDumpAllLinesAfter){
                    ofs << szInputString << std::endl;
                  }

                if (szInputString.substr(0,3) == "end" || m_nLineNumber == MAXLINES){
                    break;
                  }
              }
            ofs.close();
          } while(!bDone);
      }



  }

  void AssyGen:: ComputePinCentroid(int nTempPin, CMatrix<std::string> MAssembly,
                                    int m, int n, double &dX, double &dY, double &dZ)
  // ---------------------------------------------------------------------------
  // Function: computes the centroid in the whole assembly of rectangular or hexagonal pincell
  // Input:    number and location of the pincell
  // Output:   coordinates of pin in assembly
  // ---------------------------------------------------------------------------
  {
    int nTempPin1 = -1, nTempPin2 = 0, nInputLines;
    std::string szVolId, szVolAlias;
    if(m_szGeomType == "hexagonal"){
        double dP, dZ;
        m_Pincell(nTempPin).GetPitch(dP, dZ);

        if (m < m_nPin){
            dX = (m_nPin - n + 1)*dP/2.0 + n*dP/2.0 + (n-1)*dP - (m-1)*dP/2.0;
            dY = (m-1)*(0.5*dP/sin(pi/3.0) + 0.5*dP*sin(pi/6.0)/sin(pi/3.0));
          }
        else{
            dX = (m_nPin - n + 1)*dP/2.0 + n*dP/2.0 + (n-1)*dP - (2*m_nPin - m -1)*dP/2.0;
            dY = (m-1)*(0.5*dP/sin(pi/3.0) + 0.5*dP*sin(pi/6.0)/sin(pi/3.0));
          }
      }
    if(m_szGeomType == "rectangular"){
        double dPX, dPY, dPZ, dPX1, dPY1, dPZ1, dPX2, dPY2, dPZ2;
        if((m_Assembly(m,n)=="x")||(m_Assembly(m,n)=="xx"))
          m_Pincell(1).GetPitch(dPX, dPY, dPZ);
        else
          m_Pincell(nTempPin).GetPitch(dPX, dPY, dPZ);

        if (n==1){
            dX = 0;
            if(m==1)
              dY = 0;
          }
        else{
            dX+= dPX/2.0;
            // find the previous pincell type
            // check if it's dummy
            if((m_Assembly(m,n-1)=="x")||(m_Assembly(m,n-1)=="xx")){
                dX+=dPX/2.0;
              }
            else{
                for(int b=1; b<=m_nPincells; b++){
                    m_Pincell(b).GetLineOne(szVolId, szVolAlias, nInputLines);
                    if(m_Assembly(m,n-1) == szVolAlias)
                      nTempPin1 = b;
                  }
                if(nTempPin1 == -1){
                    std::cout << "Unknown pincell, pincell " << m_Assembly(m,n-1) << " not declared" << std::endl;
                    exit(1);
                  }
                m_Pincell(nTempPin1).GetPitch(dPX1, dPY1, dPZ1);
                // now add half of X pitch to the previous cells pitch
                dX+= dPX1/2.0;
              }
          }
        if (m > 1 && n==1){
            dY+= dPY/2.0;
            // check if it's dummy
            if((m_Assembly(m-1,n)=="x")||(m_Assembly(m-1,n)=="xx")){
                dY+= dPY/2.0;
              }
            else{
                for(int c=1; c<=m_nPincells; c++){
                    m_Pincell(c).GetLineOne(szVolId, szVolAlias, nInputLines);
                    if(m_Assembly(m-1,n) == szVolAlias)
                      nTempPin2 = c;
                  }
                m_Pincell(nTempPin2).GetPitch(dPX2, dPY2, dPZ2);
                dY+= dPY2/2.0;
              }
          }
        dZ = 0.0; // moving in XY plane only
      }//if rectangular ends

  }

  void AssyGen::IOErrorHandler (ErrorStates ECode) const
  // ---------------------------------------------------------------------------
  // Function: displays error messages related to input data
  // Input:    error code
  // Output:   none
  // ---------------------------------------------------------------------------
  {
    std::cerr << '\n';

    if (ECode == PINCELLS) // invalid number of pincells
      std::cerr << "Number of pincells must be >= 0.";
    else if (ECode == INVALIDINPUT) // invalid input
      std::cerr << "Invalid input.";
    else if (ECode == EMAT) // invalid input
      std::cerr << "Invalid Material Data.";
    else if (ECode == EGEOMTYPE) // invalid input
      std::cerr << "Invalid GeomType Data.";
    else if (ECode == EGEOMENGINE) // invalid input
      std::cerr << "Invalid Geometry Engine.";
    else if (ECode == EALIAS) // invalid input
      std::cerr << "Error Reading Aliases.";
    else if (ECode == ENEGATIVE) // invalid input
      std::cerr << "Unexpected negative value.";
    else if (ECode == EPIN) // invalid input
      std::cerr << "Invalid pinCell specs.";
    else if (ECode == EUNEQUAL) // invalid input
      std::cerr << "Number of cyliders and ducts in Z don't match, check .inp file.";
    else
      std::cerr << "Unknown error ...?";

    std::cerr << '\n' << "Error in input file line : " << m_nLineNumber;
    std::cerr << std::endl;
    exit (1);
  }

  void AssyGen:: Name_Faces(const std::string sMatName, const iBase_EntityHandle body,  iBase_TagHandle this_tag )
  // ---------------------------------------------------------------------------
  // Function: names all the faces in the body
  // Input:    none
  // Output:   none
  // ---------------------------------------------------------------------------
  {
    double dTol = 1e-4, ttol = 1e-2;
    if(m_dMAssmPitch.GetRows()!=0 && m_dMAssmPitch.GetColumns()!=0){
        if(m_szGeomType == "hexagonal")
          ttol = m_dMAssmPitch(1, 1);
        else if (m_szGeomType == "rectangular")
          ttol = m_dMAssmPitchX(1,1);
      }
    // set tolerance for surface identification
    if (ttol != 0){
        dTol=ttol*1.0e-2;
      }
    double dZTemp = 0.0;
    int flag = 0, locTemp = 0;
    iBase_EntityHandle max_surf = NULL, min_surf = NULL, side_surf =NULL;
    SimpleArray<iBase_EntityHandle> surfs;
    int nSide = 0;
    std::ostringstream os;
    std::string sMatName0=sMatName+"_top";
    std::string sMatName1=sMatName+"_bot";
    std::string sSideName;
    iGeom_getEntAdj( igeomImpl->instance(), body, iBase_FACE, ARRAY_INOUT(surfs), &err );
    //CHECK( "Problems getting max surf for rotation." );

    SimpleArray<double> max_corn, min_corn;
    iGeom_getArrBoundBox( igeomImpl->instance(), ARRAY_IN(surfs), iBase_INTERLEAVED,
                          ARRAY_INOUT( min_corn ),
                          ARRAY_INOUT( max_corn ),
                          &err );
    //CHECK( "Problems getting max surf for rotation." );
    for (int i = 0; i < surfs.size(); ++i){
        // first find the max z-coordinate
        if( (fabs(min_corn[3*i+2]-max_corn[3*i+2])) < dTol ) {
            if(flag == 0){
                dZTemp = min_corn[3*i+2];
                locTemp = i;
                flag = 1;
              }
            else if(dZTemp > min_corn[3*i+2]){
                // we have a bot surface
                min_surf = surfs[i];
                // the top surface is dZTemp
                max_surf = surfs[locTemp];
              }
            else{
                //we have a top surface
                min_surf = surfs[locTemp];
                // the top surface is dZTemp
                max_surf = surfs[i];
              }
          }
        // see if max or min set name
        if(max_surf !=0){
            iGeom_setData(igeomImpl->instance(), max_surf, this_tag,
                          sMatName0.c_str(), sMatName0.size(), &err);
            ////CHECK("setData failed");

            std::cout << sMatName0 << ",  ";
            max_surf = NULL;

          }
        if(min_surf !=0){
            iGeom_setData(igeomImpl->instance(), min_surf, this_tag,
                          sMatName1.c_str(), sMatName1.size(), &err);
            ////CHECK("setData failed");
            std::cout << sMatName1 << ",  ";
            min_surf = NULL;

          }
      }
    for (int i = 0; i < surfs.size(); ++i){
        if( (fabs(min_corn[3*i+2]-max_corn[3*i+2])) < dTol ) {
            continue; // its a max of min surface
          }
        else{ // its a side surface
            side_surf = surfs[i];
          }
        //set name for the sidesurf now
        if(side_surf !=0){
            ++nSide;
            sSideName = sMatName + "_side";
            if(m_szGeomType == "hexagonal") {
                if(nSide <= 6)
                  os << sSideName << nSide;
                else
                  os << sSideName << "_" << nSide;
              }

            if(m_szGeomType == "rectangular"){
                if(nSide <= 4)
                  os << sSideName << nSide;
                else
                  os << sSideName << "_" << nSide;
              }
            sSideName = os.str();
            iGeom_setData(igeomImpl->instance(), side_surf, this_tag,
                          sSideName.c_str(), sMatName1.size(), &err);
            ////CHECK("setData failed");
            std::cout << sSideName << ",  " ;
            sSideName = "";
            os.str("");
            side_surf = NULL;
          }
        else {
            std::cerr << "Couldn't find surface for naming" << std::endl;
          }
      }
    std::cout <<"\n";

  }


  void AssyGen::Center_Assm (char &rDir)
  // ---------------------------------------------------------------------------
  // Function: centers all the entities along x and y axis
  // Input:    none
  // Output:   none
  // ---------------------------------------------------------------------------
  {
    double xmin, xmax, ymin, ymax, zmin, zmax, xcenter = 0.0, ycenter = 0.0, zcenter = 0.0;
    // position the assembly such that origin is at the center before sa
    iGeom_getBoundBox(igeomImpl->instance(),&xmin,&ymin,&zmin,
                      &xmax,&ymax,&zmax, &err);
    ////CHECK("Failed getting bounding box");

    // moving all geom entities to center


    if( rDir =='x'){
        xcenter = (xmin+xmax)/2.0;
      }
    else if( rDir =='y'){
        ycenter = (ymin+ymax)/2.0;
      }
    else if ( rDir =='z'){
        zcenter = (zmin+zmax)/2.0;
      }
    else{
        // assume that it is centered along x and y and not z direction
        xcenter = (xmin+xmax)/2.0;
        ycenter = (ymin+ymax)/2.0;
      }

    SimpleArray<iBase_EntityHandle> all;
    iGeom_getEntities( igeomImpl->instance(), root_set, iBase_REGION,ARRAY_INOUT(all),&err );
    ////CHECK("Failed to get all entities");

    for(int i=0; i<all.size(); i++){
        iGeom_moveEnt(igeomImpl->instance(),all[i],-xcenter,-ycenter,-zcenter,&err);
      }
      // OCC bounding box computation is buggy, better to compute bounding box in python and supply to the script.
      m_PyCubGeomFile << "vol = cubit.get_entities(\"volume\")" << std::endl;
      m_PyCubGeomFile << "vl = cubit.get_total_bounding_box(\"volume\", vol)\nzcenter = 0.0" << std::endl;

    if( rDir =='x'){
        m_PyCubGeomFile << "xcenter = (vl[0]+vl[1])/2.0" << std::endl;
      }
    else if( rDir =='y'){
        m_PyCubGeomFile << "ycenter = (vl[3]+vl[4])/2.0" << std::endl;
      }
    else if ( rDir =='z'){
        m_PyCubGeomFile << "zcenter = (vl[6]+vl[7])/2.0" << std::endl;
      }
    else{
        // assume that it is centered along x and y and not z direction
        m_PyCubGeomFile << "xcenter = (vl[0]+vl[1])/2.0" << std::endl;
        m_PyCubGeomFile << "ycenter = (vl[3]+vl[4])/2.0" << std::endl;
      }

      m_PyCubGeomFile << "cubit.cmd('move vol all x -{0} y -{1} z -{2}'.format(xcenter, ycenter, zcenter) )" <<  std::endl;
  }

  void AssyGen::Section_Assm (char &cDir, double &dOffset, const std::string szReverse)
  // ---------------------------------------------------------------------------
  // Function: sections the assembly about the cutting plane
  // Input:    none
  // Output:   none
  // ---------------------------------------------------------------------------
  {
    double xmin, xmax, ymin, ymax, zmin, zmax, yzplane = 0.0, xzplane = 0.0;
    iBase_EntityHandle sec = NULL;
    int nReverse = 0;
    // check if reverse side is needed
    if(szReverse == "reverse"){
        nReverse = 1;
      }
    if( cDir =='x'){
        yzplane = 1.0;
        xzplane = 0.0;
      }
    if( cDir =='y'){
        yzplane = 0.0;
        xzplane = 1.0;
      }
    SimpleArray<iBase_EntityHandle> all;
    iGeom_getEntities( igeomImpl->instance(), root_set, iBase_REGION,ARRAY_INOUT(all),&err );
    ////CHECK("Failed to get all entities");
    // loop and section/delete entities
    for(int i=0; i < all.size(); i++){
        //get the bounding box to decide
        iGeom_getEntBoundBox(igeomImpl->instance(),all[i],&xmin,&ymin,&zmin,
                             &xmax,&ymax,&zmax, &err);
        ////CHECK("Failed get bound box");
        if(xmin > dOffset && yzplane ==1 && nReverse ==1){
            iGeom_deleteEnt(igeomImpl->instance(),all[i],&err);
            ////CHECK("Failed delete entities");
            continue;
          }
        if(ymin > dOffset && xzplane == 1 && nReverse ==1){
            iGeom_deleteEnt(igeomImpl->instance(),all[i],&err);
            ////CHECK("Failed delete entities");
            continue;
          }
        if(xmax < dOffset && yzplane ==1 && nReverse ==0){
            iGeom_deleteEnt(igeomImpl->instance(),all[i],&err);
            ////CHECK("Failed delete entities");
            continue;
          }
        if(ymax < dOffset && xzplane == 1 && nReverse ==0){
            iGeom_deleteEnt(igeomImpl->instance(),all[i],&err);
            ////CHECK("Failed delete entities");
            continue;
          }
        else{
            if(xzplane ==1 && ymax >dOffset && ymin < dOffset){
                iGeom_sectionEnt(igeomImpl->instance(), all[i],yzplane,xzplane,0, dOffset, nReverse,&sec,&err);
                ////CHECK("Failed to section ent");
              }
            if(yzplane ==1 && xmax >dOffset && xmin < dOffset){
                iGeom_sectionEnt(igeomImpl->instance(), all[i],yzplane,xzplane,0, dOffset,nReverse,&sec,&err);
                ////CHECK("Failed to section ent");
              }
          }
      }

      m_PyCubGeomFile << "reserve = " << szReverse << "\ndOffset = " << dOffset << std::endl;
      if(xzplane ==1 && ymax >dOffset && ymin < dOffset){
        m_PyCubGeomFile << "cubit.cmd('section vol all with xplane offset {0} {1}'.format(dOffset, reverse))" << std::endl;
      }
      if(yzplane ==1 && xmax >dOffset && xmin < dOffset){
        m_PyCubGeomFile << "cubit.cmd('section vol all with yplane offset {0} {1}'.format(dOffset, reverse))" << std::endl;
      }

  }

  void AssyGen::Rotate_Assm (char &cDir, double &dAngle)
  // ---------------------------------------------------------------------------
  // Function: rotates the whole assembly
  // Input:    none
  // Output:   none
  // ---------------------------------------------------------------------------
  {
    double dX = 0.0, dY=0.0, dZ=0.0;
    if( cDir =='x'){
        dX = 1.0;
      }
    if( cDir =='y'){
        dY = 1.0;
      }
    if( cDir =='z'){
        dZ = 1.0;
      }
    SimpleArray<iBase_EntityHandle> all;
    iGeom_getEntities( igeomImpl->instance(), root_set, iBase_REGION,ARRAY_INOUT(all),&err );

    m_PyCubGeomFile << "cDir = '" << cDir << "'" << "\ncubit.cmd('rotate vol all angle 30 about {0}'.format(cDir))" << std::endl;

    // loop and rotate all entities
    for(int i=0; i<all.size(); i++){
        //get the bounding box to decide
        iGeom_rotateEnt(igeomImpl->instance(),all[i],dAngle,
                        dX, dY, dZ, &err);
      }

  }

  void AssyGen::Move_Assm (double &dX,double &dY, double &dZ)
  // ---------------------------------------------------------------------------
  // Function: move's the model by dX, dY and dZ
  // Input:    none
  // Output:   none
  // ---------------------------------------------------------------------------
  {
    SimpleArray<iBase_EntityHandle> all;
    iGeom_getEntities( igeomImpl->instance(), root_set, iBase_REGION,ARRAY_INOUT(all),&err );

    m_PyCubGeomFile << "cubit.cmd('move vol all x " << dX << " y " << dY << " dZ " << dZ << "')" << std::endl;        

    // loop and rotate all entities
    for(int i=0; i<all.size(); i++){
        //get the bounding box to decide
        iGeom_moveEnt(igeomImpl->instance(),all[i],
                      dX, dY, dZ, &err);
      }

  }

  void AssyGen::Create_HexAssm(std::string &szInputString)
  // ---------------------------------------------------------------------------
  // Function: read and create the assembly for hexagonal lattice
  // Input:    error code
  // Output:   none
  // ---------------------------------------------------------------------------
  {
    CParser Parse;
    std::string card, szVolId, szVolAlias;
    int nInputLines, nTempPin = 1, t, nIFlag = 0, total_pincells = 0;
    double dX = 0.0, dY =0.0, dZ=0.0;
    double  dP, dH, dSide, dHeight;
    iBase_EntityHandle assm = NULL;
    std::cout << "\ngetting Assembly data and creating ..\n"<< std::endl;
    std::istringstream szFormatString (szInputString);

    szFormatString >> card >> m_nPin;
    if(m_nPin <=0 )
      IOErrorHandler (INVALIDINPUT);

    // width of the hexagon is n*n+1/2
    int nWidth =2*m_nPin -1;

    // creating a square array of size width
    m_Assembly.SetSize(nWidth, nWidth);
    if (m_nJouFlag == 1)
      return;

    for(int m=1; m<=nWidth; m++){
        if (!Parse.ReadNextLine (m_FileInput, m_nLineNumber, szInputString,
                                 MAXCHARS, szComment))
          IOErrorHandler (INVALIDINPUT);
        if(m>m_nPin)
          t = 2*m_nPin - m;
        else
          t = m;
        std::istringstream szFormatString1 (szInputString);

        for(int n=1; n<=(m_nPin + t - 1); n++){
            ++total_pincells;
            nTempPin = -1;
            szFormatString1 >> m_Assembly(m,n);
            if(szFormatString1.fail())
              IOErrorHandler (INVALIDINPUT);
            // if dummy pincell skip and continue
            if((m_Assembly(m,n)=="x")||(m_Assembly(m,n)=="xx")){
                continue;
              }
            // find that pincell
            ++m_nTotalPincells;
            for(int b=1; b<=m_nPincells; b++){
                m_Pincell(b).GetLineOne(szVolId, szVolAlias, nInputLines);
                if(m_Assembly(m,n) == szVolAlias)
                  nTempPin = b;
              }
            if(nTempPin == -1){
                std::cout << "Unknown pincell, pincell " << m_Assembly(m,n) << " not declared" << std::endl;
                exit(1);
              }

            // now compute the location and create it
            ComputePinCentroid(nTempPin, m_Assembly, m, n, dX, dY, dZ);
            //ERRORR("Error in function ComputePinCentroid", err);

            // now create the pincell in the location found
            std::cout << "\n--------------------------------------------------"<<std::endl;
            std::cout << " m = " << m <<" n = " << n << std::endl;
            std::cout << "creating pin: " << nTempPin;
            std::cout << " at X Y Z " << dX << " " << dY << " " << dZ << std::endl;

            if(strcmp(m_szInfo.c_str(),"on") == 0)
              m_AssmInfo << nTempPin  << " \t" << m << " \t" << n << " \t" << dX << " \t" << dY << " \t" << dZ << std::endl;

            m_Pincell(nTempPin).GetIntersectFlag(nIFlag);
            if(nIFlag){
                CreatePinCell_Intersect(nTempPin, dX, -dY, dZ);
                //ERRORR("Error in function CreatePinCell_Intersect", err);
              }
            else{
                CreatePinCell(nTempPin, dX, -dY, dZ);
                //ERRORR("Error in function CreatePinCell", err);
              }
          }
      }

    // get all the entities (in pins)defined so far, in an entity set - for subtraction later
    iGeom_getEntities( igeomImpl->instance(), root_set, iBase_REGION, ARRAY_INOUT(in_pins),&err );
    //CHECK( "ERROR : getRootSet failed!" );
    std::cout << "Expected pin definitions: " << total_pincells << "\n\nCreating surrounding outer hexes .." << std::endl;

    for (int nTemp = 1; nTemp <= m_nDuct; nTemp ++){
        if(m_nDimensions >0){

            // create outermost hexes
            for(int n=1;n<=m_nDimensions; n++){
                dSide = m_dMAssmPitch(nTemp, n)/(sqrt(3));
                dHeight = m_dMZAssm(nTemp, 2) - m_dMZAssm(nTemp, 1);

                // creating coverings
                iGeom_createPrism(igeomImpl->instance(), dHeight, 6,
                                  dSide, dSide,
                                  &assm, &err);
                m_PyCubGeomFile << "##\nassm = cubit.prism(" << dHeight << ", 6, " << dSide << ", " << dSide << ")" << std::endl;

                // rotate the prism to match the pins
                iGeom_rotateEnt (igeomImpl->instance(), assm, 30, 0, 0, 1, &err);
                m_PyCubGeomFile << "cubit.cmd('rotate body {0} angle 30 about z'.format(assm.id()) )" << std::endl;

                if(0 != m_Pincell.GetSize()){
                    m_Pincell(1).GetPitch(dP, dH);
                    dX = m_nPin*dP;
                    dY = -(m_nPin-1)*dP*sqrt(3.0)/2.0;
                  }
                else{
                    dX = 0.0;
                    dY = 0.0;
                  }
                dZ = (m_dMZAssm(nTemp, 2) + m_dMZAssm(nTemp, 1))/2.0;

                // position the prism
                iGeom_moveEnt(igeomImpl->instance(), assm, dX,dY,dZ, &err);
                m_PyCubGeomFile << "vector = [" << dX << ", " << dY << ", " << dZ << "]" << std::endl;
                m_PyCubGeomFile << "cubit.move(assm, vector)" << std::endl;

                // populate the coverings array
                int loc = (nTemp-1)*m_nDimensions + n -1;
                assms[(nTemp-1)*m_nDimensions + n -1]=assm;
                m_PyCubGeomFile << "assms.insert(" << loc << ", assm)" << std::endl;
              }
          }
      }

  }

  void AssyGen::Create_CartAssm(std::string &szInputString)
  // ---------------------------------------------------------------------------
  // Function: read and create the assembly for rectangular lattice
  // Input:    error code
  // Output:   none
  // ---------------------------------------------------------------------------
  {
    CParser Parse;
    std::string card, szVolId, szVolAlias;
    int nInputLines, nTempPin = 1, nIFlag = 0.0;
    double dX = 0.0, dY =0.0, dZ=0.0, dMoveX = 0.0, dMoveY = 0.0, dHeight = 0, dPX=0.0, dPY=0.0, dPZ=0.0;
    iBase_EntityHandle assm = NULL;

    std::istringstream szFormatString (szInputString);
    szFormatString >> card >> m_nPinX >> m_nPinY;
    if(m_nPinX <=0 || m_nPinY <=0)
      IOErrorHandler (INVALIDINPUT);
    m_Assembly.SetSize(m_nPinY,m_nPinX);

    if (m_nJouFlag == 1)
      return;

    //read the next line to get assembly info &store assembly info
    if(0 != m_Pincell.GetSize()){
        for(int m=1; m<=m_nPinY; m++){
            if (!Parse.ReadNextLine (m_FileInput, m_nLineNumber, szInputString,
                                     MAXCHARS, szComment))
              IOErrorHandler (INVALIDINPUT);
            std::istringstream szFormatString1 (szInputString);

            //store the line read in Assembly array and create / position the pin in the core
            for(int n=1; n<=m_nPinX; n++){
                szFormatString1 >> m_Assembly(m,n);
                if(szFormatString1.fail())
                  IOErrorHandler (INVALIDINPUT);


                // loop thro' all pins to get the type of pin
                for(int b=1; b<=m_nPincells; b++){
                    m_Pincell(b).GetLineOne(szVolId, szVolAlias, nInputLines);
                    if(m_Assembly(m,n) == szVolAlias)
                      nTempPin = b;
                  }

                //now compute the location where the pin needs to be placed
                ComputePinCentroid(nTempPin, m_Assembly, m, n, dX, dY, dZ);
                //ERRORR("Error in function ComputePinCentroid", err);

                // if dummy pincell skip and continue
                if((m_Assembly(m,n)=="x")||(m_Assembly(m,n)=="xx")){
                    m_Pincell(1).GetPitch(dPX, dPY, dPZ);
                    // dMoveX and dMoveY are stored for positioning the outer squares later
                    if(m == m_nPinY && n ==m_nPinX){
                        dMoveX = dX/2.0;
                        dMoveY = -dY/2.0;
                      }
                    continue;
                  }
                ++m_nTotalPincells;
                // now create the pincell in the location found
                std::cout << "\n--------------------------------------------------"<<std::endl;
                std::cout << " m = " << m <<" n = " << n << std::endl;
                std::cout << "creating pin: " << nTempPin;
                std::cout << " at X Y Z " << dX << " " << -dY << " " << dZ << std::endl;

                if(strcmp(m_szInfo.c_str(),"on") == 0)
                  m_AssmInfo << nTempPin  << " \t" << m << " \t" << n << " \t" << dX << " \t" << dY << " \t" << dZ << std::endl;

                m_Pincell(nTempPin).GetIntersectFlag(nIFlag);
                if(nIFlag){
                    CreatePinCell_Intersect(nTempPin, dX, -dY, dZ);
                    //ERRORR("Error in function CreatePinCell_Intersect", err);
                  }
                else{
                    CreatePinCell(nTempPin, dX, -dY, dZ);
                    //ERRORR("Error in function CreatePinCell", err);
                  }
                // dMoveX and dMoveY are stored for positioning the outer squares later
                if(m == m_nPinY && n ==m_nPinX){
                    dMoveX = dX/2.0;
                    dMoveY = -dY/2.0;
                  }
              }
          }
      }
    std::cout << "\n--------------------------------------------------"<<std::endl;

    // get all the entities (in pins)defined so far, in an entity set - for subtraction later
    //  iGeom_getEntities( igeomImpl->instance(), root_set, iBase_REGION, ARRAY_INOUT(in_pins),&err );
    //  //CHECK( "ERROR : getRootSet failed!" );


    if(m_nDimensions > 0){

        // create outermost rectangular blocks
        std::cout << "\nCreating surrounding outer blocks .." << std::endl;
        int nCount = -1;
        for(int nTemp = 1; nTemp <= m_nDuct; nTemp ++){
            for(int n=1;n<=m_nDimensions; n++){
                ++nCount;
                dHeight = m_dMZAssm(nTemp, 2) - m_dMZAssm(nTemp, 1);
                iGeom_createBrick(igeomImpl->instance(), m_dMAssmPitchX(nTemp, n),  m_dMAssmPitchY(nTemp, n), dHeight,
                                  &assm, &err);
                m_PyCubGeomFile << "assm = cubit.brick( " << m_dMAssmPitchX(nTemp, n) << ", " << m_dMAssmPitchY(nTemp, n) << ", " << dHeight << ")" << std::endl;

                // position the outer block to match the pins
                dZ = (m_dMZAssm(nTemp, 2) + m_dMZAssm(nTemp, 1))/2.0;
                std::cout << "Move " <<   dMoveX << " " << dMoveY <<std::endl;
                iGeom_moveEnt(igeomImpl->instance(), assm, dMoveX,dMoveY,dZ, &err);
                m_PyCubGeomFile << "vector = [" << dMoveX << ", " << dMoveY << ", " << dZ << "]" << std::endl;
                m_PyCubGeomFile << "cubit.move(assm, vector)" << std::endl;

                // populate the outer covering array squares
                assms[nCount]=assm;
                m_PyCubGeomFile << "assms.insert(" << nCount << ", assm)" << std::endl;

              }
          }
      }

  }

  void AssyGen::CreateOuterCovering ()
  // ---------------------------------------------------------------------------
  // Function: this function sets the names of the coverings
  // Input:    error code
  // Output:   none
  // ---------------------------------------------------------------------------
  {
    double xmin, xmax, ymin, ymax, zmin, zmax;
    iBase_TagHandle this_tag;
    char* tag_name =(char *)"NAME";
    std::string sMatName = "";
    std::string sMatName0 = "";
    std::string sMatName1 = "";

    // get tag handle for 'NAME' tag, already created as iGeom instance is created
    iGeom_getTagHandle(igeomImpl->instance(), tag_name, &this_tag, &err, 4);
    ////CHECK("getTagHandle failed");
    iBase_EntityHandle tmp_vol= NULL, tmp_new= NULL;

    // name the innermost outer covering common for both rectangular and hexagonal assembliees
    if(m_nDimensions >0){
        //    int tag_no = 0;
        for (int nTemp1 = 1; nTemp1 <=m_nDuct; nTemp1++){
            for(int p=1;p<=m_szAssmMatAlias.GetSize();p++){
                if(strcmp ( m_szMMAlias(nTemp1, 1).c_str(), m_szAssmMatAlias(p).c_str()) == 0){
                    sMatName =  m_szAssmMat(p);
                    //    tag_no=p;
                  }
              }

            std::cout << "\ncreated innermost block: " << sMatName << std::endl;

            tmp_vol = assms[(nTemp1 - 1)*m_nDimensions];
            iGeom_setData(igeomImpl->instance(), tmp_vol, this_tag,
                          sMatName.c_str(), sMatName.size(), &err);
            ////CHECK("setData failed");

            Name_Faces(sMatName, tmp_vol, this_tag);
            m_PyCubGeomFile << "lid = assms[" << (nTemp1 - 1)*m_nDimensions << "].id()" << std::endl;
            m_PyCubGeomFile << "cubit.set_entity_name(\"body\", lid, \""  << sMatName <<  "\" )" << std::endl;
            m_PyCubGeomFile << "name_faces(\"" << sMatName << "\", assms[" << (nTemp1 - 1)*m_nDimensions << "])" << std::endl;
          }

        int count =0;//index for edge names
        for (int nTemp = 1; nTemp <= m_nDuct; nTemp++){
            //  Naming outermost block edges - sidesets in cubit journal file
            std::cout << "Naming outermost block edges" << std::endl;
            SimpleArray<iBase_EntityHandle> edges;

            iGeom_getEntAdj( igeomImpl->instance(), assms[nTemp*m_nDimensions-1] , iBase_EDGE,ARRAY_INOUT(edges),
                &err );
            //CHECK( "ERROR : getEntAdj failed!" );

            // get the top corner edges of the outer most covering
            //m_PyCubGeomFile << "lid=assms<<["<< nTemp*m_nDimensions-1 << "].id()" << std::endl;
            m_PyCubGeomFile << "cubit.cmd('group \"g1\" equals curve in vol {0} '.format(assms[" << nTemp*m_nDimensions-1 << "].id()))" << std::endl;
            m_PyCubGeomFile << "cubit.cmd('group \"g2\" equals curve with z_max<> z_min in g1')\ncubit.cmd('group  \"g3\" subtract g2 from g1')" << std::endl;
            m_PyCubGeomFile << "cubit.cmd('curve in g3 name \"side_edge\"')" << std::endl;
            std::ostringstream os;
            for (int i = 0; i < edges.size(); ++i){
                iGeom_getEntBoundBox(igeomImpl->instance(), edges[i],&xmin,&ymin,&zmin,
                                     &xmax,&ymax,&zmax, &err);
                ////CHECK("getEntBoundBox failed.");
                double dTol = 1e-5; // tolerance for comparing coordinates

                if(fabs(zmax - m_dMZAssm(nTemp, 2)) <  dTol){
                    if(fabs(zmax-zmin) < dTol){

                        //we have a corner edge - name it
                        sMatName="side_edge";
                        ++count;
                        os << sMatName << count;
                        sMatName=os.str();
                        tmp_vol=edges[i];
                        iGeom_setData(igeomImpl->instance(), tmp_vol, this_tag,
                                      sMatName.c_str(), sMatName.size(), &err);
                        ////CHECK("setData failed");
                        std::cout << "created: " << sMatName << std::endl;
                        os.str("");
                        sMatName="";
                      }
                  }
              }
          }
        // now subtract the outermost hexes and name them
        std::cout << "Subtract outermost hexes and naming them" << std::endl;
        int nCount = 0;
        for(int nTemp=1; nTemp<=m_nDuct; nTemp++){
            for(int n=m_nDimensions; n>1 ; n--){
                if(n>1){
                    ++nCount;
                    // copy cyl before subtract
                    iGeom_copyEnt(igeomImpl->instance(), assms[(nTemp-1)*m_nDimensions + n-2], &tmp_vol, &err);
                    m_PyCubGeomFile << "tmp_vol = cubit.copy_body(assms[" << (nTemp-1)*m_nDimensions + n-2 << "])" << std::endl;

                    // subtract outer most cyl from brick
                    m_PyCubGeomFile << "\nsub1.append(tmp_vol)" << std::endl;
                    m_PyCubGeomFile << "\nsub2.append(assms[" << (nTemp-1)*m_nDimensions + n-1 <<"])" << std::endl;

                    iGeom_subtractEnts(igeomImpl->instance(), assms[(nTemp-1)*m_nDimensions + n-1], tmp_vol, &tmp_new, &err);

                    m_PyCubGeomFile << "tmp_new = cubit.subtract(sub1, sub2)" << std::endl;
                    m_PyCubGeomFile << "assms[" << (nTemp-1)*m_nDimensions + n-1 << "] = tmp_new[0]\n\nsub1[:]=[]\nsub2[:]=[]"   << std::endl;

                    assms[(nTemp-1)*m_nDimensions + n-1]=tmp_new;

                    // name the vols by searching for the full name of the abbreviated Cell Mat
                    for(int p=1;p<=m_szAssmMatAlias.GetSize();p++){
                        if(strcmp ( m_szMMAlias(nTemp, n).c_str(), m_szAssmMatAlias(p).c_str()) == 0){
                            sMatName =  m_szAssmMat(p);
                          }
                      }
                    std::cout << "created: " << sMatName << std::endl;

                    iGeom_setData(igeomImpl->instance(), tmp_new, this_tag,
                                  sMatName.c_str(), sMatName.size(), &err);
                    m_PyCubGeomFile << "lid = tmp_new[0].id()" << std::endl;
                    m_PyCubGeomFile << "cubit.set_entity_name(\"body\", lid, \""  << sMatName <<  "\" )" << std::endl;
                    m_PyCubGeomFile << "name_faces(\"" << sMatName << "\", tmp_new[0]) " << std::endl;
                    Name_Faces(sMatName, tmp_new, this_tag);
                  }
              }
          }
        std::cout << "\n--------------------------------------------------"<<std::endl;
      }

  }

  void AssyGen::Subtract_Pins()
  // ---------------------------------------------------------------------------
  // Function: subtract the pins from the outer block
  // Input:    none
  // Output:   none
  // ---------------------------------------------------------------------------
  {
    if (m_nDimensions >0){
        std::cout <<"Total number of pins in the model = " << m_nTotalPincells << std::endl;

        for (int k=1; k<=m_nDuct; k++){
            if(cp_inpins[k-1].size() ==0)
              continue;
            // put all the in pins in a matrix of size duct for subtraction with ducts
            std::vector <iBase_EntityHandle> pin_copy( cp_inpins[k-1].size(), NULL);
            m_PyCubGeomFile << "pin_copy=[]\n\nsub1=[]\nsub2=[]\n" << std::endl;
            for (int i=0; i< (int) cp_inpins[k-1].size();i++){
                iGeom_copyEnt(igeomImpl->instance(), cp_inpins[k-1][i], &pin_copy[i], &err);
                m_PyCubGeomFile << "tmp_vol = cubit.copy_body(cp_inpins[" << k-1 << "][" << i << "])" << std::endl;
                m_PyCubGeomFile << "pin_copy.append(tmp_vol)" << std::endl;
                m_PyCubGeomFile << "\nsub2.append(cp_inpins[" << k-1 << "]["<< i << "])" << std::endl;
              }

            iBase_EntityHandle tmp_vol = NULL;
            tmp_vol = assms[(k-1)*m_nDimensions];
            m_PyCubGeomFile << "tmp_vol = assms[" << k << "]" << std::endl;

            // subtract the innermost hex from the pins
            std::cout << "Duct no.: " << k << " subtracting " <<  cp_inpins[k-1].size() << " pins from the duct .. " << std::endl;

            //#if HAVE_ACIS
            iBase_EntityHandle unite= NULL, tmp_new1;

            // if there are more than one pins
            if( cp_inpins[k-1].size() > 1){

               iGeom_uniteEnts(igeomImpl->instance(), &cp_inpins[k-1][0], cp_inpins[k-1].size(), &unite, &err);
                //m_PyCubGeomFile << "##\nunitepins = cubit.unite(cp_inpins[" << k-1 <<"][0])" << std::endl;
                m_PyCubGeomFile << "\nsub1.append(assms[" << (k-1)*m_nDimensions <<"])" << std::endl;
                iGeom_subtractEnts(igeomImpl->instance(), tmp_vol,unite, &tmp_new1, &err);
                m_PyCubGeomFile << "tmp_new1 = cubit.subtract(sub2, sub1)" << std::endl;
                m_PyCubGeomFile << "tmp_vol = tmp_new1" << std::endl;

                tmp_vol = tmp_new1;
                unite = NULL;
                tmp_new1=NULL;
              }
            else{ // only one pin in in_pins
                iGeom_subtractEnts(igeomImpl->instance(), tmp_vol, cp_inpins[k-1][0], &tmp_new1, &err);
                m_PyCubGeomFile << "\nsub1.append(assms[" << (k-1)*m_nDimensions <<"])" << std::endl;
                m_PyCubGeomFile << "tmp_new1 = cubit.subtract(sub2, sub1)" << std::endl;
              }
            //#endif
            // This block was needed for OCE below 0.13 or OCC 6.6
            //#if HAVE_OCC
            //            iBase_EntityHandle tmp_new1 = NULL;
            //            // if there are more than one pins
            //            if( cp_inpins[k-1].size() > 1){
            //                std::cout << "Subtraction is slower in OCC, since each pin is subtracted one by one" << std::endl;
            //                for (int i=0; i< (int)cp_inpins[k-1].size(); i++){
            //                    // iGeom_copyEnt(igeomImpl->instance(), cp_inpins[k-1][i], &unite, &err);
            //                    iGeom_subtractEnts(igeomImpl->instance(), tmp_vol,cp_inpins[k-1][i], &tmp_new1, &err);
            //                    ////CHECK("Couldn't subtract pins from block.");
            //                    tmp_vol = tmp_new1;
            //                    tmp_new1=NULL;
            //                  }

            //              }
            //            else{ // only one pin in in_pins
            //                iGeom_subtractEnts(igeomImpl->instance(), tmp_vol, cp_inpins[k-1][0], &tmp_new1, &err);
            //                ////CHECK("Couldn't subtract pins from block.");
            //              }
            //#endif

          }
        std::cout << "\n--------------------------------------------------"<<std::endl;
      }
    else{
        std::cout <<"Nothing to subtract" << std::endl;
      }

  }

  void AssyGen::Imprint_Merge(bool if_merge, bool if_imprint)
  // ---------------------------------------------------------------------------
  // Function: Imprint and Merge
  // Input:    none
  // Output:   none
  // ---------------------------------------------------------------------------
  {
    // getting all entities for merge and imprint
    SimpleArray<iBase_EntityHandle> entities;
    iGeom_getEntities( igeomImpl->instance(), root_set, iBase_REGION, ARRAY_INOUT(entities),&err );
    //CHECK( "ERROR : getRootSet failed!" );

    if(if_imprint ==  true){
        //  now imprint
        std::cout << "\n\nImprinting...." << std::endl;
        clock_t s_imprint = clock();
        iGeom_imprintEnts(igeomImpl->instance(), ARRAY_IN(entities),&err);
        std::cout << "## Imprint CPU time used := " << (double) (clock() - s_imprint)/CLOCKS_PER_SEC
                  << " seconds" << std::endl;
        std::cout << "\n--------------------------------------------------"<<std::endl;
      }

    if(if_merge == true){
        // merge tolerance
        double dTol = 1e-4;
        // now  merge
        std::cout << "\n\nMerging...." << std::endl;
        clock_t s_merge = clock();
        iGeom_mergeEnts(igeomImpl->instance(), ARRAY_IN(entities), dTol, &err);
        std::cout << "## Merge CPU time used := " << (double) (clock() - s_merge)/CLOCKS_PER_SEC
                  << " seconds" << std::endl;
        std::cout <<"merging finished."<< std::endl;
        std::cout << "\n--------------------------------------------------"<<std::endl;
      }
  }

  void AssyGen::Create2DSurf ()
  // ---------------------------------------------------------------------------
  // Function: creating planar top surface with zmax
  // Input:    error code
  // Output:   none
  // ---------------------------------------------------------------------------
  {
    SimpleArray<iBase_EntityHandle>  all_geom;
    SimpleArray<iBase_EntityHandle> surfs;
    int *offset = NULL, offset_alloc = 0, offset_size;
    int t=0;
    std::cout << "Creating surface; 2D assembly specified..." << std::endl;

    // get all the entities in the model (delete after making a copy of top surface)
    iGeom_getEntities( igeomImpl->instance(), root_set, iBase_REGION,ARRAY_INOUT(all_geom),&err );
    //CHECK( "ERROR : Failed to get all geom" );

    // get all the surfaces in the model
    iGeom_getArrAdj( igeomImpl->instance(), ARRAY_IN(all_geom) , iBase_FACE, ARRAY_INOUT(surfs),
                     &offset, &offset_alloc, &offset_size, &err );
    //CHECK( "ERROR : getArrAdj failed!" );

    SimpleArray<double> max_corn, min_corn;
    iGeom_getArrBoundBox( igeomImpl->instance(), ARRAY_IN(surfs), iBase_INTERLEAVED,
                          ARRAY_INOUT( min_corn ),
                          ARRAY_INOUT( max_corn ),
                          &err );
    //CHECK( "Problems getting max surf for rotation." );

    // find the number of surfaces 't' for array allocation
    int nTemp = 1;
    double dTol = 1e-5;
    double dtop = m_dMZAssm(nTemp, 2);
    for (int i = 0; i < surfs.size(); ++i){
        if((fabs(max_corn[3*i+2] -  dtop) < dTol) && (fabs(min_corn[3*i+2] - dtop)<dTol))
          t++;
      }

    // allocate arrays
    SimpleArray<iBase_EntityHandle> max_surfs(t);
    SimpleArray<iBase_EntityHandle> new_surfs(t);
    t=0;

    // store the max surfaces in max_surfs
    for (int i = 0; i < surfs.size(); ++i){

        // locate surfaces for which max and min zcoord is same as maxz coord
        if((fabs(max_corn[3*i+2] -  dtop) < dTol) && (fabs(min_corn[3*i+2] - dtop) < dTol)){
            max_surfs[t] = surfs[i];
            t++;
          }
      }

    // make a copy of max_surfs
    for(int i = 0; i < max_surfs.size(); ++i){
        iGeom_copyEnt(igeomImpl->instance(), max_surfs[i], &new_surfs[i], &err);
        //CHECK( "Problems creating surface." );
      }

    // delete all the old ents
    for(int i=0; i<all_geom.size(); i++){
        iGeom_deleteEnt(igeomImpl->instance(), all_geom[i], &err);
        //CHECK( "Problems deleting cyls." );
      }
    // position the final assembly at the center
    // get the assembly on z=0 plane
    double zcenter = m_dMZAssm(nTemp, 2)/2.0;//move up
    SimpleArray<iBase_EntityHandle> all;
    iGeom_getEntities( igeomImpl->instance(), root_set, iBase_REGION,ARRAY_INOUT(all),&err );
    ////CHECK("Failed to get all entities");

    for(int i=0; i<all.size(); i++){
        iGeom_moveEnt(igeomImpl->instance(),all[i],0,0,-zcenter,&err);
        ////CHECK("Failed to move entities");
      }
    std::cout << "--------------------------------------------------"<<std::endl;

    free(offset);

  }
}
