/*********************************************
Reactor Geometry Generator
Argonne National Laboratory

AssyGen input o/p and functions
*********************************************/
#include <sstream>
#include "nrgen.hpp"
#include "parser.hpp"
#include <fstream>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <iterator>     // std::istream_iterator

#define STRINGIFY_(X) #X
#define STRINGIFY(X) STRINGIFY_(X)

#define DEFAULT_TEST_FILE STRINGIFY(SRCDIR) "/assygen_default"
#define TEST_FILE_NAME "assygen_default"
#define SRC_DIR STRINGIFY(SRCDIR) "/"

// NRGEN CLASS FUNCIONS:

void CNrgen::Banner (std::ostream& OF)
// ---------------------------------------------------------------------------
// Function: Prints program banner
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
  std::cout << '\n';
  std::cout << "\t\t---------------------------------------------------------" << '\n';
  std::cout << "\t\tProgram to Generate Nuclear Reactor Assembly Geometries      " << '\n';
  std::cout << "\t\t\t\tArgonne National Laboratory" << '\n';
  std::cout << "\t\t\t\t        2009-2010         " << '\n';
  std::cout << "\t\t---------------------------------------------------------" << '\n';
  std::cout << "\nsee http://press3.mcs.anl.gov/sigma/meshkit-library/rgg/ for details.\n"<< std::endl;
}


int CNrgen::PrepareIO (int argc, char *argv[])
// ---------------------------------------------------------------------------
// Function: Obtains file names and opens input/output files
// Input:    command line arguments
// Output:   none
// ---------------------------------------------------------------------------
{
  // set and open input output files
  bool bDone = false;
  do{
      if (2 == argc) {
          m_szFile = argv[1];
          m_szInFile=m_szFile+".inp";
          m_szJouFile = m_szFile+".jou";
          m_szSchFile = m_szFile+".template.jou";
          m_szAssmInfo = m_szFile + "_info.csv";
          m_szLogFile = m_szFile + ".log";
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
                        m_szLogFile = m_szFile + ".log";
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
          std::cout << "\nRunning default case:\n" << std::endl;

          m_szInFile = (char *)DEFAULT_TEST_FILE;
          m_szGeomFile = (char *)TEST_FILE_NAME;
          m_szJouFile = (char *)TEST_FILE_NAME;
          m_szFile =  (char *)TEST_FILE_NAME;
          m_szInFile+=".inp";
          m_szJouFile+=".jou";
          m_szSchFile = m_szFile+".template.jou";
          m_szAssmInfo = m_szFile + "_info.csv";
          m_szLogFile = m_szFile + ".log";
          m_szCommonFile = (std::string) SRC_DIR + "common.inp";

          std::cout <<"  No file specified.  Defaulting to: " << m_szInFile
                   << "  " << m_szJouFile << std::endl;
        }
      // open the file
      m_FileInput.open (m_szInFile.c_str(), std::ios::in);
      if (!m_FileInput){
          std::cout << "Unable to open file: " << m_szInFile << std::endl;
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
  //  }o
#endif
  std::cout << "\no/p geometry file name: " <<  m_szGeomFile <<std::endl;

  // writing schemes .jou file ends, now write the main journal file.
  // stuff common to both surface and volume
  m_FileOutput << "## This file is created by rgg program in MeshKit ##\n";
  m_FileOutput << "#User needs to specify mesh interval and schemes in this file\n#" << std::endl;
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
  m_FileOutput << "import step '" << m_szGeomFile <<"'" << std::endl;
#endif

  m_FileOutput << "#" << std::endl;
  return 0;
}


int CNrgen::ReadCommonInp ()
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
  return 0;
}

int CNrgen::ReadInputPhase1 ()
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

  iGeom_newGeom( 0, &geom, &err, 0 ); // this is default way of specifying engine used in configure line
  CHECK("Failed to set geometry engine.");

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
  //	cp_inpins.resize(m_nDuct);
  for (int j=0; j<m_nDuct ; j++)
    cp_inpins.push_back(std::vector<iBase_EntityHandle>());
  return 0;
}

int CNrgen::ReadPinCellData (int i)
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
  return 0;
}


int CNrgen::ReadAndCreate()
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

              err = ReadPinCellData(i);
              ERRORR("Error in ReadPinCellData", err);
              std::cout << "\nread pincell " << i << std::endl;
            }
        }
      if (szInputString.substr(0,8) == "assembly"){
          if(m_szGeomType =="hexagonal"){
              err = Create_HexAssm(szInputString);
              ERRORR("Error in Create_HexAssm", err);
            }
          if(m_szGeomType =="rectangular"){
              err = Create_CartAssm(szInputString);
              ERRORR("Error in Create_CartAssm", err);
            }
          if (m_nJouFlag == 0){
              err = CreateOuterCovering();
              ERRORR("Error in CreateOuterCovering", err);

              // subtract pins before save
              Subtract_Pins();
              if(m_nPlanar ==1){
                  err = Create2DSurf();
                  ERRORR("Error in Create2DSurf", err);
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
          err = Section_Assm(cDir, dOffset, szReverse);
          ERRORR("Error in Section_Assm", err);
          std::cout <<"--------------------------------------------------"<<std::endl;

        }
      if (szInputString.substr(0,4) == "move" && m_nJouFlag == 0){
          std::cout << "Moving geometry .." << std::endl;
          double dX, dY, dZ;
          std::istringstream szFormatString (szInputString);
          szFormatString >> card >> dX >> dY >> dZ;
          if(szFormatString.fail())
            IOErrorHandler(INVALIDINPUT);
          err = Move_Assm(dX, dY, dZ);
          ERRORR("Error in Move_Assm", err);
          std::cout <<"--------------------------------------------------"<<std::endl;

        }
      // ceter the assembly
      if (szInputString.substr(0,6) == "center" && m_nJouFlag == 0){

          char rDir = 'q';
          std::istringstream szFormatString (szInputString);
          szFormatString >> card >> rDir;
          if (rDir != 'q')
            std::cout << "Positioning assembly to "<< rDir << " center" << std::endl;
          else
            std::cout << "Positioning assembly to xy center" << std::endl;
          err = Center_Assm(rDir);
          ERRORR("Error in Center_Assm", err);
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
          err = Rotate_Assm(cDir, dAngle);
          ERRORR("Error in Rotate_Assm", err);
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
              Imprint_Merge();

              // save .sat file
              iGeom_save(geom, m_szGeomFile.c_str(), NULL, &err, m_szGeomFile.length() , 0);
              CHECK("Save to file failed.");
              std::cout << "Normal Termination.\n"<< "Geometry file: " << m_szGeomFile << " saved." << std::endl;
            }
          break;
        }
    }
  return 0;
}

int CNrgen::CreateAssyGenInputFiles()
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


  return 0;
}

int CNrgen::CreateCubitJournal()
//---------------------------------------------------------------------------
//Function: Create Cubit Journal File for generating mesh
//Input:    none
//Output:   none
//---------------------------------------------------------------------------
{
  if(m_szMeshScheme == "hole")
    m_FileOutput << "surf in group hole_surfaces scheme hole" << std::endl;

 if (m_nBLAssemblyMat !=0){
      // Also look for material name in BL material list
      for (int ll=1; ll<= m_nBLAssemblyMat; ll++){
          //if(szVCylMat(m) == m_szBLAssmMat(ll)) {
              m_FileOutput << "group 'tmpgrp' equals surf with name '" <<  m_szBLAssmMat(ll)  << "_top'" << std::endl;
              m_FileOutput << "surf in tmpgrp size {RADIAL_MESH_SIZE}" << std::endl;
              m_FileOutput << "group '" << m_szBLAssmMat(ll) << "_hole_surfaces' equals surf in tmpgrp"<< std::endl;
              m_FileOutput << "surface in group " << m_szBLAssmMat(ll) << "_hole_surfaces scheme hole rad_interval " << m_nBLMatIntervals(ll) << " bias " << m_dBLMatBias(ll) << std::endl;
     //         m_FileOutput << "mesh surf in group " << m_szBLAssmMat(ll) << "_hole_surfaces" << std::endl;
           // }
              m_FileOutput << "group 'bl_surfaces' add surf in tmpgrp" << std::endl; 
        }
    }
  // variables
  int nColor;
  std::string color[21] = {" ", "thistle", "grey", "deepskyblue", "red", "purple",  "green",
                           "yellow", "royalblue", "magenta", "cyan", "lightsalmon", "springgreen",
                           "gold", "orange", "brown", "pink", "khaki", "black", "aquamarine", "mediumslateblue"};

  // if creating only journal file load the geometry file to compute bounding box for automatic size specification
  if(m_nJouFlag == 1){
      iGeom_load(geom, m_szGeomFile.c_str(), NULL, &err, m_szGeomFile.length() , 0);
      CHECK("Failed to load geometry.");
    }

  // get the max and min coordinates of the geometry
  double x1, y1, z1, x2, y2, z2;
  iGeom_getBoundBox( geom, &x1, &y1, &z1, &x2, &y2, &z2, &err );
  CHECK( "Problems getting bouding box." );

  int nSideset=m_nNeumannSetId;
  std::string szGrp, szBlock, szSurfTop, szSurfBot, szSize, szSurfSide;
  double dHeight = 0.0, dMid = 0.0;
  int nTemp = 1;
  if(m_nDimensions > 0){
      dHeight= fabs(z2 - z1);
      dMid = z2 - dHeight/2.0;
    }

  // writing to template.jou
  m_SchemesFile << "## This file is created by rgg program in MeshKit ##\n";
  m_SchemesFile << "##Schemes " << std::endl  ;
  m_SchemesFile << "#{CIRCLE =\"circle interval 1 fraction 0.8\"}" << std::endl;
  m_SchemesFile << "#{HOLE = \"hole rad_interval 2 bias 0.0\"}" << std::endl;
  m_SchemesFile << "#{PAVE = \"pave\"}" << std::endl;
  m_SchemesFile << "#{MAP = \"map\"}" << std::endl;
  m_SchemesFile << "#{SWEEP = \"sweep\"}" << std::endl;
  m_SchemesFile << "#{TET = \"tetmesh\"}" << std::endl;
  m_SchemesFile << "#{TOP_EDGE_INTERVAL = " << m_edgeInterval << " }" << std::endl;
  m_SchemesFile << "## Dimensions" << std::endl;
  if(m_szGeomType == "hexagonal"){
      if(m_nDimensions > 0){
          m_SchemesFile << "#{PITCH =" << m_dMAssmPitch(nTemp, m_nDimensions) << "}" << std::endl;
        }
    }
  else if(m_szGeomType == "rectangular"){
      if(m_nDimensions > 0){
          m_SchemesFile << "#{PITCHX =" << m_dMAssmPitchX(nTemp, m_nDimensions)<< "}" << std::endl;
          m_SchemesFile << "#{PITCHY =" << m_dMAssmPitchY(nTemp, m_nDimensions) << "}" << std::endl;
        }
    }
  if( m_nPlanar ==0){
      m_SchemesFile << "#{Z_HEIGHT = " << dHeight << "}" << std::endl;
      m_SchemesFile << "#{Z_MID = " << dMid << "}" << std::endl;

    }
  m_SchemesFile << "##Set Mesh Sizes" << std::endl;

  if (m_szMeshType == "hex"){
      // volume only
      if(m_nPlanar == 0 ){
          if (m_dAxialSize.GetSize() == 0){
              m_SchemesFile << "#{AXIAL_MESH_SIZE = 0.1*Z_HEIGHT}" << std::endl;
            }
          else {
              m_SchemesFile << "#{AXIAL_MESH_SIZE = " << m_dAxialSize(1) << "}" << std::endl;
            }

          // create templates for specifying block z intervals
          if (m_nDuct > 1){
              m_SchemesFile << "## Set interval along Z direction ## " << std::endl;

              for( int p=1; p<= m_nDuct; p++){
                  if (m_dAxialSize.GetSize() != 0)
                    m_SchemesFile << "#{AXIAL_MESH_SIZE" << p << "=" << m_dAxialSize(p) << "}" << std::endl;
                  else
                    m_SchemesFile << "#{AXIAL_MESH_SIZE" << p << "= 0.1*Z_HEIGHT}" << std::endl;
                  m_SchemesFile << "#{BLOCK" << p << "_Z_INTERVAL = AXIAL_MESH_SIZE" << p << "}" << std::endl;
                  m_SchemesFile << "#{BLOCK" << p << "_ZBOT = " << m_dMZAssm(p, 1) << "}" << std::endl;
                  m_SchemesFile << "#{BLOCK" << p << "_ZTOP = " << m_dMZAssm(p, 2) << "}" << std::endl;
                }
              m_SchemesFile << "##" << std::endl;
            }
        }
      if (-1.0 == m_dRadialSize) {
          if (m_szGeomType == "hexagonal")
            m_SchemesFile << "#{RADIAL_MESH_SIZE = 0.1*PITCH}" << std::endl;
          else
            m_SchemesFile << "#{RADIAL_MESH_SIZE = 0.02*0.5*(PITCHX+PITCHY)}" << std::endl;
        }
      else
        m_SchemesFile << "#{RADIAL_MESH_SIZE = " << m_dRadialSize << "}" << std::endl;
    }
  else if (m_szMeshType == "tet"){
      if (-1.0 == m_dTetMeshSize) {
          if (m_szGeomType == "hexagonal")
            m_SchemesFile << "#{TET_MESH_SIZE = 0.1*PITCH}" << std::endl;
          else
            m_SchemesFile << "#{TET_MESH_SIZE = 0.02*0.5*(PITCHX+PITCHY)}" << std::endl;
        }
      else {
          m_SchemesFile << "#{TET_MESH_SIZE = " << m_dTetMeshSize  << "}" << std::endl;
        }
    }

  if(m_nHblock == -1){ // if more blocks are needed axially, create'em using hexes and the end
      // block creation dumps
      m_FileOutput << "#Creating blocks, Note: you might need to combine some blocks" << std::endl;
      // group creation dumps. each material has a group
      m_FileOutput << "#Creating groups" << std::endl;
      for(int p=1;p<=(m_szAssmMatAlias.GetSize() - m_nBLAssemblyMat);p++){
          szGrp = "g_"+ m_szAssmMat(p);
          m_szAssmMat(p);
          if(m_nPlanar ==1){
              m_FileOutput << "group \"" << szGrp << "\" add surface name \"" << m_szAssmMat(p) <<"\"" << std::endl;
            }
          else{
              m_FileOutput << "group \"" << szGrp << "\" add body name \"" << m_szAssmMat(p) <<"\"" << std::endl;
            }
        }
      for(int p = 1; p <=  (m_szAssmMatAlias.GetSize() - m_nBLAssemblyMat);p++){
          szBlock = "b_"+ m_szAssmMat(p);
          szGrp = "g_"+ m_szAssmMat(p);
          m_FileOutput << "#{nb" << p << " =NumInGrp('" << szGrp << "')}" << std::endl;
          m_FileOutput << "#{Ifndef(nb" << p << ")}" << "\n" << "#{else}" << std::endl;
          if(m_nPlanar ==1){
              m_FileOutput << "block " << m_nMaterialSetId + p -1 << " surface in " << szGrp  << std::endl;
              m_FileOutput << "block " << m_nMaterialSetId + p -1 << " name \"" << szBlock <<"\""<< std::endl;
            }
          else{
              m_FileOutput << "block " << m_nMaterialSetId + p -1 << " body in " << szGrp  << std::endl;
              m_FileOutput << "block " << m_nMaterialSetId + p -1 << " name \"" << szBlock <<"\""<< std::endl;
            }
          m_FileOutput << "#{endif}" << std::endl;
        }
      m_FileOutput << "#" << std::endl;
    }
  if(m_szMeshType == "hex"){
      // imprint
      m_FileOutput << "#Imprint geometry" << std::endl;
      m_FileOutput << "imprint all" << std::endl;
      m_FileOutput << "#" << std::endl;
      // merge

      m_FileOutput << "Merge Tolerance " << m_dMergeTol << std::endl;
      m_FileOutput << "#" << std::endl;

      m_FileOutput << "#Merge geometry" << std::endl;
      m_FileOutput << "merge all" << std::endl;
      m_FileOutput << "#" << std::endl;
    }

  // for info keyword
  if(strcmp(m_szInfo.c_str(),"on") == 0){
      int temp = 9700;
      m_FileOutput << "# stuff for info keyword, remove if not desired " << std::endl;
      m_FileOutput << "# putting pins in seperate blocks " << std::endl;
      m_FileOutput << "#" << std::endl;
      for (int i=0; i<m_nTotalPincells; i++){
          m_FileOutput << "group 'g"<< i+m_nStartpinid << "' add body with name '_xp" << i+m_nStartpinid << "_'" << std::endl;

          m_FileOutput << "#{nbody" << i+1 << " =NumInGrp('g" <<i+m_nStartpinid << "')}" << std::endl;
          m_FileOutput << "#{Ifndef(nbody" << i+1 << ")}" << "\n" << "#{else}" << std::endl;
          m_FileOutput << "block " << temp+i << " body in group g" << i+m_nStartpinid << std::endl;
          m_FileOutput << "block " << temp+i << " name '_xp" << i+m_nStartpinid << "'" << std::endl;
          m_FileOutput << "#{endif}" << std::endl;
        }
    }

  //surface only
  if(m_nPlanar ==1){
      m_FileOutput << "# Pointing surface normals to 0.0, 0.0, -1.0 or -ve Z or correct STARCCM+ cell-face orientation" << std::endl;
      m_FileOutput << "surface all normal opposite" << std::endl;
      m_FileOutput << "#" << std::endl;
    }
  // volume only
  else{
      if(m_szSideset == "yes"){

          // rename the skin surfaces, so that they don't appear as sidesets
          for (int p=1; p<=m_nDuct; p++){
              for(int q=1;q<=m_nSides; q++){
                  m_FileOutput << "group 'edge" << (m_nSides*(p-1) + q ) <<"' equals curve with name 'side_edge"
                               << (m_nSides*(p-1) + q ) << "@'" << std::endl;

                  m_FileOutput << "group 'vt" <<  (m_nSides*(p-1) + q )  <<"' equals vertex with z_max == z_min in curve in edge"
                               <<  (m_nSides*(p-1) + q ) << std::endl;

                }
            }

          // creating groups for vertices on the top surface of the duct
          for (int p=1; p<=m_nDuct; p++){
              for(int q=1;q<=m_nSides; q++){

                  if(q != m_nSides){
                      m_FileOutput << "group 'v" << (m_nSides*(p-1) + q ) <<"' intersect group vt" << (m_nSides*(p-1) + q )
                                   << " with group vt" <<  (m_nSides*(p-1) + q + 1 )  << std::endl;
                    }
                  else {
                      m_FileOutput << "group 'v" << (m_nSides*(p-1) + q ) <<"' intersect group vt" << (m_nSides*(p-1) + q )
                                   << " with group vt" <<  (m_nSides*(p-1) + 1 )  << std::endl;
                    }
                }
            }
          // creating temp surfaces groups
          for (int p=1; p<=m_nDuct; p++){
              for(int q=1;q<=m_nSides; q++){
                  m_FileOutput << "group 'st" << (m_nSides*(p-1) + q ) <<"' equals surface with z_max <> z_min in vert in v"
                               << (m_nSides*(p-1) + q ) << "'" << std::endl;
                }
            }

          // creating surface groups for obtaining surfaces
          for (int p=1; p<=m_nDuct; p++){
              for(int q=1;q<=m_nSides; q++){
                  if(q != 1){
                      m_FileOutput << "group 's" << (m_nSides*(p-1) + q ) <<"' intersect group st"  << (m_nSides*(p-1) + q )
                                   << " with group st" <<  (m_nSides*(p-1) + q - 1 )  << std::endl;
                    }
                  else {
                      m_FileOutput << "group 's" << (m_nSides*(p-1) + q ) <<"' intersect group st" << (m_nSides*(p-1) + q )
                                   << " with group st" <<  (m_nSides*(p-1) + m_nSides )  << std::endl;
                    }
                }
            }

          // renaming the skin side surfaces
          for (int p=1; p<=m_nDuct; p++){
              for(int q=1;q<=m_nSides; q++){
                  m_FileOutput << "surface in group s" <<  (m_nSides*(p-1) + q ) << " rename 'side_surface"
                               <<  (m_nSides*(p-1) + q ) << "'" << std::endl;

                }
            }
        }

      if(m_szMeshType == "hex"){

          //now set the sizes
          m_FileOutput << "#Set Meshing Scheme and Sizes, use template.jou to specify sizes" << std::endl;

          for(int p=1;p<=(m_szAssmMatAlias.GetSize() - m_nBLAssemblyMat);p++){
              szGrp = "g_"+ m_szAssmMat(p);
              szSize =  m_szAssmMat(p) + "_size";
              szSurfBot = m_szAssmMat(p) + "_bot";
              szSize =  m_szAssmMat(p) + "_surf_size";
              m_FileOutput << "group 'tmpgrp' equals surface name \""  << szSurfBot  << "\"" << std::endl;
              m_FileOutput << "surface in tmpgrp  size {"  << szSize <<"}" << std::endl;
            }
          m_FileOutput << "#" << std::endl;
        }
    }

  if(m_szMeshType == "hex"){
      // some more common stuff meshing top surfaces set the sizes and mesh
      m_FileOutput << "#Surfaces mesh, use template.jou to specify sizes" << std::endl;
      for(int p=1;p<=(m_szAssmMatAlias.GetSize() - m_nBLAssemblyMat);p++){
          szSurfTop = m_szAssmMat(p) + "_top";
          szGrp = "g_"+ m_szAssmMat(p);
          szSize =  m_szAssmMat(p) + "_surf_size";
          if(m_szMeshScheme == "hole"){
              m_FileOutput << "group 'tmpgrp' equals surface name \""  << szSurfTop  << "\"" << std::endl;
              m_FileOutput << "group 'remove_hole' intersect group tmpgrp with group hole_surfaces" << std::endl;
              m_FileOutput << "#{nIntersect=NumInGrp('remove_hole')}" << std::endl;
              m_FileOutput << "#{If(nIntersect==0)}" << std::endl;
              m_FileOutput << "surface in tmpgrp  size {"  << szSize <<"}" << std::endl;
              m_FileOutput << "surface in tmpgrp scheme {" << "PAVE" << "}"  << std::endl;
              m_FileOutput << "#{endif}"  << std::endl;
            }
          else{
              m_FileOutput << "group 'tmpgrp' equals surface name \""  << szSurfTop  << "\"" << std::endl;
              m_FileOutput << "surface in tmpgrp  size {"  << szSize <<"}" << std::endl;
              m_FileOutput << "surface in tmpgrp scheme {" << "PAVE" << "}"  << std::endl;
            }

          if (p==1 && m_edgeInterval != 99){
              m_FileOutput << "group 'edges" <<"' equals curve with name 'side_edge'"<< std::endl;
              m_FileOutput << "curve in edges interval {TOP_EDGE_INTERVAL}" << std::endl;
            }

          //    m_FileOutput << "mesh surface in " << szGrp << "\n#" << std::endl;

          // dumping these sizes schemes.jou also
          m_SchemesFile << "#{"  << szSize <<" = RADIAL_MESH_SIZE}" << std::endl;
        }
      m_FileOutput << "#" << std::endl;

      // mesh all command after meshing surface
      if (m_nDuct <= 1 ){
          m_FileOutput << "group 'tmpgrp' add surface name '_top'" << std::endl;
          m_FileOutput << "group 'tmpgrp1' subtract innerduct from tmpgrp" << std::endl;
          m_FileOutput << "group 'tmpgrp2' subtract bl_surfaces from tmpgrp1" << std::endl;
          m_FileOutput << "mesh tmpgrp2" << std::endl;
        }
      else {
          m_FileOutput << "#Meshing top surface" << std::endl;
          //m_FileOutput << "mesh surface with z_coord = " << z2 << std::endl;
          if(m_szMeshScheme == "hole"){
              m_FileOutput << "group 'tmpgrp' equals surface name \""  << szSurfTop  << "\"" << std::endl;
              m_FileOutput << "group 'remove_hole' intersect group tmpgrp with group hole_surfaces" << std::endl;
              m_FileOutput << "#{nIntersect=NumInGrp('remove_hole')}" << std::endl;
              m_FileOutput << "#{If(nIntersect==0)}" << std::endl;
              m_FileOutput << "surface in tmpgrp  size {"  << szSize <<"}" << std::endl;
              m_FileOutput << "surface in tmpgrp scheme {" << "PAVE" << "}"  << std::endl;
              m_FileOutput << "#{endif}"  << std::endl;
              m_FileOutput << "mesh surface with z_coord = " << z2 << std::endl;
            }
          else{
              m_FileOutput << "group 'tmpgrp' equals surface name \""  << szSurfTop  << "\"" << std::endl;
              m_FileOutput << "surface in tmpgrp  size {"  << szSize <<"}" << std::endl;
              m_FileOutput << "surface in tmpgrp scheme {" << "PAVE" << "}"  << std::endl;
              m_FileOutput << "mesh surface with z_coord = " << z2 << std::endl;
            }
        }
 if (m_nBLAssemblyMat !=0){
      // Also look for material name in BL material list
      for (int ll=1; ll<= m_nBLAssemblyMat; ll++){
              m_FileOutput << "mesh surf in group " << m_szBLAssmMat(ll) << "_hole_surfaces" << std::endl;
        }
      m_FileOutput << "mesh surf in innerduct" << std::endl;
    }

      if(m_nPlanar == 0){ // volumes only
          if (m_nDuct == 1){
              m_FileOutput << "surf with z_coord > {Z_MID -.1*Z_HEIGHT}" <<
                              " and z_coord < {Z_MID + .1*Z_HEIGHT} size {AXIAL_MESH_SIZE}" << std::endl ;
              m_FileOutput << "mesh vol all" << std::endl;
            }
          else if (m_nDuct > 1){
              m_FileOutput << "### Setting Z intervals on ducts and meshing along Z " << std::endl;
              for( int p=m_nDuct; p>= 1; p--){
                  if(dMid == 0){ // z - centered
                      m_FileOutput << "surf with z_coord  > " << m_dMZAssm(p, 1) - dHeight/2.0
                                   << " and z_coord < " << m_dMZAssm(p, 2) - dHeight/2.0 << " interval " << "{BLOCK" << p << "_Z_INTERVAL}" << std::endl;
                      m_FileOutput << "mesh vol with z_coord  > " << m_dMZAssm(p, 1) - dHeight/2.0
                                   << " and z_coord < " << m_dMZAssm(p, 2) - dHeight/2.0 << std::endl;
                    }
                  else{
                      m_FileOutput << "surf with z_coord  > " << m_dMZAssm(p, 1)
                                   << " and z_coord < " << m_dMZAssm(p, 2) << " interval " << "{BLOCK" << p << "_Z_INTERVAL}" << std::endl;
                      m_FileOutput << "mesh vol with z_coord  > " << m_dMZAssm(p, 1)
                                   << " and z_coord < " << m_dMZAssm(p, 2) << std::endl;

                      m_FileOutput << "##" << std::endl;
                    }
                }
            }
        }
    }
  else if(m_szMeshType == "tet"){
      m_FileOutput << "##"<< std::endl;
      m_FileOutput << "# groupings for creating vertex groups"<< std::endl;
      for (int p=1; p<=m_nDuct; p++){
          for(int q=1;q<=m_nSides; q++){
              m_FileOutput << "group 'edge" << (m_nSides*(p-1) + q ) <<"' equals curve with name 'side_edge"
                           << (m_nSides*(p-1) + q ) << "@'" << std::endl;

              m_FileOutput << "group 'vt" <<  (m_nSides*(p-1) + q )  <<"' equals vertex with z_max == z_min in curve in edge"
                           <<  (m_nSides*(p-1) + q ) << std::endl;

            }
        }

      // creating groups for vertices on the top surface of the duct
      for (int p=1; p<=m_nDuct; p++){
          for(int q=1;q<=m_nSides; q++){

              if(q != m_nSides){
                  m_FileOutput << "group 'v" << (m_nSides*(p-1) + q ) <<"' intersect group vt" << (m_nSides*(p-1) + q )
                               << " with group vt" <<  (m_nSides*(p-1) + q + 1 )  << std::endl;
                }
              else {
                  m_FileOutput << "group 'v" << (m_nSides*(p-1) + q ) <<"' intersect group vt" << (m_nSides*(p-1) + q )
                               << " with group vt" <<  (m_nSides*(p-1) + 1 )  << std::endl;
                }
            }
        }

      // creating temp surfaces groups
      for (int p=1; p<=m_nDuct; p++){
          for(int q=1;q<=m_nSides; q++){
              m_FileOutput << "group 'st" << (m_nSides*(p-1) + q ) <<"' equals surface with z_max <> z_min in vert in v"
                           << (m_nSides*(p-1) + q ) << "'" << std::endl;
            }
        }

      // creating side curve and surface groups
      for (int p=1; p<=m_nDuct; p++){
          for(int q=1;q<=m_nSides; q++){

              m_FileOutput << "group 'c" <<  (m_nSides*(p-1) + q )  <<"' equals curve with z_max <> z_min in vert in v"
                           <<  (m_nSides*(p-1) + q ) << std::endl;

              if(q != 1){
                  m_FileOutput << "group 's" << (m_nSides*(p-1) + q ) <<"' intersect group st"  << (m_nSides*(p-1) + q )
                               << " with group st" <<  (m_nSides*(p-1) + q - 1 )  << std::endl;
                }
              else {
                  m_FileOutput << "group 's" << (m_nSides*(p-1) + q ) <<"' intersect group st" << (m_nSides*(p-1) + q )
                               << " with group st" <<  (m_nSides*(p-1) + m_nSides )  << std::endl;
                }
            }
        }

      // renaming the side surfaces for getting the split surfaces later
      for (int p=1; p<=m_nDuct; p++){
          for(int q=1;q<=m_nSides; q++){
              m_FileOutput << "surface in group s" <<  (m_nSides*(p-1) + q ) << " rename 'side_surface"
                           <<  (m_nSides*(p-1) + q ) << "'" << std::endl;
            }
        }

      // splitting the surfaces
      for (int p=1; p<=m_nDuct; p++){
          for(int q=1;q<=m_nSides; q++){

              m_FileOutput << "split surface in group s" <<  (m_nSides*(p-1) + q )  <<" direction curve in group c"
                           <<  (m_nSides*(p-1) + q ) << std::endl;
            }
        }

      // get all the split surfaces in individual groups
      for (int p=1; p<=m_nDuct; p++){
          for(int q=1;q<=m_nSides; q++){

              m_FileOutput << "group 'sname" << (m_nSides*(p-1) + q ) <<  "' equals surface with name 'side_surface"
                           <<  (m_nSides*(p-1) + q ) << "'"<< std::endl;
              m_FileOutput << "group 'svert" << (m_nSides*(p-1) + q ) <<  "' equals surface in vert in v"
                           <<  (m_nSides*(p-1) + q ) << std::endl;
              m_FileOutput << "group 'ssplit" << (m_nSides*(p-1) + q ) <<  "' intersect group sname" <<  (m_nSides*(p-1) + q )
                           << " with group svert" << (m_nSides*(p-1) + q ) << std::endl;
            }
        }

      // get all the split surfaces in individual groups
      for (int p=1; p<=m_nDuct; p++){
          for(int q=1;q<=m_nSides; q++){

              if(q != 1){
                  m_FileOutput << "group 'ssplit_" << (m_nSides*(p-1) + q ) <<"' intersect group sname"  << (m_nSides*(p-1) + q )
                               << " with group svert" <<  (m_nSides*(p-1) + q - 1 )  << std::endl;
                }
              else {
                  m_FileOutput << "group 'ssplit_" <<  (m_nSides*(p-1) + q ) <<"' intersect group sname" << (m_nSides*(p-1) + q )
                               << " with group svert" <<  (m_nSides*(p-1) + m_nSides )  << std::endl;
                }
            }
        }
      // imprint
      m_FileOutput << "#Imprint geometry" << std::endl;
      m_FileOutput << "imprint all" << std::endl;
      m_FileOutput << "#" << std::endl;
      m_FileOutput << "Merge Tolerance " << m_dMergeTol << std::endl;
      m_FileOutput << "#" << std::endl;

      // merge
      m_FileOutput << "#Merge geometry" << std::endl;
      m_FileOutput << "merge all" << std::endl;
      m_FileOutput << "#" << std::endl;

      m_FileOutput << "#Set mesh scheme and size" << std::endl;
      m_FileOutput << "volume all scheme {TET} size {TET_MESH_SIZE}" << std::endl;

      // mesh one side of each duct, such that one is flipped mesh of the other
      for (int p=1; p<=m_nDuct; p++){

          m_FileOutput << "mesh surface in group ssplit" <<  (m_nSides*(p-1) + 1) << std::endl;

          m_FileOutput << "surface in group ssplit_" << (m_nSides*(p-1) + 1) << " scheme copy source surface in group ssplit"
                       <<  (m_nSides*(p-1) + 1)
                        << " source curve in group c" <<  (m_nSides*(p-1) + 1 ) << " target curve in group c" <<  (m_nSides*(p-1) + m_nSides )
                        << " source vertex in group v" <<  (m_nSides*(p-1) + 1) << " target vertex in group v"  <<  (m_nSides*(p-1) + m_nSides )
                        << " nosmoothing" << std::endl;

          m_FileOutput << "mesh surface in group  ssplit_" << (m_nSides*(p-1) + 1) << std::endl;
        }

      // setting the copy mesh commands on the above pair of split surfaces to have all surfaces symmetrical
      for (int p=1; p<=m_nDuct; p++){
          for(int q=1;q<=m_nSides; q++){
              if(q != m_nSides){
                  m_FileOutput << "copy mesh surface in ssplit" <<  (m_nSides*(p-1) + 1)
                               << " onto surface in ssplit" << (m_nSides*(p-1) + q + 1 )
                               << " source vertex in group v" << (m_nSides*(p-1) + 1)
                               << " target vertex in group v" <<  (m_nSides*(p-1) + q + 1) << " nosmoothing" <<  std::endl;

                  m_FileOutput << "copy mesh surface in ssplit_" << (m_nSides*(p-1) + 1 )
                               << " onto surface in ssplit_" << (m_nSides*(p-1) + q +1 )
                               << " source vertex in group v" << (m_nSides*p)
                               << " target vertex in group v" <<  (m_nSides*(p-1) + q) << " nosmoothing" << std::endl;

                }
              else{
                  // do nothing
                }
            }
        }

      m_FileOutput << "# Mesh all volumes now" << std::endl;
      m_FileOutput << "mesh vol all" << std::endl;
    }

  // create and sidesets after meshing
  m_FileOutput << "#" << std::endl;
  //    }
  if(m_szSideset == "yes"){
      // top surface sidesets
      m_FileOutput << "#Creating top surface sidesets" << std::endl;
      m_FileOutput << "create group 'surfall'" << std::endl;
      for(int p=1;p<=(m_szAssmMatAlias.GetSize() - m_nBLAssemblyMat);p++){
          ++nSideset;
          szSurfTop = m_szAssmMat(p)+"_top";
          // Avoid creation if empty sideset
          m_FileOutput << "group 'tmpgrp' equals surface with name '" << szSurfTop << "' in vol in block " << m_nMaterialSetId + p -1 << std::endl;
          m_FileOutput << "sideset " << nSideset << " surface in tmpgrp " << std::endl;
          m_FileOutput << "sideset " << nSideset << " name \"" << szSurfTop << "_ss\"" << std::endl;
        }
      m_FileOutput << "#" << std::endl;
      for(int p=1;p<=m_nBLAssemblyMat;p++){
          ++nSideset;
          szSurfTop = m_szBLAssmMat(p)+"_top";
          // Avoid creation if empty sideset
          m_FileOutput << "group 'tmpgrp' equals surface with name '" << szSurfTop << std::endl;
          m_FileOutput << "sideset " << nSideset << " surface in tmpgrp " << std::endl;
          m_FileOutput << "sideset " << nSideset << " name \"" << szSurfTop << "_ss\"" << std::endl;
        }
      m_FileOutput << "#" << std::endl;
    }


  if(m_nPlanar ==0){
      if(m_szSideset == "yes"){
          // now create bot and side sideset
          m_FileOutput << "#Creating bot/side surface sidesets" << std::endl;
          for(int p=1;p<=(m_szAssmMatAlias.GetSize() - m_nBLAssemblyMat);p++){
              szSurfTop = m_szAssmMat(p)+"_bot";
              m_FileOutput << "#" << std::endl;
              ++nSideset;
              m_FileOutput << "group 'tmpgrp' equals surface with name '" << szSurfTop << "' in vol in block " << m_nMaterialSetId + p -1 << std::endl;
              m_FileOutput << "sideset " << nSideset << " surface in tmpgrp " << std::endl;
              m_FileOutput << "sideset " << nSideset << " name \"" << szSurfTop << "_ss\"" << std::endl;
            }
          for(int p=1;p<=m_nBLAssemblyMat;p++){
              ++nSideset;
              szSurfTop = m_szBLAssmMat(p)+"_bot";
              // Avoid creation if empty sideset
              m_FileOutput << "group 'tmpgrp' equals surface with name '" << szSurfTop << std::endl;
              m_FileOutput << "sideset " << nSideset << " surface in tmpgrp " << std::endl;
              m_FileOutput << "sideset " << nSideset << " name \"" << szSurfTop << "_ss\"" << std::endl;
            }
          m_FileOutput << "#" << std::endl;

          for(int p=1;p<=(m_szAssmMatAlias.GetSize() - m_nBLAssemblyMat);p++){
              szSurfSide = m_szAssmMat(p)+"_side";
              ++nSideset;
              if(m_szGeomType == "hexagonal"){
                  for (int u=1; u<=6;u++){
                      m_FileOutput << "group 'tmpgrp" << u <<"' equals surf with name '" << szSurfSide << u << "'" << std::endl;
                    }
                  m_FileOutput << "sideset " << nSideset << " surface in tmpgrp1 tmpgrp2 tmpgrp3 tmpgrp4 tmpgrp5 tmpgrp6'" << std::endl;
                  m_FileOutput << "sideset " << nSideset << " name \"" << szSurfSide << "1_ss\"" << std::endl;
                  ++nSideset;
                }
              if(m_szGeomType == "hexagonal"){
                  for (int u=7; u<=12;u++){
                      m_FileOutput << "group 'tmpgrp" << u <<"' equals surf with name '" << szSurfSide << "_" << u << "'" << std::endl;
                    }
                  m_FileOutput << "sideset " << nSideset << " surface in tmpgrp7 tmpgrp8 tmpgrp9 tmpgrp10 tmpgrp11 tmpgrp12'" << std::endl;
                  m_FileOutput << "sideset " << nSideset << " name \"" << szSurfSide << "2_ss\"" << std::endl;
                }
              if(m_szGeomType == "rectangular"){
                  for (int u=1; u<=4;u++){
                      m_FileOutput << "group 'tmpgrp" << u <<"' equals surf with name '" << szSurfSide << u << "'" << std::endl;
                    }
                  m_FileOutput << "sideset " << nSideset << " surface in tmpgrp1 tmpgrp2 tmpgrp3 tmpgrp4'" << std::endl;
                  m_FileOutput << "sideset " << nSideset << " name \"" << szSurfSide << "2_ss\"" << std::endl;
                  ++nSideset;
                }

              if(m_szGeomType == "rectangular"){
                  for (int u=5; u<=8;u++){
                      m_FileOutput << "group 'tmpgrp" << u <<"' equals surf with name '" << szSurfSide  << "_" << u << "'" << std::endl;
                    }
                  m_FileOutput << "sideset " << nSideset << " surface in tmpgrp5 tmpgrp6 tmpgrp7 tmpgrp8'" << std::endl;
                  m_FileOutput << "sideset " << nSideset << " name \"" << szSurfSide << "2_ss\"" << std::endl;
                }
            }


          m_FileOutput << "#" << std::endl;

          m_FileOutput << "#Creating sideset for outer most side surfaces" << std::endl;
          ++nSideset;

          m_FileOutput << "group 'tmpgrp' equals surf with name 'side_surface'" << std::endl;
          m_FileOutput << "sideset " << nSideset << " surface in tmpgrp " << std::endl;
          m_FileOutput << "sideset " << nSideset << " name \"" << "outer_side_ss\"" << std::endl;
        }
    }
  if(m_nHblock != -1){ // if more blocks are needed axially, create'em using hexes and the end
      // block creation dumps
      m_FileOutput << "#Creating blocks, Note: you might need to combine some blocks" << std::endl;
      // group creation dumps. each material has a group
      m_FileOutput << "#Creating groups" << std::endl;
      if(m_szMeshType == "hex"){
          for(int p=1;p<=(m_szAssmMatAlias.GetSize() - m_nBLAssemblyMat);p++){
              szGrp = "g_"+ m_szAssmMat(p);
              m_szAssmMat(p);
              if(m_nPlanar ==1){
                  m_FileOutput << "group \"" << szGrp << "\" add surface name \"" << m_szAssmMat(p) <<"\"" << std::endl;
                }
              else{

                  m_FileOutput << "group \"" << szGrp << "\" add body name \"" << m_szAssmMat(p) <<"\"" << std::endl;
                }
            }
          for(int p = 1; p <=  (m_szAssmMatAlias.GetSize() - m_nBLAssemblyMat);p++){
              szBlock = "b_"+ m_szAssmMat(p);
              szGrp = "g_"+ m_szAssmMat(p);
              m_FileOutput << "#{nb" << p << " =NumInGrp('" << szGrp << "')}" << std::endl;
              m_FileOutput << "#{Ifndef(nb" << p << ")}" << "\n" << "#{else}" << std::endl;
              if(m_nPlanar ==1){
                  m_FileOutput << "block " << m_nMaterialSetId + p -1 << " surface in " << szGrp  << std::endl;
                  m_FileOutput << "block " << m_nMaterialSetId + p -1 << " name \"" << szBlock <<"\""<< std::endl;
                }
              else{
                  m_FileOutput << "block " << m_nMaterialSetId + p -1 << " hex in body in " << szGrp  << std::endl;
                  m_FileOutput << "block " << m_nMaterialSetId + p -1 << " name \"" << szBlock <<"\""<< std::endl;
                }
              m_FileOutput << "#{endif}" << std::endl;
            }
          m_FileOutput << "#" << std::endl;
        }
      else{
          std::cout << "Error: Terminating journal file writing. \n Hex block (Hblock keyword) is not supported for a tet mesh." << std::endl;
          exit(1);
        }
    }

  // create super blocks
  if(m_nSuperBlocks > 0){
      for(int o = 1; o <= m_nSuperBlocks; o++){
          m_FileOutput << "block " << sb(o).m_nSuperBlockId << " vol in block ";
          for (int p = 1; p <= sb(o).m_nNumSBContents; p++){
              m_FileOutput << m_nMaterialSetId + sb(o).m_nSBContents(p) << " ";
            }
          m_FileOutput << "\n" << "block " << sb(o).m_nSuperBlockId << " name '" << sb(o).m_szSuperBlockAlias << "'" << std::endl;
          m_FileOutput << "delete block " ;
          for (int q = 1; q <= sb(o).m_nNumSBContents; q++){
              m_FileOutput << m_nMaterialSetId + sb(o).m_nSBContents(q) << " ";
            }
          m_FileOutput << "\n" << std::endl;
        }
    }



  if(m_nHblock > 0){
      // now dump the commands for making hex layers as blocks and subtracting from original
      double delta = (m_dZend - m_dZstart)/m_nHblock;
      for(int i=0; i<(m_szAssmMatAlias.GetSize() - m_nBLAssemblyMat); i++){
          m_FileOutput << "## BLOCK CREATION USING HEXES" << std::endl;
          for(int j=0; j<m_nHblock; j++){
              m_FileOutput << "group 'tmpgrp" << j+1 << "' equals hex in block " <<  m_nMaterialSetId + i
                           << " with z_coord < " << m_dZstart + (j+1)*delta << " and z_coord > "
                           << m_dZstart + j*delta << std::endl;
            }
          for(int j=0; j<m_nHblock; j++){
              m_FileOutput << "block " <<  m_nMaterialSetId+i << " group tmpgrp" << j+1 << " remove" << std::endl;
            }
          for(int j=0; j<m_nHblock; j++){
              if((m_szAssmMatAlias.GetSize() - m_nBLAssemblyMat) < 10)
                m_FileOutput << "block " << j+1 <<  m_nMaterialSetId+i << " group tmpgrp" << j+1 << std::endl;
              else
                m_FileOutput << "block " << (j+1)*10 <<  m_nMaterialSetId+i << " group tmpgrp" << j+1 << std::endl;
            }
        }
    }
  if(m_nMaterialSetId != 1)
    m_FileOutput << "renumber hex all start_id " << MAXLINES*1000 << std::endl;
  // color now
  m_FileOutput << "#Set color for different parts" << std::endl;
  if(m_nPlanar == 0){ // volumes only
      for(int p=1;p<=(m_szAssmMatAlias.GetSize() - m_nBLAssemblyMat);p++){
          szGrp = "g_"+ m_szAssmMat(p);
          if(p>20)
            nColor = 1;
          else
            nColor = p;
          m_FileOutput << "color body in " << szGrp << " " << color[nColor] << std::endl;
        }
    }
  else{ //surfaces
      // color now
      for(int p=1;p<=(m_szAssmMatAlias.GetSize() - m_nBLAssemblyMat);p++){
          szGrp = "g_"+ m_szAssmMat(p);
          if(p>20)
            nColor = 1;
          else
            nColor = p;
          m_FileOutput << "color surface in " << szGrp << " " << color[nColor] << std::endl;
        }
    }

  m_FileOutput << "delete group all" << std::endl;
  // save as .cub file dump
  m_FileOutput << "#\n#Save file" << std::endl;
  if(save_exodus){
      std::string szSave = m_szFile + ".exo";
      std::transform(szSave.begin(), szSave.end(), szSave.begin(), ::tolower);
      m_FileOutput << "export mesh '"<< szSave <<"'" << " overwrite"<<std::endl;
    }
  else{
      std::string szSave = m_szFile + ".cub";
      std::transform(szSave.begin(), szSave.end(), szSave.begin(), ::tolower);
      m_FileOutput << "save as '"<< szSave <<"'" << " overwrite"<<std::endl;
    }

  std::cout << "Schemes file created: " << m_szSchFile << std::endl;
  std::cout << "Cubit journal file created: " << m_szJouFile << std::endl;
  if(strcmp(m_szInfo.c_str(),"on") == 0)
    std::cout << "Assembly info file created: " << m_szAssmInfo << std::endl;

  m_FileOutput << "Timer Stop" << std::endl;
  return 0;
}

int CNrgen:: ComputePinCentroid(int nTempPin, CMatrix<std::string> MAssembly, 
                                int m, int n, double &dX, double &dY, double &dZ)
// ---------------------------------------------------------------------------
// Function: computes the centroid in the whole assembly of rectangular or hexagonal pincell
// Input:    number and location of the pincell
// Output:   coordinates of pin in assembly
// ---------------------------------------------------------------------------
{
  int nTempPin1 = 0, nTempPin2 = 0, nInputLines;
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
  return 0;
}

void CNrgen::IOErrorHandler (ErrorStates ECode) const
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

int CNrgen::TerminateProgram ()
// ---------------------------------------------------------------------------
// Function: terminates the program steps by closing the input/output files
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
  iGeom_dtor(geom, &err);
  CHECK( "Interface destruction didn't work properly." );
  // close the input and output files
  m_FileInput.close ();
  m_FileOutput.close ();
  m_SchemesFile.close ();
  if(strcmp(m_szInfo.c_str(),"on") == 0)
    m_AssmInfo.close ();

  return 0;
}


// print error function definition (iGeom)
bool CNrgen::Print_Error( const char* desc, 
                          int err,
                          iGeom_Instance geom,
                          const char* file,
                          int line )
{
  char buffer[1024];
  iGeom_getDescription( geom, buffer, sizeof(buffer) );
  buffer[sizeof(buffer)-1] = '\0';

  std::cerr << "ERROR: " << desc << std::endl
            << "  Error code: " << err << std::endl
            << "  Error desc: " << buffer << std::endl
            << "  At        : " << file << ':' << line << std::endl
               ;

  return false; // must always return false or CHECK macro will break
}

int CNrgen:: Name_Faces(const std::string sMatName, const iBase_EntityHandle body,  iBase_TagHandle this_tag )
// ---------------------------------------------------------------------------
// Function: names all the faces in the body
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
  double dTol = 1e-4, ttol = 0.0;
  if(m_szGeomType == "hexagonal")
    ttol = m_dMAssmPitch(1, 1);
  else if (m_szGeomType == "rectangular")
    ttol = m_dMAssmPitchX(1,1);
  // set tolerance for surface identification
  if (ttol != 0){
      dTol=ttol*1.0e-4;
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
  iGeom_getEntAdj( geom, body, iBase_FACE, ARRAY_INOUT(surfs), &err );
  CHECK( "Problems getting max surf for rotation." );

  SimpleArray<double> max_corn, min_corn;
  iGeom_getArrBoundBox( geom, ARRAY_IN(surfs), iBase_INTERLEAVED,
                        ARRAY_INOUT( min_corn ),
                        ARRAY_INOUT( max_corn ),
                        &err );
  CHECK( "Problems getting max surf for rotation." );
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
          iGeom_setData(geom, max_surf, this_tag,
                        sMatName0.c_str(), sMatName0.size(), &err);
          CHECK("setData failed");

          std::cout << sMatName0 << ",  ";
          max_surf = NULL;

        }
      if(min_surf !=0){
          iGeom_setData(geom, min_surf, this_tag,
                        sMatName1.c_str(), sMatName1.size(), &err);
          CHECK("setData failed");
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
          iGeom_setData(geom, side_surf, this_tag,
                        sSideName.c_str(), sMatName1.size(), &err);
          CHECK("setData failed");
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
  return 0;
}


int CNrgen::Center_Assm (char &rDir)
// ---------------------------------------------------------------------------
// Function: centers all the entities along x and y axis
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
  double xmin, xmax, ymin, ymax, zmin, zmax, xcenter = 0.0, ycenter = 0.0, zcenter = 0.0;
  // position the assembly such that origin is at the center before sa
  iGeom_getBoundBox(geom,&xmin,&ymin,&zmin,
                    &xmax,&ymax,&zmax, &err);
  CHECK("Failed getting bounding box");

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
  iGeom_getEntities( geom, root_set, iBase_REGION,ARRAY_INOUT(all),&err );
  CHECK("Failed to get all entities");

  for(int i=0; i<all.size(); i++){
      iGeom_moveEnt(geom,all[i],-xcenter,-ycenter,-zcenter,&err);
      CHECK("Failed to move entities");
    }
  return 0;
}

int CNrgen::Section_Assm (char &cDir, double &dOffset, const std::string szReverse)
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
  iGeom_getEntities( geom, root_set, iBase_REGION,ARRAY_INOUT(all),&err );
  CHECK("Failed to get all entities");
  // loop and section/delete entities
  for(int i=0; i < all.size(); i++){
      //get the bounding box to decide
      iGeom_getEntBoundBox(geom,all[i],&xmin,&ymin,&zmin,
                           &xmax,&ymax,&zmax, &err);
      CHECK("Failed get bound box");
      if(xmin > dOffset && yzplane ==1 && nReverse ==1){
          iGeom_deleteEnt(geom,all[i],&err);
          CHECK("Failed delete entities");
          continue;
        }
      if(ymin > dOffset && xzplane == 1 && nReverse ==1){
          iGeom_deleteEnt(geom,all[i],&err);
          CHECK("Failed delete entities");
          continue;
        }
      if(xmax < dOffset && yzplane ==1 && nReverse ==0){
          iGeom_deleteEnt(geom,all[i],&err);
          CHECK("Failed delete entities");
          continue;
        }
      if(ymax < dOffset && xzplane == 1 && nReverse ==0){
          iGeom_deleteEnt(geom,all[i],&err);
          CHECK("Failed delete entities");
          continue;
        }
      else{
          if(xzplane ==1 && ymax >dOffset && ymin < dOffset){
              iGeom_sectionEnt(geom, all[i],yzplane,xzplane,0, dOffset, nReverse,&sec,&err);
              CHECK("Failed to section ent");
            }
          if(yzplane ==1 && xmax >dOffset && xmin < dOffset){
              iGeom_sectionEnt(geom, all[i],yzplane,xzplane,0, dOffset,nReverse,&sec,&err);
              CHECK("Failed to section ent");
            }
        }
    }
  return 0;
}

int CNrgen::Rotate_Assm (char &cDir, double &dAngle)
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
  iGeom_getEntities( geom, root_set, iBase_REGION,ARRAY_INOUT(all),&err );
  CHECK("Failed to get all entities");
  // loop and rotate all entities
  for(int i=0; i<all.size(); i++){
      //get the bounding box to decide
      iGeom_rotateEnt(geom,all[i],dAngle,
                      dX, dY, dZ, &err);
      CHECK("Failed rotate entities");
    }
  return 0;
}

int CNrgen::Move_Assm (double &dX,double &dY, double &dZ)
// ---------------------------------------------------------------------------
// Function: move's the model by dX, dY and dZ
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
  SimpleArray<iBase_EntityHandle> all;
  iGeom_getEntities( geom, root_set, iBase_REGION,ARRAY_INOUT(all),&err );
  CHECK("Failed to get all entities");
  // loop and rotate all entities
  for(int i=0; i<all.size(); i++){
      //get the bounding box to decide
      iGeom_moveEnt(geom,all[i],
                    dX, dY, dZ, &err);
      CHECK("Failed move entities");
    }
  return 0;
}

int CNrgen::Create_HexAssm(std::string &szInputString)
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
    return 0;

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

          // now compute the location and create it
          err = ComputePinCentroid(nTempPin, m_Assembly, m, n, dX, dY, dZ);
          ERRORR("Error in function ComputePinCentroid", err);

          // now create the pincell in the location found
          std::cout << "\n--------------------------------------------------"<<std::endl;
          std::cout << " m = " << m <<" n = " << n << std::endl;
          std::cout << "creating pin: " << nTempPin;
          std::cout << " at X Y Z " << dX << " " << dY << " " << dZ << std::endl;

          if(strcmp(m_szInfo.c_str(),"on") == 0)
            m_AssmInfo << nTempPin  << " \t" << m << " \t" << n << " \t" << dX << " \t" << dY << " \t" << dZ << std::endl;

          m_Pincell(nTempPin).GetIntersectFlag(nIFlag);
          if(nIFlag){
              err = CreatePinCell_Intersect(nTempPin, dX, -dY, dZ);
              ERRORR("Error in function CreatePinCell_Intersect", err);
            }
          else{
              err = CreatePinCell(nTempPin, dX, -dY, dZ);
              ERRORR("Error in function CreatePinCell", err);
            }
        }
    }

  // get all the entities (in pins)defined so far, in an entity set - for subtraction later
  iGeom_getEntities( geom, root_set, iBase_REGION, ARRAY_INOUT(in_pins),&err );
  CHECK( "ERROR : getRootSet failed!" );
  std::cout << "Expected pin definitions: " << total_pincells << "\n\nCreating surrounding outer hexes .." << std::endl;

  for (int nTemp = 1; nTemp <= m_nDuct; nTemp ++){
      if(m_nDimensions >0){

          // create outermost hexes
          for(int n=1;n<=m_nDimensions; n++){
              dSide = m_dMAssmPitch(nTemp, n)/(sqrt(3));
              dHeight = m_dMZAssm(nTemp, 2) - m_dMZAssm(nTemp, 1);

              // creating coverings
              iGeom_createPrism(geom, dHeight, 6,
                                dSide, dSide,
                                &assm, &err);
              CHECK("Prism creation failed.");

              // rotate the prism to match the pins
              iGeom_rotateEnt (geom, assm, 30, 0, 0, 1, &err);
              CHECK("Rotation failed failed.");

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
              iGeom_moveEnt(geom, assm, dX,dY,dZ, &err);
              CHECK("Move failed failed.");

              // populate the coverings array
              assms[(nTemp-1)*m_nDimensions + n -1]=assm;
            }
        }
    }
  return 0;
}

int CNrgen::Create_CartAssm(std::string &szInputString)
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
    return 0;

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
              err = ComputePinCentroid(nTempPin, m_Assembly, m, n, dX, dY, dZ);
              ERRORR("Error in function ComputePinCentroid", err);

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
                  err = CreatePinCell_Intersect(nTempPin, dX, -dY, dZ);
                  ERRORR("Error in function CreatePinCell_Intersect", err);
                }
              else{
                  err = CreatePinCell(nTempPin, dX, -dY, dZ);
                  ERRORR("Error in function CreatePinCell", err);
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
  //  iGeom_getEntities( geom, root_set, iBase_REGION, ARRAY_INOUT(in_pins),&err );
  //  CHECK( "ERROR : getRootSet failed!" );


  if(m_nDimensions > 0){

      // create outermost rectangular blocks
      std::cout << "\nCreating surrounding outer blocks .." << std::endl;
      int nCount = -1;
      for(int nTemp = 1; nTemp <= m_nDuct; nTemp ++){
          for(int n=1;n<=m_nDimensions; n++){
              ++nCount;
              dHeight = m_dMZAssm(nTemp, 2) - m_dMZAssm(nTemp, 1);
              iGeom_createBrick(geom, m_dMAssmPitchX(nTemp, n),  m_dMAssmPitchY(nTemp, n), dHeight,
                                &assm, &err);
              CHECK("Prism creation failed.");

              // position the outer block to match the pins
              dX = m_dMAssmPitchX(nTemp, n)/4.0;
              dY =  m_dMAssmPitchY(nTemp, n)/4.0;
              dZ = (m_dMZAssm(nTemp, 2) + m_dMZAssm(nTemp, 1))/2.0;
              std::cout << "Move " <<   dMoveX << " " << dMoveY <<std::endl;
              iGeom_moveEnt(geom, assm, dMoveX,dMoveY,dZ, &err);
              CHECK("Move failed failed.");

              // populate the outer covering array squares
              assms[nCount]=assm;
            }
        }
    }
  return 0;
}

int CNrgen::CreateOuterCovering () 
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
  iGeom_getTagHandle(geom, tag_name, &this_tag, &err, 4);
  CHECK("getTagHandle failed");
  iBase_EntityHandle tmp_vol= NULL, tmp_new= NULL;

  // name the innermost outer covering common for both rectangular and hexagonal assembliees
  if(m_nDimensions >0){
      //    int tag_no = 0;
      for (int nTemp1 = 1; nTemp1 <=m_nDuct; nTemp1++){
          for(int p=1;p<=m_szAssmMatAlias.GetSize();p++){
              if(strcmp ( m_szMMAlias(nTemp1, 1).c_str(), m_szAssmMatAlias(p).c_str()) == 0){
                  sMatName =  m_szAssmMat(p);
                  //	  tag_no=p;
                }
            }

          std::cout << "\ncreated innermost block: " << sMatName << std::endl;

          tmp_vol = assms[(nTemp1 - 1)*m_nDimensions];
          iGeom_setData(geom, tmp_vol, this_tag,
                        sMatName.c_str(), sMatName.size(), &err);
          CHECK("setData failed");

          err = Name_Faces(sMatName, tmp_vol, this_tag);
          ERRORR("Error in function Name_Faces", err);
        }

      int count =0;//index for edge names
      for (int nTemp = 1; nTemp <= m_nDuct; nTemp++){
          //  Naming outermost block edges - sidesets in cubit journal file
          std::cout << "Naming outermost block edges" << std::endl;
          SimpleArray<iBase_EntityHandle> edges;

          iGeom_getEntAdj( geom, assms[nTemp*m_nDimensions-1] , iBase_EDGE,ARRAY_INOUT(edges),
              &err );
          CHECK( "ERROR : getEntAdj failed!" );

          // get the top corner edges of the outer most covering
          std::ostringstream os;
          for (int i = 0; i < edges.size(); ++i){
              iGeom_getEntBoundBox(geom, edges[i],&xmin,&ymin,&zmin,
                                   &xmax,&ymax,&zmax, &err);
              CHECK("getEntBoundBox failed.");
              double dTol = 1e-5; // tolerance for comparing coordinates

              if(fabs(zmax - m_dMZAssm(nTemp, 2)) <  dTol){
                  if(fabs(zmax-zmin) < dTol){

                      //we have a corner edge - name it
                      sMatName="side_edge";
                      ++count;
                      os << sMatName << count;
                      sMatName=os.str();
                      tmp_vol=edges[i];
                      iGeom_setData(geom, tmp_vol, this_tag,
                                    sMatName.c_str(), sMatName.size(), &err);
                      CHECK("setData failed");
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
                  iGeom_copyEnt(geom, assms[(nTemp-1)*m_nDimensions + n-2], &tmp_vol, &err);
                  CHECK("Couldn't copy inner duct wall prism.");

                  // subtract outer most cyl from brick
                  iGeom_subtractEnts(geom, assms[(nTemp-1)*m_nDimensions + n-1], tmp_vol, &tmp_new, &err);
                  CHECK("Subtract of inner from outer failed.");

                  assms[(nTemp-1)*m_nDimensions + n-1]=tmp_new;

                  // name the vols by searching for the full name of the abbreviated Cell Mat
                  for(int p=1;p<=m_szAssmMatAlias.GetSize();p++){
                      if(strcmp ( m_szMMAlias(nTemp, n).c_str(), m_szAssmMatAlias(p).c_str()) == 0){
                          sMatName =  m_szAssmMat(p);
                        }
                    }
                  std::cout << "created: " << sMatName << std::endl;

                  iGeom_setData(geom, tmp_new, this_tag,
                                sMatName.c_str(), sMatName.size(), &err);
                  CHECK("setData failed");
                  err = Name_Faces(sMatName, tmp_new, this_tag);
                }
            }
        }
      std::cout << "\n--------------------------------------------------"<<std::endl;
    }
  return 0;
}

int CNrgen::Subtract_Pins()
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
          for (int i=0; i< (int) cp_inpins[k-1].size();i++){
              iGeom_copyEnt(geom, cp_inpins[k-1][i], &pin_copy[i], &err);
              CHECK("Couldn't copy inner duct wall prism.");
            }

          iBase_EntityHandle tmp_vol = NULL;
          tmp_vol = assms[(k-1)*m_nDimensions];

          // subtract the innermost hex from the pins
          std::cout << "Duct no.: " << k << " subtracting " <<  cp_inpins[k-1].size() << " pins from the duct .. " << std::endl;

#if HAVE_ACIS
          iBase_EntityHandle unite= NULL, tmp_new1;

          // if there are more than one pins
          if( cp_inpins[k-1].size() > 1){

              iGeom_uniteEnts(geom, &cp_inpins[k-1][0], cp_inpins[k-1].size(), &unite, &err);
              CHECK( "uniteEnts failed!" );

              iGeom_subtractEnts(geom, tmp_vol,unite, &tmp_new1, &err);
              CHECK("Couldn't subtract pins from block.");

              tmp_vol = tmp_new1;
              unite = NULL;
              tmp_new1=NULL;
            }
          else{ // only one pin in in_pins
              iGeom_subtractEnts(geom, tmp_vol, cp_inpins[k-1][0], &tmp_new1, &err);
              CHECK("Couldn't subtract pins from block.");
            }
#endif
#if HAVE_OCC
          iBase_EntityHandle tmp_new1 = NULL;
          // if there are more than one pins
          if( cp_inpins[k-1].size() > 1){
              std::cout << "Subtraction is slower in OCC, since each pin is subtracted one by one" << std::endl;
              for (int i=0; i< (int)cp_inpins[k-1].size(); i++){
                  // iGeom_copyEnt(geom, cp_inpins[k-1][i], &unite, &err);
                  iGeom_subtractEnts(geom, tmp_vol,cp_inpins[k-1][i], &tmp_new1, &err);
                  CHECK("Couldn't subtract pins from block.");
                  tmp_vol = tmp_new1;
                  tmp_new1=NULL;
                }

            }
          else{ // only one pin in in_pins
              iGeom_subtractEnts(geom, tmp_vol, cp_inpins[k-1][0], &tmp_new1, &err);
              CHECK("Couldn't subtract pins from block.");
            }
#endif

        }
      std::cout << "\n--------------------------------------------------"<<std::endl;
    }
  else{
      std::cout <<"Nothing to subtract" << std::endl;
    }
  return 0;
}

int CNrgen::Imprint_Merge()
// ---------------------------------------------------------------------------
// Function: Imprint and Merge 
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
   // getting all entities for merge and imprint
   SimpleArray<iBase_EntityHandle> entities;
   iGeom_getEntities( geom, root_set, iBase_REGION, ARRAY_INOUT(entities),&err );
   CHECK( "ERROR : getRootSet failed!" );

    //  now imprint
    std::cout << "\n\nImprinting...." << std::endl;
    iGeom_imprintEnts(geom, ARRAY_IN(entities),&err);
    CHECK("Imprint failed.");
    std::cout << "\n--------------------------------------------------"<<std::endl;

     // merge tolerance
     double dTol = 1e-4;
     // now  merge
     std::cout << "\n\nMerging...." << std::endl;
     iGeom_mergeEnts(geom, ARRAY_IN(entities), dTol, &err);
     CHECK("Merge failed.");
     std::cout <<"merging finished."<< std::endl;
     std::cout << "\n--------------------------------------------------"<<std::endl;
  return 0;
}

int CNrgen::Create2DSurf ()
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
  iGeom_getEntities( geom, root_set, iBase_REGION,ARRAY_INOUT(all_geom),&err );
  CHECK( "ERROR : Failed to get all geom" );

  // get all the surfaces in the model
  iGeom_getArrAdj( geom, ARRAY_IN(all_geom) , iBase_FACE, ARRAY_INOUT(surfs),
                   &offset, &offset_alloc, &offset_size, &err );
  CHECK( "ERROR : getArrAdj failed!" );

  SimpleArray<double> max_corn, min_corn;
  iGeom_getArrBoundBox( geom, ARRAY_IN(surfs), iBase_INTERLEAVED,
                        ARRAY_INOUT( min_corn ),
                        ARRAY_INOUT( max_corn ),
                        &err );
  CHECK( "Problems getting max surf for rotation." );

  // find the number of surfaces 't' for array allocation
  int nTemp = 1;
  double dTol = 1e-5;
  double dtop = m_dMZAssm(nTemp, 2);
  for (int i = 0; i < surfs.size(); ++i){
      if((abs(max_corn[3*i+2] -  dtop) < dTol) && (abs(min_corn[3*i+2] - dtop)<dTol))
        t++;
    }

  // allocate arrays
  SimpleArray<iBase_EntityHandle> max_surfs(t);
  SimpleArray<iBase_EntityHandle> new_surfs(t);
  t=0;

  // store the max surfaces in max_surfs
  for (int i = 0; i < surfs.size(); ++i){

      // locate surfaces for which max and min zcoord is same as maxz coord
      if((abs(max_corn[3*i+2] -  dtop) < dTol) && (abs(min_corn[3*i+2] - dtop) < dTol)){
          max_surfs[t] = surfs[i];
          t++;
        }
    }

  // make a copy of max_surfs
  for(int i = 0; i < max_surfs.size(); ++i){
      iGeom_copyEnt(geom, max_surfs[i], &new_surfs[i], &err);
      CHECK( "Problems creating surface." );
    }

  // delete all the old ents
  for(int i=0; i<all_geom.size(); i++){
      iGeom_deleteEnt(geom, all_geom[i], &err);
      CHECK( "Problems deleting cyls." );
    }
  // position the final assembly at the center
  // get the assembly on z=0 plane
  double zcenter = m_dMZAssm(nTemp, 2)/2.0;//move up
  SimpleArray<iBase_EntityHandle> all;
  iGeom_getEntities( geom, root_set, iBase_REGION,ARRAY_INOUT(all),&err );
  CHECK("Failed to get all entities");

  for(int i=0; i<all.size(); i++){
      iGeom_moveEnt(geom,all[i],0,0,-zcenter,&err);
      CHECK("Failed to move entities");
    }
  std::cout << "--------------------------------------------------"<<std::endl;

  free(offset);
  return 0;
}


int CNrgen::CreatePinCell(int i, double dX, double dY, double dZ)
//---------------------------------------------------------------------------
//Function: Create pincell i in location dX dY and dZ
//Input:    none
//Output:   none
//---------------------------------------------------------------------------
{
  int nRadii=0, nCyl=0, nCells = 0;
  double dCylMoveX = 0.0, dCylMoveY = 0.0, dHeightTotal = 0.0;
  double dHeight =0.0,dZMove = 0.0, PX = 0.0,PY = 0.0,PZ = 0.0, dP=0.0;
  CVector<double> dVCylZPos(2), dVCylXYPos(2), dVStartZ, dVEndZ;;
  CVector<std::string> szVMatName, szVMatAlias, szVCellMat;
  iBase_EntityHandle cell = NULL, cyl= NULL, tmp_vol= NULL,tmp_vol1= NULL, tmp_new= NULL;
  std::vector<iBase_EntityHandle> cp_in;
  // name tag handle
  iBase_TagHandle this_tag= NULL;
  char* tag_name = (char*)"NAME";

  std::string sMatName = "";
  std::string sMatName1 = "";
  int nDuctIndex = -1;

  if(strcmp(m_szInfo.c_str(),"on") == 0){
      std::ostringstream os;
      pin_name = "_xp";
      os << (m_nTotalPincells + m_nStartpinid - 1);
      os << "_";
      std::string pid = os.str(); //retrieve as a string
      pin_name+=pid;
    }

  // get tag handle for 'NAME' tag, already created as iGeom instance is created
  iGeom_getTagHandle(geom, tag_name, &this_tag, &err, 4);
  CHECK("getTagHandle failed");

  // get cell material
  m_Pincell(i).GetCellMatSize(nCells);
  SimpleArray<iBase_EntityHandle> cells(nCells);

  // branch when cells are present
  if(nCells > 0){
      dVStartZ.SetSize(nCells);
      dVEndZ.SetSize(nCells);
      szVCellMat.SetSize(nCells);
      m_Pincell(i).GetCellMat(dVStartZ, dVEndZ, szVCellMat);

      // get cylinder data
      m_Pincell(i).GetNumCyl(nCyl);

      for(int n=1;n<=nCells; n++){
          // get cylinder locations
          m_Pincell(i).GetCylZPos(n, dVCylZPos);
          nDuctIndex = -1;
          dHeight = fabs(dVEndZ(n) - dVStartZ(n));
          // get the index for cp_inpins based on Z-heights
          for (int dd = 1; dd <= m_nDuct; dd++){
              if((m_dMZAssm(dd, 2)) >= (dVCylZPos(2)) && (m_dMZAssm(dd, 1)) >= (dVCylZPos(1)))
                nDuctIndex = dd;
              if (nDuctIndex != -1)
                break;
            }
          if(m_szGeomType =="hexagonal"){

              m_Pincell(i).GetPitch(dP, dHeightTotal); // this dHeight is not used in creation
              double dSide = dP/(sqrt(3));

              if(nCells >0){
                  // create prism
                  iGeom_createPrism(geom, dHeight, 6,
                                    dSide, dSide,
                                    &cell, &err);
                  CHECK("Prism creation failed.");
                }
            }
          // if rectangular geometry
          if(m_szGeomType =="rectangular"){

              m_Pincell(i).GetPitch(PX, PY, PZ);

              if(nCells >0){
                  // create brick
                  iGeom_createBrick( geom,PX,PY,dHeight,&cell,&err );
                  CHECK("Couldn't create pincell.");
                }
            }

          dZMove = (dVStartZ(n)+dVEndZ(n))/2.0;
          if(nCells > 0){
              // position the brick in assembly
              iGeom_moveEnt(geom, cell, dX, dY, dZMove, &err);
              CHECK("Couldn't move cell.");
              cells[n-1]=cell;

              //search for the full name of the abbreviated Cell Mat and set name
              for(int p=1;p<= m_szAssmMatAlias.GetSize();p++){
                  if(strcmp (szVCellMat(n).c_str(), m_szAssmMatAlias(p).c_str()) == 0){
                      sMatName = m_szAssmMat(p);
                    }
                }
              std::cout << "created: " << sMatName << std::endl;
              iGeom_setData(geom, cell, this_tag,
                            sMatName.c_str(), sMatName.size(), &err);
              CHECK("setData failed");

              if(strcmp(m_szInfo.c_str(),"on") == 0){
                  iGeom_setData(geom, cell, this_tag,
                                pin_name.c_str(), pin_name.size(), &err);
                  std::cout << "Naming pin body :" <<  pin_name << std::endl;
                }


              Name_Faces(sMatName, cell, this_tag);
              CHECK("Name_Faces failed");
            }
          // loop and create cylinders
          if(nCyl > 0){
              m_Pincell(i).GetCylSizes(n, nRadii);
              SimpleArray<iBase_EntityHandle> cyls(nRadii);

              //declare variables
              CVector<double> dVCylRadii(2*nRadii);
              CVector<std::string> szVMat(nRadii);
              CVector<std::string> szVCylMat(nRadii);
              int nType = 0;
              //get values
              m_Pincell(i).GetCylRadii(n, dVCylRadii);
              m_Pincell(i).GetCylPos(n, dVCylXYPos);
              m_Pincell(i).GetCylMat(n, szVCylMat);
              m_Pincell(i).GetCylZPos(n, dVCylZPos);
              m_Pincell(i).GetCellType(n, nType);

              dHeight = dVCylZPos(2)-dVCylZPos(1);

              for (int m=1; m<=nRadii; m++){

                  if (nType == 0){
                      iGeom_createCylinder(geom, dHeight, dVCylRadii(m), dVCylRadii(m),
                                           &cyl, &err);
                      CHECK("Couldn't create fuel rod.");
                      std::cout << m << ": Creating cylinder with radii " << dVCylRadii(m) << std::endl;
                    }
                  else{
                      iGeom_createCone(geom, dHeight, dVCylRadii(2*m-1), dVCylRadii(2*m-1), dVCylRadii(2*m),
                                       &cyl, &err);
                      CHECK("Couldn't create fuel rod.");
                    }
                  // move their centers and also move to the assembly location  ! Modify if cyl is outside brick
                  dCylMoveX = dVCylXYPos(1)+dX;
                  dCylMoveY = dVCylXYPos(2)+dY;
                  dZMove = (dVCylZPos(1)+dVCylZPos(2))/2.0;

                  iGeom_moveEnt(geom, cyl, dCylMoveX,dCylMoveY,dZMove, &err);
                  CHECK("Couldn't move cyl.");
                  cyls[m-1] = cyl;
                }

              if(nCells > 0){
                  // copy cyl before subtract
                  iGeom_copyEnt(geom, cyls[nRadii-1], &tmp_vol, &err);
                  CHECK("Couldn't copy inner duct wall prism.");

                  // subtract outer most cyl from brick
                  iGeom_subtractEnts(geom, cells[n-1], tmp_vol, &tmp_new, &err);
                  CHECK("Subtract of inner from outer failed.");

                  // copy the new into the cyl array
                  cells[n-1] = tmp_new; cell = tmp_new;

                }
              cp_in.push_back(tmp_new);

              //set tag on inner most cylinder, search for the full name of the abbreviated Cell Mat
              for(int p=1;p<=m_szAssmMatAlias.GetSize();p++){
                  if(strcmp (szVCylMat(1).c_str(), m_szAssmMatAlias(p).c_str()) == 0){
                      sMatName = m_szAssmMat(p);
                    }
                }
              tmp_vol1=cyls[0]; //inner most cyl

              cp_in.push_back(tmp_vol1);
              iGeom_setData(geom, tmp_vol1, this_tag,
                            sMatName.c_str(), 10, &err);
              CHECK("setData failed");
              if(strcmp(m_szInfo.c_str(),"on") == 0){
                  iGeom_setData(geom, tmp_vol1, this_tag,
                                pin_name.c_str(), pin_name.size(), &err);
                  std::cout << "Naming pin body :" <<  pin_name<< std::endl;
                }

              err= Name_Faces(sMatName, tmp_vol1, this_tag);
              ERRORR("Error in function Name_Faces", err);

              // other cyl annulus after substraction
              for (int b=nRadii; b>1; b--){

                  iGeom_copyEnt(geom, cyls[b-2], &tmp_vol, &err);
                  CHECK("Couldn't copy inner duct wall prism.");

                  //subtract tmp vol from the outer most
                  iGeom_subtractEnts(geom, cyls[b-1], tmp_vol, &tmp_new, &err);
                  CHECK("Subtract of inner from outer failed.");

                  // now search for the full name of the abbreviated Cell Mat
                  //	  int tag_no;
                  for(int p=1;p<=m_szAssmMatAlias.GetSize();p++){
                      if(strcmp (szVCylMat(b).c_str(), m_szAssmMatAlias(p).c_str()) == 0){
                          //	      tag_no = p;
                          sMatName =  m_szAssmMat(p);
                        }
                    }
                  std::cout << "created: " << sMatName << std::endl;
                  cp_in.push_back(tmp_new);
                  // set the name of the annulus
                  iGeom_setData(geom, tmp_new, this_tag,
                                sMatName.c_str(),sMatName.size(), &err);
                  CHECK("setData failed");

                  if(strcmp(m_szInfo.c_str(),"on") == 0){
                      iGeom_setData(geom, tmp_new, this_tag,
                                    pin_name.c_str(), pin_name.size(), &err);
                      std::cout << "Naming pin body :" <<  pin_name<< std::endl;
                    }
                  err= Name_Faces(sMatName, tmp_new, this_tag);
                  ERRORR("Error in function Name_Faces", err);

                  // copy the new into the cyl array
                  cyls[b-1] = tmp_new;
                  tmp_vol=NULL;
                }
            }
          if(nDuctIndex > 0){
              for (int count = 0; count < (int) cp_in.size(); count++)
                cp_inpins[nDuctIndex-1].push_back(cp_in[count]);
            }
          cp_in.clear();
        }
    }
  // this branch of the routine is responsible for creating cylinders with '0' cells
  if(nCells == 0){

      // get cylinder data
      m_Pincell(i).GetNumCyl(nCyl);
      nCells = nCyl;

      for(int n=1;n<=nCells; n++){
          nDuctIndex = -1;
          if(m_szGeomType =="hexagonal"){

              m_Pincell(i).GetPitch(dP, dHeightTotal); // this dHeight is not used in creation
            }
          // if rectangular geometry
          if(m_szGeomType =="rectangular"){

              m_Pincell(i).GetPitch(PX, PY, PZ);
            }

          // loop and create cylinders
          if(nCyl > 0){
              m_Pincell(i).GetCylSizes(n, nRadii);
              SimpleArray<iBase_EntityHandle> cyls(nRadii);

              //declare variables
              CVector<double> dVCylRadii(2*nRadii);
              CVector<std::string> szVMat(nRadii);
              CVector<std::string> szVCylMat(nRadii);
              int nType = 0;
              //get values
              m_Pincell(i).GetCylRadii(n, dVCylRadii);
              m_Pincell(i).GetCylPos(n, dVCylXYPos);
              m_Pincell(i).GetCylMat(n, szVCylMat);
              m_Pincell(i).GetCylZPos(n, dVCylZPos);
              m_Pincell(i).GetCellType(n, nType);

              dHeight = dVCylZPos(2)-dVCylZPos(1);

              // get the index for cp_inpins based on Z-heights
              for (int dd = 1; dd <= m_nDuct; dd++){
                  if((m_dMZAssm(dd, 2)) >= (dVCylZPos(2)) && (m_dMZAssm(dd, 1)) >= (dVCylZPos(1)))
                    nDuctIndex = dd;
                  if (nDuctIndex != -1)
                    break;
                }

              for (int m=1; m<=nRadii; m++){
                  if (nType == 0){
                      iGeom_createCylinder(geom, dHeight, dVCylRadii(m), dVCylRadii(m),
                                           &cyl, &err);
                      CHECK("Couldn't create fuel rod.");
                    }
                  else{
                      iGeom_createCone(geom, dHeight, dVCylRadii(2*m - 1), dVCylRadii(2*m - 1), dVCylRadii(2*m),
                                       &cyl, &err);
                      CHECK("Couldn't create fuel rod.");
                    }

                  // move their centers and also move to the assembly location  ! Modify if cyl is outside brick
                  dCylMoveX = dVCylXYPos(1)+dX;
                  dCylMoveY = dVCylXYPos(2)+dY;
                  dZMove = (dVCylZPos(1)+dVCylZPos(2))/2.0;

                  iGeom_moveEnt(geom, cyl, dCylMoveX,dCylMoveY,dZMove, &err);
                  CHECK("Couldn't move cyl.");
                  cyls[m-1] = cyl;
                }

              //set tag on inner most cylinder, search for the full name of the abbreviated Cell Mat
              for(int p=1;p<=m_szAssmMatAlias.GetSize();p++){
                  if(strcmp (szVCylMat(1).c_str(), m_szAssmMatAlias(p).c_str()) == 0){
                      sMatName = m_szAssmMat(p);
                    }
                }
              std::cout << "created: " << sMatName << std::endl;
              tmp_vol1=cyls[0]; //inner most cyl

              cp_in.push_back(tmp_vol1);

              iGeom_setData(geom, tmp_vol1, this_tag,
                            sMatName.c_str(), 10, &err);
              CHECK("setData failed");

              if(strcmp(m_szInfo.c_str(),"on") == 0){
                  iGeom_setData(geom, tmp_vol1, this_tag,
                                pin_name.c_str(), pin_name.size(), &err);
                  std::cout << "Naming pin body :" <<  pin_name<< std::endl;
                }

              err= Name_Faces(sMatName, tmp_vol1, this_tag);
              ERRORR("Error in function Name_Faces", err);

              // other cyl annulus after substraction
              for (int b=nRadii; b>1; b--){

                  iGeom_copyEnt(geom, cyls[b-2], &tmp_vol, &err);
                  CHECK("Couldn't copy inner duct wall prism.");

                  //subtract tmp vol from the outer most
                  iGeom_subtractEnts(geom, cyls[b-1], tmp_vol, &tmp_new, &err);
                  CHECK("Subtract of inner from outer failed.");

                  // now search for the full name of the abbreviated Cell Mat
                  for(int p=1;p<=m_szAssmMatAlias.GetSize();p++){
                      if(strcmp (szVCylMat(b).c_str(), m_szAssmMatAlias(p).c_str()) == 0){
                          sMatName =  m_szAssmMat(p);
                        }
                    }
                  std::cout <<"created: " << sMatName << std::endl;

                  cp_in.push_back(tmp_new);

                  // set the name of the annulus
                  iGeom_setData(geom, tmp_new, this_tag,
                                sMatName.c_str(),sMatName.size(), &err);
                  CHECK("setData failed");

                  if(strcmp(m_szInfo.c_str(),"on") == 0){
                      iGeom_setData(geom, tmp_new, this_tag,
                                    pin_name.c_str(), pin_name.size(), &err);
                      std::cout << "Naming pin body :" <<  pin_name<< std::endl;
                    }

                  err= Name_Faces(sMatName, tmp_new, this_tag);
                  ERRORR("Error in function Name_Faces", err);

                  // copy the new into the cyl array
                  cyls[b-1] = tmp_new;
                }
            }
          if(nDuctIndex > 0){
              for (int count = 0; count < (int) cp_in.size(); count++)
                cp_inpins[nDuctIndex-1].push_back(cp_in[count]);
            }
          cp_in.clear();
        }
    }
  return 0;
}


int CNrgen::CreatePinCell_Intersect(int i, double dX, double dY, double dZ)
//---------------------------------------------------------------------------
//Function: Create pincell i in location dX dY and dZ
//Input:    none
//Output:   none
//---------------------------------------------------------------------------
{
  int nRadii=0, nCyl=0, nCells = 0;
  double dCylMoveX = 0.0, dCylMoveY = 0.0, dHeightTotal = 0.0;
  double dHeight =0.0,dZMove = 0.0, PX = 0.0,PY = 0.0,PZ = 0.0, dP=0.0;
  CVector<double> dVCylZPos(2), dVCylXYPos(2), dVEndZ, dVStartZ;
  CVector<std::string> szVMatName, szVMatAlias, szVCellMat;
  iBase_EntityHandle cell = NULL, cyl= NULL, tmp_vol1= NULL, tmp_new= NULL, cell_copy = NULL, intersec = NULL;
  std::vector<iBase_EntityHandle> cp_in;

  // name tag handle
  iBase_TagHandle this_tag= NULL;
  char* tag_name = (char*)"NAME";

  std::string sMatName = "";
  std::string sMatName0 = "";
  std::string sMatName1 = "";
  int nDuctIndex = -1;

  if(strcmp(m_szInfo.c_str(),"on") == 0){
      std::ostringstream os;
      pin_name = "_xp";
      os << (m_nTotalPincells + m_nStartpinid - 1);
      os << "_";
      std::string pid = os.str(); //retrieve as a string
      pin_name+=pid;
    }

  // get tag handle for 'NAME' tag, already created as iGeom instance is created
  iGeom_getTagHandle(geom, tag_name, &this_tag, &err, 4);
  CHECK("getTagHandle failed");

  // get cell material
  m_Pincell(i).GetCellMatSize(nCells);
  m_Pincell(i).GetNumCyl(nCyl);
  SimpleArray<iBase_EntityHandle> cells(nCells);

  // branch when cells are present
  if(nCells > 0){
      dVStartZ.SetSize(nCells);
      dVEndZ.SetSize(nCells);
      szVCellMat.SetSize(nCells);
      m_Pincell(i).GetCellMat(dVStartZ, dVEndZ, szVCellMat);

      // get cylinder data
      m_Pincell(i).GetNumCyl(nCyl);

      for(int n=1;n<=nCells; n++){

          dHeight = fabs(dVEndZ(n) - dVStartZ(n));
          // get the index for cp_inpins based on Z-heights
          for (int dd = 1; dd <= m_nDuct; dd++){
              if((m_dMZAssm(dd, 2)) >= (dVCylZPos(2)) && (m_dMZAssm(dd, 1)) >= (dVCylZPos(1)))
                nDuctIndex = dd;
              if (nDuctIndex != -1)
                break;
            }
          if(m_szGeomType =="hexagonal"){

              m_Pincell(i).GetPitch(dP, dHeightTotal); // this dHeight is not used in creation

              double dSide = dP/(sqrt(3));

              if(nCells >0){
                  // create prism
                  iGeom_createPrism(geom, dHeight, 6,
                                    dSide, dSide,
                                    &cell, &err);
                  CHECK("Prism creation failed.");
                }
            }
          // if rectangular geometry
          if(m_szGeomType =="rectangular"){

              m_Pincell(i).GetPitch(PX, PY, PZ);

              if(nCells >0){
                  // create brick
                  iGeom_createBrick( geom,PX,PY,dHeight,&cell,&err );
                  CHECK("Couldn't create pincell.");
                }
            }

          dZMove = (dVStartZ(n)+dVEndZ(n))/2.0;

          if(nCells > 0){
              // position the brick in assembly
              iGeom_moveEnt(geom, cell, dX, dY, dZMove, &err);
              CHECK("Couldn't move cell.");
              cells[n-1]=cell;
            }
          // loop and create cylinders
          if(nCyl > 0){
              m_Pincell(i).GetCylSizes(n, nRadii);
              SimpleArray<iBase_EntityHandle> cyls(nRadii);
              SimpleArray<iBase_EntityHandle> cell_copys(nRadii);
              SimpleArray<iBase_EntityHandle> intersec_main(nRadii);
              iBase_EntityHandle  tmp_intersec;
              //declare variables
              CVector<double> dVCylRadii(2*nRadii);
              CVector<std::string> szVMat(nRadii);
              CVector<std::string> szVCylMat(nRadii);
              int nType = 0;
              //get values
              m_Pincell(i).GetCylRadii(n, dVCylRadii);
              m_Pincell(i).GetCylPos(n, dVCylXYPos);
              m_Pincell(i).GetCylMat(n, szVCylMat);
              m_Pincell(i).GetCylZPos(n, dVCylZPos);
              m_Pincell(i).GetCellType(n, nType);

              dHeight = dVCylZPos(2)-dVCylZPos(1);

              for (int m=1; m<=nRadii; m++){
                  if (nType == 0){
                      iGeom_createCylinder(geom, dHeight, dVCylRadii(m), dVCylRadii(m),
                                           &cyl, &err);
                      CHECK("Couldn't create fuel rod.");
                    }
                  else{
                      iGeom_createCone(geom, dHeight, dVCylRadii(2*m-1), dVCylRadii(2*m-1), dVCylRadii(2*m),
                                       &cyl, &err);
                      CHECK("Couldn't create fuel rod.");
                    }

                  // move their centers and also move to the assembly location  ! Modify if cyl is outside brick
                  dCylMoveX = dVCylXYPos(1)+dX;
                  dCylMoveY = dVCylXYPos(2)+dY;
                  dZMove = (dVCylZPos(1)+dVCylZPos(2))/2.0;

                  iGeom_moveEnt(geom, cyl, dCylMoveX,dCylMoveY,dZMove, &err);
                  CHECK("Couldn't move cyl.");
                  cyls[m-1] = cyl;


                  //copy cell nRadii  times for intersection with cylinders
                  iGeom_copyEnt(geom, cells[n-1], &cell_copy, &err);
                  CHECK("Couldn't copy inner duct wall prism.");
                  cell_copys[m-1] = cell_copy;

                  iGeom_intersectEnts(geom, cell_copys[m-1], cyls[m-1],&intersec,&err);
                  CHECK("intersection failed");
                  intersec_main[m-1] = intersec;
                  intersec = NULL;
                }

              //set tag on inner most cylinder, search for the full name of the abbreviated Cell Mat
              tmp_vol1=intersec_main[0];
              for(int p=1;p<=m_szAssmMatAlias.GetSize();p++){
                  if(strcmp (szVCylMat(1).c_str(), m_szAssmMatAlias(p).c_str()) == 0){
                      sMatName = m_szAssmMat(p);
                    }
                }

              cp_in.push_back(tmp_vol1);

              iGeom_setData(geom, tmp_vol1, this_tag,
                            sMatName.c_str(), 10, &err);
              CHECK("setData failed");

              if(strcmp(m_szInfo.c_str(),"on") == 0){
                  iGeom_setData(geom, tmp_vol1, this_tag,
                                pin_name.c_str(), pin_name.size(), &err);
                  std::cout << "Naming pin body :" <<  pin_name<< std::endl;
                }
              err= Name_Faces(sMatName, tmp_vol1, this_tag);
              ERRORR("Error in function Name_Faces", err);

              // copy the outermost cyl
              iGeom_copyEnt(geom, intersec_main[nRadii-1], &tmp_intersec, &err);
              CHECK("Couldn't copy inner duct wall prism.");

              // subtract the outermost cyl from the cell
              iGeom_subtractEnts(geom, cells[n-1], tmp_intersec, &tmp_new, &err);
              CHECK("Subtract of inner from outer failed.");
              // now search for the full name of the abbreviated Cell Mat
              //	int tag_no;
              for(int p=1;p<=m_szAssmMatAlias.GetSize();p++){
                  if(strcmp (szVCellMat(n).c_str(), m_szAssmMatAlias(p).c_str()) == 0){
                      //	    tag_no = p;
                      sMatName =  m_szAssmMat(p);
                    }
                }
              std::cout << "created: " << sMatName << std::endl;

              cp_in.push_back(tmp_new);

              // set the name of the annulus
              iGeom_setData(geom, tmp_new, this_tag,
                            sMatName.c_str(),sMatName.size(), &err);
              CHECK("setData failed");

              iGeom_setData(geom, tmp_new, this_tag,
                            pin_name.c_str(), pin_name.size(), &err);
              std::cout << "Naming pin body :" <<  pin_name<< std::endl;

              err = Name_Faces(sMatName, tmp_new, this_tag);
              ERRORR("Error in function Name_Faces", err);

              for (int b=nRadii; b>1; b--){

                  iGeom_copyEnt(geom, intersec_main[b-2], &tmp_intersec, &err);
                  CHECK("Couldn't copy inner duct wall prism.");
                  //subtract tmp vol from the outer most
                  iGeom_subtractEnts(geom, intersec_main[b-1], tmp_intersec, &tmp_new, &err);
                  CHECK("Subtract of inner from outer failed.");

                  // now search for the full name of the abbreviated Cell Mat
                  //	  int tag_no;
                  for(int p=1;p<=m_szAssmMatAlias.GetSize();p++){
                      if(strcmp (szVCylMat(b).c_str(), m_szAssmMatAlias(p).c_str()) == 0){
                          //	      tag_no = p;
                          sMatName =  m_szAssmMat(p);
                        }
                    }
                  std::cout << "created: " << sMatName << std::endl;

                  cp_in.push_back(tmp_new);

                  // set the name of the annulus
                  iGeom_setData(geom, tmp_new, this_tag,
                                sMatName.c_str(),sMatName.size(), &err);
                  CHECK("setData failed");

                  if(strcmp(m_szInfo.c_str(),"on") == 0){
                      iGeom_setData(geom, tmp_new, this_tag,
                                    pin_name.c_str(), pin_name.size(), &err);
                      std::cout << "Naming pin body :" <<  pin_name<< std::endl;
                    }
                  err = Name_Faces(sMatName, tmp_new, this_tag);
                  ERRORR("Error in function Name_Faces", err);

                  // copy the new into the cyl array
                  cyls[b-1] = tmp_new;

                }
            }
          if(nDuctIndex > 0){
              for (int count = 0; count < (int) cp_in.size(); count++)
                cp_inpins[nDuctIndex-1].push_back(cp_in[count]);
            }
          cp_in.clear();
        }
    }
  // this branch of the routine is responsible for creating cylinders with '0' cells
  if(nCells == 0){

      // get cylinder data
      m_Pincell(i).GetNumCyl(nCyl);
      nCells = nCyl;
      cells.resize(nCells);

      for(int n=1;n<=nCells; n++){

          // get some cylinder parameters to create the cell material for intersection
          m_Pincell(i).GetCylZPos(n, dVCylZPos);
          dHeight = dVCylZPos(2)-dVCylZPos(1);
          dZMove = (dVCylZPos(1)+dVCylZPos(2))/2.0;

          if(m_szGeomType =="hexagonal"){

              m_Pincell(i).GetPitch(dP, dHeightTotal); // this dHeight is not used in creation
              double dSide = dP/(sqrt(3));

              iGeom_createPrism(geom, dHeight, 6,
                                dSide, dSide,
                                &cell, &err);
              CHECK("Prism creation failed.");

            }
          // if rectangular geometry
          if(m_szGeomType =="rectangular"){

              m_Pincell(i).GetPitch(PX, PY, PZ);
              // create brick
              iGeom_createBrick( geom,PX,PY,PZ, &cell,&err );
              CHECK("Couldn't create pincell.");
            }

          iGeom_moveEnt(geom, cell, dX, dY, dZMove, &err);
          CHECK("Couldn't move cell.");

          cells[n-1]=cell;
          // loop and create cylinders
          if(nCyl > 0){
              m_Pincell(i).GetCylSizes(n, nRadii);

              //declare variables
              SimpleArray<iBase_EntityHandle> cyls(nRadii), cell_copys(nRadii), intersec_main(nRadii), intersec_copy(nRadii);
              iBase_EntityHandle  tmp_intersec;
              CVector<double> dVCylRadii(2*nRadii);
              CVector<std::string> szVMat(nRadii), szVCylMat(nRadii);
              int nType = 0;
              //get values
              m_Pincell(i).GetCylRadii(n, dVCylRadii);
              m_Pincell(i).GetCylPos(n, dVCylXYPos);
              m_Pincell(i).GetCylMat(n, szVCylMat);
              m_Pincell(i).GetCellType(n, nType);


              // get the index for cp_inpins based on Z-heights
              for (int dd = 1; dd <= m_nDuct; dd++){
                  if((m_dMZAssm(dd, 2)) >= (dVCylZPos(2)) && (m_dMZAssm(dd, 1)) >= (dVCylZPos(1)))
                    nDuctIndex = dd;
                  if (nDuctIndex != -1)
                    break;
                }

              for (int m=1; m<=nRadii; m++){
                  if (nType == 0){
                      iGeom_createCylinder(geom, dHeight, dVCylRadii(m), dVCylRadii(m),
                                           &cyl, &err);
                      CHECK("Couldn't create fuel rod.");
                    }
                  else{
                      iGeom_createCone(geom, dHeight, dVCylRadii(2*m-1), dVCylRadii(2*m-1), dVCylRadii(2*m),
                                       &cyl, &err);
                      CHECK("Couldn't create fuel rod.");
                    }

                  // move their centers and also move to the assembly location  ! Modify if cyl is outside brick
                  dCylMoveX = dVCylXYPos(1)+dX;
                  dCylMoveY = dVCylXYPos(2)+dY;

                  iGeom_moveEnt(geom, cyl, dCylMoveX,dCylMoveY,dZMove, &err);
                  CHECK("Couldn't move cyl.");
                  cyls[m-1] = cyl;

                  //copy cell nRadii  times for intersection with cylinders
                  iGeom_copyEnt(geom, cells[n-1], &cell_copy, &err);
                  CHECK("Couldn't copy inner duct wall prism.");
                  //	  cell_copys[m-1] = cell_copy;

                  iGeom_intersectEnts(geom, cell_copy, cyls[m-1],&intersec,&err);
                  CHECK("intersection failed");
                  intersec_main[m-1] = intersec;
                  intersec = NULL;
                }

              //set tag on inner most cylinder, search for the full name of the abbreviated Cell Mat
              tmp_vol1=intersec_main[0];
              for(int p=1;p<=m_szAssmMatAlias.GetSize();p++){
                  if(strcmp (szVCylMat(1).c_str(), m_szAssmMatAlias(p).c_str()) == 0){
                      sMatName = m_szAssmMat(p);
                    }
                }

              cp_in.push_back(tmp_vol1);

              iGeom_setData(geom, tmp_vol1, this_tag,
                            sMatName.c_str(), 10, &err);
              CHECK("setData failed");

              if(strcmp(m_szInfo.c_str(),"on") == 0){
                  iGeom_setData(geom, tmp_vol1, this_tag,
                                pin_name.c_str(), pin_name.size(), &err);
                  std::cout << "Naming pin body :" <<  pin_name<< std::endl;
                }

              err= Name_Faces(sMatName, tmp_vol1, this_tag);
              ERRORR("Error in function Name_Faces", err);

              // delete the cell as this is the case when no. cell material is specified
              iGeom_deleteEnt(geom, cells[n-1], &err);
              CHECK("Entity deletion failed");


              // other cyl annulus after substraction
              for (int b=nRadii; b>1; b--){

                  iGeom_copyEnt(geom, intersec_main[b-2], &tmp_intersec, &err);
                  CHECK("Couldn't copy inner duct wall prism.");
                  //subtract tmp vol from the outer most
                  iGeom_subtractEnts(geom, intersec_main[b-1], tmp_intersec, &tmp_new, &err);
                  CHECK("Subtract of inner from outer failed.");

                  // now search for the full name of the abbreviated Cell Mat
                  //	  int tag_no;
                  for(int p=1;p<=m_szAssmMatAlias.GetSize();p++){
                      if(strcmp (szVCylMat(b).c_str(), m_szAssmMatAlias(p).c_str()) == 0){
                          //	      tag_no = p;
                          sMatName =  m_szAssmMat(p);
                        }
                    }
                  std::cout << "created: " << sMatName << std::endl;

                  cp_in.push_back(tmp_new);

                  // set the name of the annulus
                  iGeom_setData(geom, tmp_new, this_tag,
                                sMatName.c_str(),sMatName.size(), &err);
                  CHECK("setData failed");
                  if(strcmp(m_szInfo.c_str(),"on") == 0){
                      iGeom_setData(geom, tmp_new, this_tag,
                                    pin_name.c_str(), pin_name.size(), &err);
                      std::cout << "Naming pin body :" <<  pin_name<< std::endl;
                    }

                  err = Name_Faces(sMatName, tmp_new, this_tag);
                  ERRORR("Error in function Name_Faces", err);

                  // copy the new into the cyl array
                  cyls[b-1] = tmp_new;

                }
            }
          if(nDuctIndex > 0){
              for (int count = 0; count < (int) cp_in.size(); count++)
                cp_inpins[nDuctIndex-1].push_back(cp_in[count]);
            }
          cp_in.clear();
        }
    }
  return 0;
}
