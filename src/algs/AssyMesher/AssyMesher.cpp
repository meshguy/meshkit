#include <stdio.h>
#include "meshkit/AssyMesher.hpp"
#include "meshkit/MKCore.hpp"
#include "meshkit/SizingFunction.hpp"
#include "meshkit/RegisterMeshOp.hpp"
#include "meshkit/LocalSet.hpp"
#include "meshkit/ModelEnt.hpp"

#include "iMesh_extensions.h"
#include "MBCN.h"

#include <set>


namespace MeshKit
{
// static registration of this  mesh scheme
moab::EntityType AssyMesher_tps[] = { moab::MBTET,
                                      moab::MBQUAD,
                                      moab::MBTRI,
                                      moab::MBHEX,
                                      moab::MBMAXTYPE};
const moab::EntityType* AssyMesher::output_types()
{ return AssyMesher_tps; }

AssyMesher::AssyMesher(MKCore *mk, const MEntVector &me_vec)
: MeshScheme( mk, me_vec),
  igeom(mk->igeom_instance()), imesh(mk->imesh_instance()),
  mb (mk->moab_instance())
{
  m_nPlanar = 0; //default is 3D
  m_nLineNumber = 0;
  szComment = "!";
  MAXCHARS = 300;
  pi = M_PI;
  m_dRadialSize = -1.0;
  m_dAxialSize = -1.0;
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
  m_dMergeTol = 1e-4;
  m_GeomFile = "";
}

AssyMesher::~AssyMesher()
{}

bool AssyMesher::add_modelent(ModelEnt *model_ent)
{
  return MeshOp::add_modelent(model_ent);
}

void AssyMesher::setup_this()
{
  // populate the model entities based on the geometry
  mk_core()->populate_model_ents();

  // create a set that will hold names of pin materials
  std::set<std::string> pinMtrlsSet;

  // collect names of pin materials based on input file
  for (int pci = 1; pci <= m_nPincells; ++pci)
  {
    int numCyl = 0;
    m_Pincell(pci).GetNumCyl(numCyl);
    for (int cyli = 1; cyli <= numCyl; ++cyli)
    {
      int cylSizes = 0;
       m_Pincell(pci).GetCylSizes(cyli, cylSizes);
      if (cylSizes > 0)
      {
        CVector<std::string> cylMat(cylSizes);
        m_Pincell(pci).GetCylMat(cyli, cylMat);
        for (int cmai = 1; cmai <= cylMat.GetSize(); ++cmai)
        {
          for (int mtrlIndx = 1; mtrlIndx <= m_nAssemblyMat; ++mtrlIndx)
          {
            if (m_szAssmMatAlias(mtrlIndx) == cylMat(cmai))
            {
              pinMtrlsSet.insert(m_szAssmMat(mtrlIndx));
            }
          }
        }
      }
    }
  }

  // get a vector of all surfaces
  std::vector<iGeom::EntityHandle> allSurfs;
  iGeom::EntitySetHandle rootSetHandle = igeom->getRootSet();
  igeom->getEntities(rootSetHandle, iBase_FACE, allSurfs);

  // get the name tag
  iGeom::TagHandle nameTag;
  int tagSize;
  igeom->getTagHandle("NAME", nameTag);
  igeom->getTagSizeBytes(nameTag, tagSize);

  // collect surfaces that have names identifying them as the tops of pins
  std::vector<iGeom::EntityHandle> pinTopSurfs;
  for (unsigned int ei = 0; ei < allSurfs.size(); ++ei)
  {
    char* entName = new char[tagSize];
    entName[0] = 0;
    igeom->getData(allSurfs[ei], nameTag, entName);
    // TODO: error check for result == iBase_SUCCESS
    // it should always be true for current code of CGM
    // Other calls to igeom should do that too.  iGeom::Error result = ...
    size_t enLen = strlen(entName);
    char* endMatchChar = strchr(entName, '@');
    if (endMatchChar == NULL)
    {
      endMatchChar = &entName[enLen];
    }
    if ((endMatchChar - entName) > 4 &&
        strncmp(endMatchChar - 4, "_top", 4) == 0)
    {
      std::string matNameStr(entName, endMatchChar - 4);
      if (pinMtrlsSet.find(matNameStr) != pinMtrlsSet.end())
      {
        pinTopSurfs.push_back(allSurfs[ei]);
      }
    }

    delete[] entName;
  }

  // get the tag on the geometry that identifies the associated model entity
  iGeom::TagHandle meTag = mk_core()->igeom_model_tag();

  // sizing function for radial mesh size . . . MeshKit core will delete it
  SizingFunction* radialMeshSizePtr;
  if (m_dRadialSize <= 0)
  {
    radialMeshSizePtr = new SizingFunction(mk_core(), -1, 1);
  }
  else
  {
    radialMeshSizePtr = new SizingFunction(mk_core(), -1, m_dRadialSize);
  }
  int radialSizeIndex = radialMeshSizePtr->core_index();
  std::cout << "Radial mesh size: " << m_dRadialSize << std::endl;

  // gather pointers to model entities for pin top surfaces and
  // set the mesh size on the pin top surfaces
  MEntVector pinTopSurfMes;
  for (unsigned int ptsi = 0; ptsi < pinTopSurfs.size(); ++ptsi)
  {
    // Remark: iterator pattern may perform better that indexed loop here
    iGeom::EntityHandle geoPinTopHandle = pinTopSurfs[ptsi];
    ModelEnt* mePinTop;
    igeom->getData(geoPinTopHandle, meTag, &mePinTop);
    // TODO: check success
    pinTopSurfMes.push_back(mePinTop);
    mePinTop->sizing_function_index(radialSizeIndex);
  }

  mk_core()->insert_node(mk_core()->construct_meshop("CAMALPaver",
      pinTopSurfMes), (MeshOp*) this);
}

void AssyMesher::execute_this()
{
  std::cout << "Execute : start meshing the assembly" << std::endl;
//  Start doing the steps in .jou file: /MeshKit/rgg/io.cpp:routine:CreateCubitJournal()
//  AssyMesher
//  1. Find surfaces with names <pin_material>_top and set radial mesh size, also set mesh scheme to CAMALTriMesher or GRUMMP trimesher
//  2. Sweep the volumes with top surfaces of pins, use bottom surfaces of the pins if required: bottom surfaces are name as <pin_material_bot>
//  3. Find edges in surfaces with names <material_side> set size equal to edge length
//  4. After meshing assign block names for all pins, by filtering using volumes
//  5. Now mesh the top cutout portion using  CAMALTriMesher or GRUMMP trimesher (Report if this fails, if this fails change in radial mesh size and edge interval might be needed)
//  6. Now sweep and name blocks.
//  7. Create Neumann Sets

  // step 1


}

void AssyMesher::PrepareIO (int argc, char *argv[], std::string  TestDir)
// ---------------------------------------------------------------------------
// Function: Obtains file names and opens input/output files
// Input:    command line arguments
// Output:   none
// ---------------------------------------------------------------------------
{
  // set and open input output files
  bool bDone = false;
#ifdef HAVE_ACIS
#define EXTENSION ".sat";
#endif
#ifdef HAVE_OCC
#define EXTENSION ".brep";
#endif
  do{
      if (2 == argc) {
          m_InputFile = (std::string)argv[1] + ".inp";
          m_LogName = m_InputFile + ".log";
          m_GeomFile = (std::string)argv[1] + EXTENSION;
      }
      else if (1 == argc){
          m_LogFile << "\nRunning default case:\n" << std::endl;
          m_GeomFile = TestDir + "/" + (char *)DEFAULT_TEST_AM+ EXTENSION;
          m_InputFile = TestDir + "/" + (char *)DEFAULT_TEST_AM + ".inp";
          m_LogName = (std::string)DEFAULT_TEST_AM + ".log";
      }
      // open input file for reading
      m_FileInput.open (m_InputFile.c_str(), std::ios::in);
      if (!m_FileInput){
          m_LogFile << "Unable to open file: " << m_InputFile << std::endl;
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
      m_LogFile <<  "\t\t         Tool to generate assembly mesh      " << '\n';
      m_LogFile <<  "\t\t\t\tArgonne National Laboratory" << '\n';
      m_LogFile <<  "\t\t\t\t        2015         " << '\n';
      m_LogFile <<  "\t\t---------------------------------------------------------" << '\n';
      m_LogFile <<  "\nsee README file for using the program and details on various cards.\n"<< std::endl;

  }while (!bDone);

// Read AssyGen input file
  CParser Parse;
  int nCyl =0, nCellMat=0, nInputLines=0;
  std::string card, szVolId, szVolAlias;
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
          if (szInputString.substr(0,8) == "cylinder"){
            ++nCyl;
          }
          if (szInputString.substr(0,12) == "cellmaterial"){
            ++nCellMat;
          }
        }

	// set the sizes
	if (nCyl>0) {
	  if (nCellMat > 0) {
	    m_Pincell(i).SetCellMatSize(nCellMat);
	  }
	  m_Pincell(i).SetNumCyl(nCyl);
	}
	else if (nCyl ==0) {
          // used to be nInputLines > 0 . . . is it an error if
          // neither a cylinder nor a material line in the pin cell input?
	  if (nCellMat > 0)
	    m_Pincell(i).SetCellMatSize(nCellMat);
	}
	nCyl = 0;
	nCellMat = 0;
      }
    }
    // breaking condition
    if(szInputString.substr(0,3) == "end"){
      std::istringstream szFormatString (szInputString);
      break;
    }
  }

//  //ACIS ENGINE
//#ifdef HAVE_ACIS
//  //  if(m_szEngine == "acis"){
//  m_szGeomFile = m_InputFile+".sat";
//  //  }
//#elif defined(HAVE_OCC)
//  //  OCC ENGINE
//  //  if (m_szEngine == "occ"){
//  m_szG= m_InputFile+".stp";
//  //  }
//#endif
//  std::cout << "\no/p geometry file name: " <<  m_szGeomFile <<std::endl;

  //Rewind the input file
  m_FileInput.clear (std::ios_base::goodbit);
  m_FileInput.seekg (0L, std::ios::beg);
  m_nLineNumber = 0;
//  CParser Parse;
//  std::string card;

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
    else if (szInputString.substr(0,8) == "geometry"){
      std::string outfile;
      std::istringstream szFormatString (szInputString);
      szFormatString >> card >> outfile;
      if(strcmp (outfile.c_str(), "surface") == 0 || szFormatString.fail())
        m_nPlanar=1;
    }
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

	ReadPinCellData( i);
	std::cout << "\nread pincell " << i << std::endl;
      }
    }

    // 'yes' or 'no' for creating sidesets
    if (szInputString.substr(0,13) == "createsideset"){
      std::istringstream szFormatString (szInputString);
      szFormatString >> card >> m_szSideset;
      std::cout <<"--------------------------------------------------"<<std::endl;
    }
    // specify a merge tolerance value for cubit journal file
    if (szInputString.substr(0,14) == "mergetolerance"){
      std::istringstream szFormatString (szInputString);
      szFormatString >> card >> m_dMergeTol;
      std::cout <<"--------------------------------------------------"<<std::endl;
    }  // Handle mesh size inputs
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
      szFormatString >> card >> m_dAxialSize;
      if(m_dAxialSize < 0 || szFormatString.fail())
        IOErrorHandler(ENEGATIVE);
      std::cout <<"--------------------------------------------------"<<std::endl;

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
    if (szInputString.substr(0,3) == "end"){


//      if ( m_nJouFlag == 0){
//        // impring merge before saving
//        // Imprint_Merge();

//	// save .sat file
//	IBERRCHK(igeom->save(m_szGeomFile.c_str()), *igeom);
//	std::cout << "Normal Termination.\n"<< "Geometry file: " << m_szGeomFile << " saved." << std::endl;
//      }
      break;
    }
  }

  // Done reading now load file

  IBERRCHK(igeom->load(m_GeomFile.c_str()), *igeom);

}

///////////////////////////////////////////////////////////////////////////////////////////

void AssyMesher::ReadPinCellData ( int i)
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

	//set local array
	dVCylRadii.SetSize(nRadii);
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

	//set local array
	dVCylRadii.SetSize(nRadii);
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


void AssyMesher::IOErrorHandler (ErrorStates ECode) const
// ---------------------------------------------------------------------------
//! Function: Displays error messages related to input data \n
//! Input:    Error code \n
//! Output:   none \n
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
  else
    std::cerr << "Unknown error ...?";

  std::cerr << '\n' << "Error in input file line : " << m_nLineNumber;
  std::cerr << std::endl;
  exit (1);
}

} // namespace MeshKit
