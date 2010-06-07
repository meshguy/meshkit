/*********************************************
Feb,10
Reactor Geometry Generator
Argonne National Laboratory

file with input o/p funtions
*********************************************/
#include <sstream>
#include "nrgen.hpp"
#include "parser.hpp"
#include <fstream>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#define ARRAY_INOUT( A ) A.ptr(), &A.capacity(), &A.size()
#define ARRAY_IN( A ) &A[0], A.size()

#define DEFAULT_TEST_FILE "assygen_default"

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
  std::cout << "\nsee README file for using the program and details on various cards.\n"<< std::endl;
}

int CNrgen::PrepareIO (int argc, char *argv[])
// ---------------------------------------------------------------------------
// Function: Obtains file names and opens input/output files
// Input:    command line arguments
// Output:   none
// ---------------------------------------------------------------------------
{
  //set geometry engine 
  std::string engine;
  //ACIS ENGINE
  engine = "ACIS";
  iGeom_newGeom( 0, &geom, &err, 0 ); // this is default way of specifying ACIS engine

  // OCC ENGINE
  // engine = "OCC;
  //  iGeom_newGeom( ";engine=OCC", &geom, &err, 12 );

  // set and open input output files
  bool bDone = false;
  do{
    if (2 == argc) {
      m_szFile = argv[1];
      m_szInFile=m_szFile+".inp";
      m_szGeomFile = m_szFile+".sat";
      if(engine =="OCC"){
	m_szGeomFile = m_szFile+".brep";
      }
      m_szJouFile = m_szFile+".jou";
      m_szSchFile = m_szFile+".template.jou";
    }
    else {
      std::cerr << "Usage: " << argv[0] << " <input file> WITHOUT EXTENSION"<< std::endl;   
      m_szInFile = (char *)DEFAULT_TEST_FILE;
      m_szGeomFile = (char *)DEFAULT_TEST_FILE;
      m_szJouFile = (char *)DEFAULT_TEST_FILE;
      m_szFile =  (char *)DEFAULT_TEST_FILE;
      m_szInFile+=".inp";
      m_szGeomFile+=".sat";
      if(engine =="OCC"){
	m_szGeomFile = m_szFile+".brep";
      }
      m_szJouFile+=".jou";
      m_szSchFile = m_szFile+".template.jou";
      std::cout <<"  No file specified.  Defaulting to: " << m_szInFile
		<< "  " << m_szGeomFile << "  " << m_szJouFile << std::endl;
    }
    // open the file
    m_FileInput.open (m_szInFile.c_str(), std::ios::in); 
    if (!m_FileInput){
      std::cout << "Unable to open file" << std::endl;
      std::cout << "Usage: assygen <input filename WITHOUT EXTENSION>"<< std::endl;
      m_FileInput.clear ();
      exit(1);
    }
    else
      bDone = true; // file opened successfully
  } while (!bDone);
  std::cout << "\nEntered input file name: " <<  m_szInFile <<std::endl;

  // open the file
  do{
    m_FileOutput.open (m_szJouFile.c_str(), std::ios::out); 
    if (!m_FileOutput){
      std::cout << "Unable to open o/p file" << std::endl;
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
      std::cout << "Unable to open o/p file" << std::endl;
      m_SchemesFile.clear ();
      exit(1);
    }
    else
      bDone = true; // file opened successfully
  } while (!bDone);

  std::cout<<"o/p geometry file name: "<<m_szGeomFile << "\no/p Cubit journal file name: "<< m_szJouFile    
	   << "\nInfo: o/p file extension must correspond to geometry engine used to build cgm " << std::endl;
  return 0;
}

int CNrgen::CountPinCylinders ()
// -------------------------------------------------------------------------------------------
// Function: reads the input file to count the no. of cyl in a pincell, before the actual read
// Input:    none
// Output:   none
// -------------------------------------------------------------------------------------------
{
  CParser Parse;
  int nCyl =0, nCellMat=0, nInputLines=0;
  std::string card, szVolId, szVolAlias;
  // count the total number of cylinder commands in each pincell
  for(;;){
    if (!Parse.ReadNextLine (m_FileInput, m_nLineNumber, szInputString, 
			     MAXCHARS, szComment))
      IOErrorHandler (INVALIDINPUT);
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
    // breaking condition
    if(szInputString.substr(0,3) == "end"){
      std::istringstream szFormatString (szInputString);
      break;
    }
  }
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
      if( (strcmp (m_szGeomType.c_str(), "hexagonal") != 0) &&
	  (strcmp (m_szGeomType.c_str(), "rectangular") != 0) || szFormatString.fail())
	IOErrorHandler(EGEOMTYPE);
    }

    if (szInputString.substr(0,8) == "geometry"){
      std::string outfile;
      std::istringstream szFormatString (szInputString);
      szFormatString >> card >> outfile;
      if(strcmp (outfile.c_str(), "surface") == 0 || szFormatString.fail())
	m_nPlanar=1;
    }   
    if (szInputString.substr(0,9) == "materials"){

      std::cout << "getting assembly material data" << std::endl;
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
      }
    }   
    if( (szInputString.substr(0,10) == "dimensions") || 
	(szInputString.substr(0,4) == "duct") ){

      std::cout << "getting assembly dimensions" << std::endl;
      if(m_szGeomType =="hexagonal"){
	std::istringstream szFormatString (szInputString);
	m_dVXYAssm.SetSize(2); m_dVZAssm.SetSize(2);

	szFormatString >> card >> m_nDimensions 
		       >> m_dVXYAssm(1) >> m_dVXYAssm(2)
		       >> m_dVZAssm(1) >> m_dVZAssm(2);

	m_dVAssmPitch.SetSize(m_nDimensions); m_szMAlias.SetSize(m_nDimensions);

	assms.setSize(m_nDimensions); // setup while reading the problem size


	for (int i=1; i<=m_nDimensions; i++){
	  szFormatString >> m_dVAssmPitch(i);
	  if( m_dVAssmPitch(i) < 0 )
	    IOErrorHandler(ENEGATIVE);
	}

	for (int i=1; i<=m_nDimensions; i++){
	  szFormatString >> m_szMAlias(i);
	  if(strcmp (m_szMAlias(i).c_str(), "") == 0)  
	    IOErrorHandler(EALIAS);
	}
      }   
      if(m_szGeomType =="rectangular"){
	std::istringstream szFormatString (szInputString);
	m_dVXYAssm.SetSize(2); 
	m_dVZAssm.SetSize(2);
	szFormatString >> card >> m_nDimensions 
		       >> m_dVXYAssm(1) >> m_dVXYAssm(2)
		       >> m_dVZAssm(1) >> m_dVZAssm(2);
	if (szFormatString.fail())
	  IOErrorHandler(INVALIDINPUT);	  
	m_dVAssmPitchX.SetSize(m_nDimensions);
	m_dVAssmPitchY.SetSize(m_nDimensions);
	m_szMAlias.SetSize(m_nDimensions);
	assms.setSize(m_nDimensions); 

	for (int i=1; i<=m_nDimensions; i++){
	  szFormatString >> m_dVAssmPitchX(i) >> m_dVAssmPitchY(i);
	  if( m_dVAssmPitchX(i) < 0 || m_dVAssmPitchY(i) < 0 || szFormatString.fail())
	    IOErrorHandler(ENEGATIVE);
	}

	for (int i=1; i<=m_nDimensions; i++){
	  szFormatString >> m_szMAlias(i);
	  if(strcmp (m_szMAlias(i).c_str(), "") == 0 || szFormatString.fail())  
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
      if(m_nDimensions > 0){
	dTotalHeight = m_dVZAssm(2)-m_dVZAssm(1);
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
      err = CreateOuterCovering();
      ERRORR("Error in CreateOuterCovering", err);

      // subtract pins before save
      Subtract_Pins();   
      if(m_nPlanar ==1){
	err = Create2DSurf();
	ERRORR("Error in Create2DSurf", err);
      }
    }

    // section the assembly as described in section card
    if (szInputString.substr(0,7) == "section"){
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
    if (szInputString.substr(0,4) == "move"){
      std::cout << "Moving geometry .." << std::endl;
      double dX, dY, dZ;
      std::istringstream szFormatString (szInputString);
      szFormatString >> card >> dX >> dY >> dZ;
      err = Move_Assm(dX, dY, dZ);
      ERRORR("Error in Move_Assm", err);
      std::cout <<"--------------------------------------------------"<<std::endl;

    }
    // ceter the assembly 
    if (szInputString.substr(0,6) == "center"){
      std::cout << "Positioning assembly to center" << std::endl;
      err = Center_Assm();
      ERRORR("Error in Center_Assm", err);
      m_Centered = true;
      std::cout <<"--------------------------------------------------"<<std::endl;
    }
    // rotate the assembly if rotate card is specified
    if (szInputString.substr(0,6) == "rotate"){
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
    // Handle mesh size inputs
    if (szInputString.substr(0,14) == "radialmeshsize"){
      std::istringstream szFormatString (szInputString);
      szFormatString >> card >> m_dRadialSize;
      if(m_dRadialSize < 0 || szFormatString.fail())
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
    if (szInputString.substr(0,3) == "end"){
 
      // impring merge before saving
      Imprint_Merge();
 
      // save .sat file
      iGeom_save(geom, m_szGeomFile.c_str(), NULL, &err, m_szGeomFile.length() , 0);
      CHECK("Save to file failed.");
      std::cout << "Normal Termination.\n"<< "Geometry file: " << m_szGeomFile << " saved." << std::endl;
      break;
    }
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
  // variables
  int nColor;
  std::string color[21] = {" ", "thistle", "grey", "deepskyblue", "red", "purple",  "green",
                          "yellow", "royalblue", "magenta", "cyan", "lightsalmon", "springgreen",
			  "gold", "orange", "brown", "pink", "khakhi", "black", "aquamurine", "mediumslateblue"};
  int nSideset=0, i, j;
  std::string szGrp, szBlock, szSurfTop, szSurfBot, szSize, szSurfSide;
  double dHeight = 0.0, dMid = 0.0;
  if(m_nDimensions > 0){
    dHeight=  m_dVZAssm(2)-m_dVZAssm(1);
    dMid = (m_Centered ? 0.0 : m_dVZAssm(1) + dHeight/2.0);
  }

  // writing to schemes .jou 
  m_SchemesFile << "## This file is created by rgg program in MeshKit ##\n";
  m_SchemesFile << "##Schemes " << std::endl  ;
  m_SchemesFile << "#{CIRCLE =\"circle interval 1 fraction 0.8\"}" << std::endl;
  m_SchemesFile << "#{HOLE = \"hole rad_interval 2 bias 0.0\"}" << std::endl;
  m_SchemesFile << "#{PAVE = \"pave\"}" << std::endl;
  m_SchemesFile << "#{MAP = \"map\"}" << std::endl;
  m_SchemesFile << "#{SWEEP = \"sweep\"}" << std::endl;  
  m_SchemesFile << "## Dimensions" << std::endl;
  if(m_szGeomType == "hexagonal"){
    if(m_nDimensions > 0){
      m_SchemesFile << "#{PITCH =" << m_dVAssmPitch(m_nDimensions) << "}" << std::endl;
    }
  }
  else if(m_szGeomType == "rectangular"){
    if(m_nDimensions > 0){    
      m_SchemesFile << "#{PITCHX =" << m_dVAssmPitchX(m_nDimensions)<< "}" << std::endl;
      m_SchemesFile << "#{PITCHY =" << m_dVAssmPitchY(m_nDimensions) << "}" << std::endl;
    }
  }
  if( m_nPlanar ==0){
    m_SchemesFile << "#{Z_HEIGHT = " << dHeight << "}" << std::endl;
    m_SchemesFile << "#{Z_MID = " << dMid << "}" << std::endl;
    
  }
  m_SchemesFile << "##Set Mesh Sizes" << std::endl;

  // volume only
  if(m_nPlanar == 0){
    if (-1.0 == m_dAxialSize)
      m_SchemesFile << "#{AXIAL_MESH_SIZE = 0.1*Z_HEIGHT}" << std::endl;
    else
      m_SchemesFile << "#{AXIAL_MESH_SIZE = " << m_dAxialSize << "}" << std::endl;
  }

  if (-1.0 == m_dRadialSize) {
    if (m_szGeomType == "hexagonal")
      m_SchemesFile << "#{RADIAL_MESH_SIZE = 0.1*PITCH}" << std::endl;
    else
      m_SchemesFile << "#{RADIAL_MESH_SIZE = 0.02*0.5*(PITCHX+PITCHY)}" << std::endl;
  }
  else
    m_SchemesFile << "#{RADIAL_MESH_SIZE = " << m_dRadialSize << "}" << std::endl;

  // stuff common to both surface and volume
  m_FileOutput << "## This file is created by rgg program in MeshKit ##\n";
  m_FileOutput << "#User needs to specify mesh interval and schemes in this file\n#" << std::endl;
  m_FileOutput << "{include(\"" << m_szSchFile << "\")}" <<std::endl;
  m_FileOutput << "#" << std::endl;
  
  // import the geometry file
  m_FileOutput << "# Import geometry file " << std::endl;
  m_FileOutput << "import '" << m_szGeomFile <<"'" <<std::endl;
  m_FileOutput << "#" << std::endl;

  // merge
  m_FileOutput << "#merge geometry" << std::endl; 
  m_FileOutput << "merge all" << std::endl;
  m_FileOutput << "#" << std::endl;

  // top surface sidesets
  m_FileOutput << "#Creating top surface sidesets" << std::endl; 
  for(int p=1;p<=m_szAssmMatAlias.GetSize();p++){
    szSurfTop = m_szAssmMat(p)+"_top";
    ++nSideset;
    m_FileOutput << "group 'tmpgrp' equals surface name \""  << szSurfTop  << "\"" << std::endl;
    m_FileOutput << "sideset " << nSideset << " surface in tmpgrp" << std::endl;    
  }
  m_FileOutput << "#" << std::endl;

  //surface only
  if(m_nPlanar ==1){ 

    // sideset curves on top surface creation dumps
    m_FileOutput << "#Creating side curve sidesets" << std::endl; 
    if(m_szGeomType =="hexagonal"){
      for(i=1; i<=6; i++){
	++nSideset;
	m_FileOutput << "group 'tmpgrp' equals curve name \"side_edge"  << i  << "\"" << std::endl;
	m_FileOutput << "sideset " << nSideset << " curve in tmpgrp" << std::endl;    
      }
    }
    if(m_szGeomType =="rectangular"){
      for(j=1; j<=4; j++){
	++nSideset;
	m_FileOutput << "group 'tmpgrp' equals curve name \"side_edge"  << j  << "\"" << std::endl;
	m_FileOutput << "sideset " << nSideset << " curve in tmpgrp" << std::endl;
      }
      m_FileOutput << "#" << std::endl;
    }

    // group creation dumps. each material surface  has a group
    m_FileOutput << "#Creating groups" << std::endl;  
    for(int p=1;p<=m_szAssmMatAlias.GetSize();p++){
      szGrp = "g_"+ m_szAssmMat(p);
      m_szAssmMat(p);
      m_FileOutput << "group \"" << szGrp << "\" add surface name \"" << m_szAssmMat(p) <<"\"" << std::endl;
    }
    m_FileOutput << "#" << std::endl;

    // block creation dumps
    m_FileOutput << "#Creating blocks, Note: you might need to combine some blocks" << std::endl; 
    for(int p=1;p<=m_szAssmMatAlias.GetSize();p++){
      szBlock = "b_"+ m_szAssmMat(p);
      szGrp = "g_"+ m_szAssmMat(p);
      m_FileOutput << "block " << p << " surface in " << szGrp  << std::endl;
      m_FileOutput << "block " << p << " name \"" << szBlock <<"\""<< std::endl;
    }
    m_FileOutput << "#" << std::endl;
  }

  // volume only
  else{ 

    // bottom surface sidesets
    m_FileOutput << "#Creating bot and side surface sidesets" << std::endl; 
    for(int p=1;p<=m_szAssmMatAlias.GetSize();p++){
      szSurfTop = m_szAssmMat(p)+"_bot";
      szSurfSide = m_szAssmMat(p)+"_side";  

      m_FileOutput << "#" << std::endl;

      ++nSideset;
      m_FileOutput << "group 'tmpgrp' equals surface name \""  << szSurfTop  << "\"" << std::endl;
      m_FileOutput << "sideset " << nSideset << " surface in tmpgrp" << std::endl;

      ++nSideset;
      m_FileOutput << "group 'tmpgrp' equals surface name \""  << szSurfSide  << "\"" << std::endl;
      m_FileOutput << "sideset " << nSideset << " surface in tmpgrp" << std::endl;
    }
    m_FileOutput << "#" << std::endl;

    // group creation dumps. each material surface  has a group
    for(int p=1;p<=m_szAssmMatAlias.GetSize();p++){
      szGrp = "g_"+ m_szAssmMat(p);
      m_szAssmMat(p);
      m_FileOutput << "group \"" << szGrp << "\" add body name \"" << m_szAssmMat(p) <<"\"" << std::endl;
    }
    m_FileOutput << "#" << std::endl;

    // block creation dumps
    m_FileOutput << "#Creating blocks, Note: you might need to combine some blocks" << std::endl; 
    for(int p=1;p<=m_szAssmMatAlias.GetSize();p++){
      szBlock = "b_"+ m_szAssmMat(p);
      szGrp = "g_"+ m_szAssmMat(p);
      m_FileOutput << "block " << p << " body in " << szGrp  << std::endl;
      m_FileOutput << "block " << p << " name \"" << szBlock <<"\""<< std::endl;
    }
    m_FileOutput << "#" << std::endl;

    //now set the sizes
    m_FileOutput << "#Set Meshing Scheme and Sizes, use template.jou to specify sizes" << std::endl; 

    for(int p=1;p<=m_szAssmMatAlias.GetSize();p++){
      szGrp = "g_"+ m_szAssmMat(p);
      szSize =  m_szAssmMat(p) + "_size";
      szSurfBot = m_szAssmMat(p) + "_bot";
      szSize =  m_szAssmMat(p) + "_surf_size";
      m_FileOutput << "group 'tmpgrp' equals surface name \""  << szSurfBot  << "\"" << std::endl;
      m_FileOutput << "surface in tmpgrp  size {"  << szSize <<"}" << std::endl;

      //       m_FileOutput << "vol in body in " << szGrp << " size  {"  << szSize <<"}" << std::endl;
      //       m_FileOutput << "vol in body in " << szGrp << " scheme {" << "SWEEP}"  << std::endl;
      //       m_FileOutput << "mesh vol in body in " << szGrp << "\n#" << std::endl;   

      // dumping these sizes schemes.jou also
      m_SchemesFile << "#{"  << szSize <<" = RADIAL_MESH_SIZE}" << std::endl;
    }  
    m_FileOutput << "#" << std::endl;
  }

  // some more common stuff meshing top surfaces set the sizes and mesh
  m_FileOutput << "#surfaces mesh, use template.jou to specify sizes" << std::endl; 
  for(int p=1;p<=m_szAssmMatAlias.GetSize();p++){
    szSurfTop = m_szAssmMat(p) + "_top";
    szGrp = "g_"+ m_szAssmMat(p);
    szSize =  m_szAssmMat(p) + "_surf_size";
    m_FileOutput << "group 'tmpgrp' equals surface name \""  << szSurfTop  << "\"" << std::endl;
    m_FileOutput << "surface in tmpgrp  size {"  << szSize <<"}" << std::endl;
    m_FileOutput << "surface in tmpgrp scheme {" << "PAVE" << "}"  << std::endl;
    //    m_FileOutput << "mesh surface in " << szGrp << "\n#" << std::endl;   

    // dumping these sizes schemes.jou also
    m_SchemesFile << "#{"  << szSize <<" = RADIAL_MESH_SIZE}" << std::endl;
  }  
  m_FileOutput << "#" << std::endl;

  // mesh all command after meshing surface
  m_FileOutput << "group 'tmpgrp' add surface name '_top'" << std::endl;
  m_FileOutput << "mesh tmpgrp" << std::endl;
  m_FileOutput << "#" << std::endl;
   
  if(m_nPlanar == 0){ // volumes only
    m_FileOutput << "surface with z_coord > {Z_MID -.1*Z_HEIGHT}" <<
      " and z_coord < {Z_MID + .1*Z_HEIGHT} size {AXIAL_MESH_SIZE}" << std::endl ;
    m_FileOutput << "mesh vol all" << std::endl;
  // color now
     for(int p=1;p<=m_szAssmMatAlias.GetSize();p++){
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
     for(int p=1;p<=m_szAssmMatAlias.GetSize();p++){
      szGrp = "g_"+ m_szAssmMat(p);
      if(p>20)
	nColor = 1;
      else 
	nColor = p;
      m_FileOutput << "color surface in " << szGrp << " " << color[nColor] << std::endl;
    }
  }   


  // save as .cub file dump
  m_FileOutput << "#\n#Save file" << std::endl; 
  std::string szSave = m_szFile + ".cub";
  m_FileOutput << "save as '"<< szSave <<"'" << " overwrite"<<std::endl; 


  std::cout << "Schemes file created: " << m_szSchFile << std::endl;
  std::cout << "Cubit journal file created: " << m_szJouFile << std::endl;
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
  int err2 = err;
  iGeom_getDescription( geom, buffer, &err2, sizeof(buffer) );
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
  // get the surface with max z
  double dTol=1.0e-3, dZTemp;
  int flag = 0, locTemp;
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
    if( (abs(min_corn[3*i+2]-max_corn[3*i+2])) < dTol ) {
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
    if( (abs(min_corn[3*i+2]-max_corn[3*i+2])) < dTol ) {
      continue; // its a max of min surface
    }
    else{ // its a side surface
      side_surf = surfs[i];
    }
    //set name for the sidesurf now
    if(side_surf !=0){
      ++nSide;
      sSideName = sMatName + "_side";
      os << sSideName << nSide;
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


int CNrgen::Center_Assm ()
// ---------------------------------------------------------------------------
// Function: centers all the entities along x and y axis
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
  double xmin, xmax, ymin, ymax, zmin, zmax;
  // position the assembly such that origin is at the center before sa
  iGeom_getBoundBox(geom,&xmin,&ymin,&zmin,
		    &xmax,&ymax,&zmax, &err);
  CHECK("Failed getting bounding box");

  // moving all geom entities to center      
  double xcenter = (xmin+xmax)/2.0;
  double ycenter = (ymin+ymax)/2.0;
  double zcenter = (zmin+zmax)/2.0;
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
// Function: rotates the whole assembly
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
  int nInputLines, nTempPin, t, nIFlag = 0;
  double dX = 0.0, dY =0.0, dZ=0.0;
  double  dP, dH, dSide, dHeight;
  iBase_EntityHandle assm = NULL;
  std::cout << "\ngetting Assembly data and creating ..\n"<< std::endl;
  std::istringstream szFormatString (szInputString);

  szFormatString >> card >> m_nPin;

  // width of the hexagon is n*n+1/2
  int nWidth =2*m_nPin -1;

  // creating a square array of size width
  m_Assembly.SetSize(nWidth, nWidth);

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
      szFormatString1 >> m_Assembly(m,n);
      if(szFormatString1.fail())
	IOErrorHandler (INVALIDINPUT);
      // if dummy pincell skip and continue
      if((m_Assembly(m,n)=="x")||(m_Assembly(m,n)=="xx")){
	continue;
      }
      // find that pincell
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

      m_Pincell(nTempPin).GetIntersectFlag(nIFlag);
      if(nIFlag){
	err = CreatePinCell_Intersect(nTempPin, dX, dY, dZ);
	ERRORR("Error in function CreatePinCell_Intersect", err);
      }
      else{
	err = CreatePinCell(nTempPin, dX, dY, dZ);
	ERRORR("Error in function CreatePinCell", err);
      }
    }
  }

  // get all the entities (in pins)defined so far, in an entity set - for subtraction later
  iGeom_getEntities( geom, root_set, iBase_REGION, ARRAY_INOUT(in_pins),&err );
  CHECK( "ERROR : getRootSet failed!" );

  if(m_nDimensions >0){

    // create outermost hexes
    std::cout << "\n\nCreating surrounding outer hexes .." << std::endl;
    for(int n=1;n<=m_nDimensions; n++){
      dSide = m_dVAssmPitch(n)/(sqrt(3));
      dHeight = m_dVZAssm(2) - m_dVZAssm(1);

      // creating coverings
      iGeom_createPrism(geom, dHeight, 6, 
			dSide, dSide,
			&assm, &err); 
      CHECK("Prism creation failed.");
	  
      // rotate the prism to match the pins
      iGeom_rotateEnt (geom, assm, 30, 0, 0, 1, &err);
      CHECK("Rotation failed failed.");

      m_Pincell(1).GetPitch(dP, dH);
      dX = m_nPin*dP;
      dY = (m_nPin-1)*dP*sqrt(3.0)/2.0;
      dZ = (m_dVZAssm(2) + m_dVZAssm(1))/2.0;

      // position the prism
      iGeom_moveEnt(geom, assm, dX,dY,dZ, &err);
      CHECK("Move failed failed.");  
      
      // populate the coverings array
      assms[n-1]=assm;
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
  int nInputLines, nTempPin, nIFlag = 0.0;
  double dX = 0.0, dY =0.0, dZ=0.0, dMoveX = 0.0, dMoveY = 0.0, dHeight = 0, dPX, dPY, dPZ;
  iBase_EntityHandle assm = NULL;

  std::istringstream szFormatString (szInputString);
  szFormatString >> card >> m_nPinX >> m_nPinY;
  m_Assembly.SetSize(m_nPinX,m_nPinY);

  //read the next line to get assembly info &store assembly info
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

      // if dummy pincell skip and continue
      if((m_Assembly(m,n)=="x")||(m_Assembly(m,n)=="xx")){
	m_Pincell(1).GetPitch(dPX, dPY, dPZ);
	dX+=dPX;
	continue;
      }	      
      // loop thro' all pins to get the type of pin
      for(int b=1; b<=m_nPincells; b++){
	m_Pincell(b).GetLineOne(szVolId, szVolAlias, nInputLines);
	if(m_Assembly(m,n) == szVolAlias)
	  nTempPin = b;
      }

      //now compute the location where the pin needs to be placed
      err = ComputePinCentroid(nTempPin, m_Assembly, m, n, dX, dY, dZ);
      ERRORR("Error in function ComputePinCentroid", err);

      // now create the pincell in the location found
      std::cout << "\n--------------------------------------------------"<<std::endl;
      std::cout << " m = " << m <<" n = " << n << std::endl;
      std::cout << "creating pin: " << nTempPin;
      std::cout << " at X Y Z " << dX << " " << dY << " " << dZ << std::endl;

      m_Pincell(nTempPin).GetIntersectFlag(nIFlag);
      if(nIFlag){
	err = CreatePinCell_Intersect(nTempPin, dX, dY, dZ);
	ERRORR("Error in function CreatePinCell_Intersect", err);
      }
      else{
	err = CreatePinCell(nTempPin, dX, dY, dZ);
	ERRORR("Error in function CreatePinCell", err);
      }
      // dMoveX and dMoveY are stored for positioning the outer squares later
      if(m == m_nPinY && n ==m_nPinX){
	dMoveX = dX/2.0;
	dMoveY = dY/2.0;
      }
    }
  }
  std::cout << "\n--------------------------------------------------"<<std::endl;

  // get all the entities (in pins)defined so far, in an entity set - for subtraction later
  iGeom_getEntities( geom, root_set, iBase_REGION, ARRAY_INOUT(in_pins),&err );
  CHECK( "ERROR : getRootSet failed!" );

	
  if(m_nDimensions > 0){
	
    // create outermost rectangular blocks
    std::cout << "\nCreating surrounding outer blocks .." << std::endl;
    for(int n=1;n<=m_nDimensions; n++){
      dHeight = m_dVZAssm(2) - m_dVZAssm(1);
      iGeom_createBrick(geom, m_dVAssmPitchX(n),  m_dVAssmPitchY(n),dHeight, 
			&assm, &err); 
      CHECK("Prism creation failed.");	  
      
      // position the outer block to match the pins	  
      dX = m_dVAssmPitchX(n)/4.0;
      dY =  m_dVAssmPitchY(n)/4.0;
      dZ = (m_dVZAssm(2) + m_dVZAssm(1))/2.0;
      std::cout << "Move " <<   dMoveX << " " << dMoveY <<std::endl;
      iGeom_moveEnt(geom, assm, dMoveX,dMoveY,dZ, &err);
      CHECK("Move failed failed.");  

      // populate the outer covering array squares
      assms[n-1]=assm;
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
    int tag_no;
    for(int p=1;p<=m_szAssmMatAlias.GetSize();p++){
      if(strcmp ( m_szMAlias(1).c_str(), m_szAssmMatAlias(p).c_str()) == 0){
	sMatName =  m_szAssmMat(p);
	tag_no=p;
      }
    }
    std::cout << "\ncreated innermost block: " << sMatName << std::endl;
	
    tmp_vol = assms[0];
    iGeom_setData(geom, tmp_vol, this_tag,
		  sMatName.c_str(), sMatName.size(), &err);
    CHECK("setData failed");

    err = Name_Faces(sMatName, tmp_vol, this_tag);
    ERRORR("Error in function Name_Faces", err);

    //  Naming outermost block edges - sidesets in cubit journal file
    std::cout << "Naming outermost block edges" << std::endl;	
    SimpleArray<iBase_EntityHandle> edges;
    iGeom_getEntAdj( geom, assms[m_nDimensions-1] , iBase_EDGE,ARRAY_INOUT(edges),
		     &err );
    CHECK( "ERROR : getEntAdj failed!" );

    // get the top corner edges of the outer most covering
    int count =0;//index for edge names
    std::ostringstream os;
    for (int i = 0; i < edges.size(); ++i){   
      iGeom_getEntBoundBox(geom, edges[i],&xmin,&ymin,&zmin,
			   &xmax,&ymax,&zmax, &err);
      CHECK("getEntBoundBox failed."); 
      if(zmax==zmin){
	if(zmax==m_dVZAssm(2)){

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

    // now subtract the outermost hexes and name them
    for(int n=m_nDimensions; n>1 ; n--){
      if(n>1){ 

	// copy cyl before subtract 
	iGeom_copyEnt(geom, assms[n-2], &tmp_vol, &err);
	CHECK("Couldn't copy inner duct wall prism.");
	    
	// subtract outer most cyl from brick
	iGeom_subtractEnts(geom, assms[n-1], tmp_vol, &tmp_new, &err);
	CHECK("Subtract of inner from outer failed.");	  
	    
	assms[n-1]=tmp_new;
	    
	// name the vols by searching for the full name of the abbreviated Cell Mat
	for(int p=1;p<=m_szAssmMatAlias.GetSize();p++){
	  if(strcmp ( m_szMAlias(n).c_str(), m_szAssmMatAlias(p).c_str()) == 0){
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
  return 0;
}

int CNrgen::Subtract_Pins()
// ---------------------------------------------------------------------------
// Function: subtract the pins from the outer block
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
  if (m_nDimensions >0 && in_pins.size()>0){
    iBase_EntityHandle unite= NULL, tmp_vol = NULL, tmp_new1 = NULL;
    SimpleArray<iBase_EntityHandle> copy_inpins(in_pins.size());
    tmp_vol = assms[0];

    // subtract the innermost hex from the pins
    std::cout << "Subtracting all pins from assembly .. " << std::endl;

    // make a copy of the pins
    if(in_pins.size() > 1){
      for (int i=0; i<in_pins.size();i++){
	iGeom_copyEnt(geom, in_pins[i], &copy_inpins[i], &err);
	CHECK("Couldn't copy inner duct wall prism.");		  
      }					
      // unite all pins
      iGeom_uniteEnts(geom,ARRAY_IN(in_pins), &unite, &err);
      CHECK( "uniteEnts failed!" );	  

      // subtract inner covering with united pins
      iGeom_subtractEnts(geom, tmp_vol,unite, &tmp_new1, &err);
      CHECK("Couldn't subtract pins from block.");
    }
    else{ // only one pin in in_pins
      iGeom_copyEnt(geom, in_pins[0], &copy_inpins[0], &err);
      iGeom_subtractEnts(geom, tmp_vol, in_pins[0], &tmp_new1, &err);
    }
  }
  std::cout << "\n--------------------------------------------------"<<std::endl;
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
  
  //   // merge tolerance
  //   double dTol = 1e-4;
  //   // now  merge
  //   std::cout << "\n\nMerging...." << std::endl;
  //   iGeom_mergeEnts(geom, ARRAY_IN(entities), dTol, &err);
  //   CHECK("Merge failed.");
  //   std::cout <<"merging finished."<< std::endl;
  //   std::cout << "\n--------------------------------------------------"<<std::endl;
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
  double dtol = m_dVZAssm(2);
  for (int i = 0; i < surfs.size(); ++i){ 
    if((max_corn[3*i+2] ==  dtol) && (min_corn[3*i+2] ==  dtol))
      t++;
  }
  
  // allocate arrays
  SimpleArray<iBase_EntityHandle> max_surfs(t);
  SimpleArray<iBase_EntityHandle> new_surfs(t);
  t=0;

  // store the max surfaces in max_surfs
  for (int i = 0; i < surfs.size(); ++i){

    // locate surfaces for which max and min zcoord is same as maxz coord
    if((max_corn[3*i+2] ==  dtol) && (min_corn[3*i+2] ==  dtol)){
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
  double zcenter = m_dVZAssm(2)/2.0;//move up
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

  // name tag handle
  iBase_TagHandle this_tag= NULL;
  char* tag_name = (char*)"NAME";

  std::string sMatName = "";
  std::string sMatName1 = "";

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

      dHeight = dVEndZ(n) - dVStartZ(n);  

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
      
      dZMove = (dVEndZ(n)+dVEndZ(n-1))/2.0;

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

	Name_Faces(sMatName, cell, this_tag);
	CHECK("Name_Faces failed");
      }
      // loop and create cylinders
      if(nCyl > 0){
	m_Pincell(i).GetCylSizes(n, nRadii);
	SimpleArray<iBase_EntityHandle> cyls(nRadii);

	//declare variables
	CVector<double> dVCylRadii(nRadii);
	CVector<std::string> szVMat(nRadii);
	CVector<std::string> szVCylMat(nRadii);
 
	//get values
	m_Pincell(i).GetCylRadii(n, dVCylRadii);
	m_Pincell(i).GetCylPos(n, dVCylXYPos);
	m_Pincell(i).GetCylMat(n, szVCylMat);
	m_Pincell(i).GetCylZPos(n, dVCylZPos);
	dHeight = dVCylZPos(2)-dVCylZPos(1);

	for (int m=1; m<=nRadii; m++){ 
 
	  iGeom_createCylinder(geom, dHeight, dVCylRadii(m), dVCylRadii(m),
			       &cyl, &err);
	  CHECK("Couldn't create fuel rod.");

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
	  iGeom_subtractEnts(geom, cells[n-1], cyls[nRadii-1], &tmp_new, &err);
	  CHECK("Subtract of inner from outer failed.");
	
	  // copy the new into the cyl array
	  cells[n-1] = tmp_new; cell = tmp_new;
	  cyls[nRadii-1]=tmp_vol;

	}
	//set tag on inner most cylinder, search for the full name of the abbreviated Cell Mat
	for(int p=1;p<=m_szAssmMatAlias.GetSize();p++){
	  if(strcmp (szVCylMat(1).c_str(), m_szAssmMatAlias(p).c_str()) == 0){
	    sMatName = m_szAssmMat(p);
	  }
	}
	tmp_vol1=cyls[0]; //inner most cyl

	iGeom_setData(geom, tmp_vol1, this_tag,
		      sMatName.c_str(), 10, &err);
	CHECK("setData failed");
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
	  int tag_no;
	  for(int p=1;p<=m_szAssmMatAlias.GetSize();p++){
	    if(strcmp (szVCylMat(b).c_str(), m_szAssmMatAlias(p).c_str()) == 0){
	      tag_no = p;
	      sMatName =  m_szAssmMat(p);
	    }
	  }
	  std::cout << "created: " << sMatName << std::endl;
	  // set the name of the annulus
	  iGeom_setData(geom, tmp_new, this_tag,
			sMatName.c_str(),sMatName.size(), &err);
	  CHECK("setData failed");

	  err= Name_Faces(sMatName, tmp_new, this_tag);
	  ERRORR("Error in function Name_Faces", err);

	  // copy the new into the cyl array
	  cyls[b-1] = tmp_new;
	}
      }
    }
  }
  // this branch of the routine is responsible for creating cylinders with '0' cells
  if(nCells == 0){ 
 
    // get cylinder data
    m_Pincell(i).GetNumCyl(nCyl);
    nCells = nCyl;


    for(int n=1;n<=nCells; n++){

 
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
	CVector<double> dVCylRadii(nRadii);
	CVector<std::string> szVMat(nRadii);
	CVector<std::string> szVCylMat(nRadii);
 
	//get values
	m_Pincell(i).GetCylRadii(n, dVCylRadii);
	m_Pincell(i).GetCylPos(n, dVCylXYPos);
	m_Pincell(i).GetCylMat(n, szVCylMat);
	m_Pincell(i).GetCylZPos(n, dVCylZPos);

	dHeight = dVCylZPos(2)-dVCylZPos(1);

	for (int m=1; m<=nRadii; m++){ 
 
	  iGeom_createCylinder(geom, dHeight, dVCylRadii(m), dVCylRadii(m),
			       &cyl, &err);
	  CHECK("Couldn't create fuel rod.");

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

	iGeom_setData(geom, tmp_vol1, this_tag,
		      sMatName.c_str(), 10, &err);
	CHECK("setData failed");

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
	  std::cout << "created: " << sMatName << std::endl;
	  // set the name of the annulus
	  iGeom_setData(geom, tmp_new, this_tag,
			sMatName.c_str(),sMatName.size(), &err);
	  CHECK("setData failed");
	  err= Name_Faces(sMatName, tmp_new, this_tag);
	  ERRORR("Error in function Name_Faces", err);

	  // copy the new into the cyl array
	  cyls[b-1] = tmp_new;
	}
      }
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

  // name tag handle
  iBase_TagHandle this_tag= NULL;
  char* tag_name = (char*)"NAME";

  std::string sMatName = "";
  std::string sMatName0 = "";
  std::string sMatName1 = "";

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

      dHeight = dVEndZ(n) - dVStartZ(n);  

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
      
      dZMove = (dVEndZ(n)+dVEndZ(n-1))/2.0;

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
	CVector<double> dVCylRadii(nRadii);
	CVector<std::string> szVMat(nRadii);
	CVector<std::string> szVCylMat(nRadii);
	
	//get values
	m_Pincell(i).GetCylRadii(n, dVCylRadii);
	m_Pincell(i).GetCylPos(n, dVCylXYPos);
	m_Pincell(i).GetCylMat(n, szVCylMat);
	m_Pincell(i).GetCylZPos(n, dVCylZPos);
	dHeight = dVCylZPos(2)-dVCylZPos(1);

	for (int m=1; m<=nRadii; m++){ 
 
	  iGeom_createCylinder(geom, dHeight, dVCylRadii(m), dVCylRadii(m),
			       &cyl, &err);
	  CHECK("Couldn't create fuel rod.");

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
	iGeom_setData(geom, tmp_vol1, this_tag,
		      sMatName.c_str(), 10, &err);
	CHECK("setData failed");
	    
	err= Name_Faces(sMatName, tmp_vol1, this_tag);
	ERRORR("Error in function Name_Faces", err);
	 
	// copy the outermost cyl
	iGeom_copyEnt(geom, intersec_main[nRadii-1], &tmp_intersec, &err);
	CHECK("Couldn't copy inner duct wall prism."); 

	// subtract the outermost cyl from the cell
	iGeom_subtractEnts(geom, cells[n-1], tmp_intersec, &tmp_new, &err);
	CHECK("Subtract of inner from outer failed.");
	// now search for the full name of the abbreviated Cell Mat
	int tag_no;
	for(int p=1;p<=m_szAssmMatAlias.GetSize();p++){
	  if(strcmp (szVCellMat(n).c_str(), m_szAssmMatAlias(p).c_str()) == 0){
	    tag_no = p;
	    sMatName =  m_szAssmMat(p);
	  }
	}
	std::cout << "created: " << sMatName << std::endl;
	// set the name of the annulus
	iGeom_setData(geom, tmp_new, this_tag,
		      sMatName.c_str(),sMatName.size(), &err);
	CHECK("setData failed");
	err = Name_Faces(sMatName, tmp_new, this_tag);
	ERRORR("Error in function Name_Faces", err);

 	for (int b=nRadii; b>1; b--){  

 	  iGeom_copyEnt(geom, intersec_main[b-2], &tmp_intersec, &err);
 	  CHECK("Couldn't copy inner duct wall prism."); 
 	  //subtract tmp vol from the outer most
	  iGeom_subtractEnts(geom, intersec_main[b-1], tmp_intersec, &tmp_new, &err);
	  CHECK("Subtract of inner from outer failed.");

	  // now search for the full name of the abbreviated Cell Mat
	  int tag_no;
	  for(int p=1;p<=m_szAssmMatAlias.GetSize();p++){
	    if(strcmp (szVCylMat(b).c_str(), m_szAssmMatAlias(p).c_str()) == 0){
	      tag_no = p;
	      sMatName =  m_szAssmMat(p);
	    }
	  }
	  std::cout << "created: " << sMatName << std::endl;
	  // set the name of the annulus
	  iGeom_setData(geom, tmp_new, this_tag,
			sMatName.c_str(),sMatName.size(), &err);
	  CHECK("setData failed");
	  err = Name_Faces(sMatName, tmp_new, this_tag);
	  ERRORR("Error in function Name_Faces", err);

	  // copy the new into the cyl array
	  cyls[b-1] = tmp_new;

	}
      }
    }
  }
  // this branch of the routine is responsible for creating cylinders with '0' cells
  if(nCells == 0){ 
 
    // get cylinder data
    m_Pincell(i).GetNumCyl(nCyl);
    nCells = nCyl;

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
	CVector<double> dVCylRadii(nRadii);
	CVector<std::string> szVMat(nRadii), szVCylMat(nRadii);
	
	//get values
	m_Pincell(i).GetCylRadii(n, dVCylRadii);
	m_Pincell(i).GetCylPos(n, dVCylXYPos);
	m_Pincell(i).GetCylMat(n, szVCylMat);

	for (int m=1; m<=nRadii; m++){ 
 
	  iGeom_createCylinder(geom, dHeight, dVCylRadii(m), dVCylRadii(m),
			       &cyl, &err);
	  CHECK("Couldn't create fuel rod.");

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
	iGeom_setData(geom, tmp_vol1, this_tag,
		      sMatName.c_str(), 10, &err);
	CHECK("setData failed");
	    
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
	  int tag_no;
	  for(int p=1;p<=m_szAssmMatAlias.GetSize();p++){
	    if(strcmp (szVCylMat(b).c_str(), m_szAssmMatAlias(p).c_str()) == 0){
	      tag_no = p;
	      sMatName =  m_szAssmMat(p);
	    }
	  }
	  std::cout << "created: " << sMatName << std::endl;
	  // set the name of the annulus
	  iGeom_setData(geom, tmp_new, this_tag,
			sMatName.c_str(),sMatName.size(), &err);
	  CHECK("setData failed");
	  err = Name_Faces(sMatName, tmp_new, this_tag);
	  ERRORR("Error in function Name_Faces", err);

	  // copy the new into the cyl array
	  cyls[b-1] = tmp_new;

	}
      }

    }
  }
  return 0;
}
