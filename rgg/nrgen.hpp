/*********************************************
Dec,09
Reactor Geometry Generator
Argonne National Laboratory

CNrgen class definition.
*********************************************/
#ifndef __RGG_FRAME_H__
#define __RGG_FRAME_H__

#include <fstream>
#include <iostream>
#include <math.h>
#include "vectortemplate.hpp"
#include "matrixtemplate.hpp"
#include "pincell.hpp"
#include "SimpleArray.hpp"
#include "iGeom.h"

// helper macro for assygen
#define ERRORR(a,b) {if (0 != err) {std::cerr << a << std::endl; return b;}}

// helper macro for igeom
#define CHECK( STR ) if (err != iBase_SUCCESS) return Print_Error( STR, err, geom, __FILE__, __LINE__ )

class CNrgen
{
public:
  CNrgen ();    // ctor
  ~CNrgen ();   // dtor
  enum ErrorStates {PINCELLS, INVALIDINPUT, EMAT, EGEOMTYPE, EGEOMENGINE, ENEGATIVE, EALIAS, EPIN};

  void Banner (std::ostream& OF);
  int PrepareIO (int argc, char *argv[]);
  int ReadInputPhase1 ();
  int ReadAndCreate ();
  int Name_Faces(const std::string sMatName, const iBase_EntityHandle body, 
		   iBase_TagHandle this_tag);
  int Center_Assm();
  int Section_Assm (char&, double&, const std::string);
  int Rotate_Assm (char&, double&);
  int Move_Assm (double&,double&,double&);
  int Create_HexAssm(std::string &);
  int Create_CartAssm(std::string &);
  int CreateOuterCovering();
  int Imprint_Merge ();
  int Subtract_Pins ();
  int Create2DSurf();
  int ReadPinCellData(int i);
  int CreatePinCell_Intersect(int i, double dX,
		    double dY, double dZ);
  int CreatePinCell(int i, double dX,
		    double dY, double dZ);  
  int CreateCubitJournal();
  int ComputePinCentroid(int, CMatrix<std::string>, int, int,
			 double&, double&, double&);
  bool Print_Error(const char* desc, 
			 int err,
			 iGeom_Instance geom,
			 const char* file,
			 int line );
  int TerminateProgram ();

private:

  // file Input
  std::ifstream m_FileInput;  
    
  // journal file Output
  std::ofstream m_FileOutput, m_SchemesFile;

  // string for file names
  std::string m_szFile, m_szInFile, m_szGeomFile,m_szJouFile, m_szSchFile;    

  // matrix for holding pincell arrangement
  CMatrix<std::string> m_Assembly; 

  // vector for duct specification 
  CVector<double> m_dVAssmPitch, m_dVAssmPitchX, m_dVAssmPitchY, m_dVXYAssm, m_dVZAssm, m_dAssmSize;
  
  // vector for material names
  CVector<std::string> m_szAssmMat, m_szAssmMatAlias, m_szMAlias;  

  // vector holding a pincell
  CVector<CPincell> m_Pincell; 

  // string for geomtype, engine
  std::string m_szEngine;
  std::string m_szGeomType;       

  // integers for vectors sizes, err etc
  int m_nAssemblyMat, m_nDimensions, m_nPincells , m_nAssmVol, m_nPin, m_nPinX, m_nPinY, err, m_nLineNumber, m_nPlanar, m_nNeumannSetId, m_nMaterialSetId; 

  // doubles for pincell pitch, pi and mesh sizes resp.
  double m_dPitch, pi, m_dRadialSize, m_dAxialSize;      
 
  // bool for checking if assembly is centered
  bool m_Centered;

  // igeom related
  SimpleArray<iBase_EntityHandle> assms, in_pins;
  iGeom_Instance geom;
  iBase_EntitySetHandle root_set;
  
  // parsing related
  std::string szInputString;
  std::string szComment;
  int MAXCHARS;

  // error handlers
  void IOErrorHandler (ErrorStates) const;
  friend class CPincell;
};

#endif
