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
#include <vector>
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
  enum ErrorStates {PINCELLS, INVALIDINPUT, EMAT, EGEOMTYPE, EGEOMENGINE, ENEGATIVE, EALIAS, EPIN, EUNEQUAL};

  void Banner (std::ostream& OF);
  int PrepareIO (int argc, char *argv[]);
  int ReadInputPhase1 ();
  int ReadCommonInp ();
  int ReadAndCreate ();
  int Name_Faces(const std::string sMatName, const iBase_EntityHandle body, 
		 iBase_TagHandle this_tag);
  int Center_Assm(char&);
  int Section_Assm (char&, double&, const std::string);
  int Rotate_Assm (char&, double&);
  int Move_Assm (double&,double&,double&);
  int Create_HexAssm(std::string &);
  int Create_CartAssm(std::string &);
  int CreateOuterCovering();
  int CreateAssyGenInputFiles();
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
  bool have_common;

private:
  
  // number of sides in the geometry
  int com_run_count;
  int m_nSides;
  bool m_bCreateMatFiles, save_exodus;
  // file Input
  std::ifstream m_FileInput, m_FileCommon;
    
  // journal file Output
  std::ofstream m_FileOutput, m_SchemesFile, m_AssmInfo;

  // string for file names
  std::string m_szFile, m_szInFile, m_szCommonFile, m_szGeomFile,m_szJouFile, m_szSchFile, m_szAssmInfo, m_szInfo, m_szSmooth, m_szLogFile, m_szMeshScheme;

   std::vector< std::vector<iBase_EntityHandle> > cp_inpins;

   std::vector<std::string> m_szDuctMats;

// matrix for holding pincell arrangement
  CMatrix<std::string> m_Assembly; 

  // matrix for holding verts coordinates used in tet-meshing
  CMatrix<double> m_dMTopSurfCoords; 

  // vector for duct specification 
  CMatrix<double> m_dMAssmPitch, m_dMAssmPitchX, m_dMAssmPitchY, m_dMXYAssm, m_dMZAssm;
  
  // vector for material names
  CVector<std::string> m_szAssmMat, m_szAssmMatAlias;
  CVector<std::string> m_szBLAssmMat;

  CVector<double> m_dAxialSize, m_dBLMatBias;
  CVector<int> m_nListMatSet, m_nListNeuSet, m_nBLMatIntervals;

  CMatrix<std::string> m_szMMAlias;  

  // vector holding a pincell
  CVector<CPincell> m_Pincell; 

  struct superblocks{
      int m_nSuperBlockId;
      std::string m_szSuperBlockAlias;
      int m_nNumSBContents;
      CVector<int> m_nSBContents;
  };
  int m_nSuperBlocks;
  CVector<superblocks> sb;
  int tmpSB;
  // string for geomtype, engine, meshtype
  std::string m_szEngine, m_szInnerDuct;
  std::string m_szGeomType;       
  std::string m_szMeshType;
  std::string m_szSideset; 
  int m_nAssyGenInputFiles;
  std::string pin_name;
  // integers for vectors sizes, err etc
  int m_nBLAssemblyMat, m_nAssemblyMat, m_nDimensions, m_nPincells , m_nAssmVol, m_nPin, m_nPinX, m_nPinY, err, m_nLineNumber, m_nPlanar,
    m_nNeumannSetId, m_nMaterialSetId, m_nDuct, m_nDuctNum, m_nJouFlag, m_nTotalPincells; 
  int m_edgeInterval, m_nStartpinid, m_nHblock;
  // doubles for pincell pitch, pi and mesh sizes resp.
  double m_dPitch, pi, m_dRadialSize, m_dTetMeshSize, m_dMergeTol, m_dZstart, m_dZend;
 
  // igeom related
  SimpleArray<iBase_EntityHandle> assms, in_pins;
  iGeom_Instance geom;
  iBase_EntitySetHandle root_set;
  
  // parsing related
  std::string szInputString;
  std::string szComment;
  int MAXCHARS, MAXLINES;
  
  // error handlers
  void IOErrorHandler (ErrorStates) const;
  friend class CPincell;
};

#endif
