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
#include "constants.hpp"
#include "vectortemplate.hpp"
#include "matrixtemplate.hpp"
#include "pincell.hpp"

#include "iGeom.h"

class CNrgen
{
public:
  CNrgen ();    // ctor
  ~CNrgen ();   // dtor
  enum ErrorStates {PINCELLS, INVALIDINPUT};

  // helper functions
  void Banner (std::ostream& OF);
  void PrepareIO (int argc, char *argv[]);
  void CountPinCylinders ();
  int  ReadAndCreate ();
  int Get_Max_Surf (iBase_EntityHandle vol, iBase_EntityHandle* max_surf);
  void ReadPinCellData(int i);

  int CreatePinCell(int i, double dX,
		    double dY, double dZ);
  void ComputePinCentroid(int, CMatrix<std::string>, int, int,
			  double&, double&, double&);
  int TerminateProgram ();

  // modifier functions
  void SetSize ();

  // variables
  int m_nLineNumber, m_nPlanar;	   // current line number in input file
  std::ifstream m_FileInput;	// File Input
  std::ofstream m_FileOutput;	// File Output
  std::string m_szFileNameOut;
  std::string m_szFileNameIn;
private:
  CMatrix<std::string> m_Assembly;
  std::string m_szGeomType;
  int m_nAssemblyMat, m_nDimensions, m_nPincells , m_nAssmVol;
  CVector<double> m_dVAssmPitch, m_dVAssmPitchX, m_dVAssmPitchY;;
  CVector<double> m_dVXYAssm, m_dVZAssm, m_dAssmSize;
  CVector<std::string> m_szAssmMat, m_szAssmMatAlias, m_szMAlias;

  CVector<CPincell> m_Pincell;
  int m_nPin, m_nPinX, m_nPinY;
  double m_dPitch;


  // error handlers
  void ErrorHandler (int) const; 
  void IOErrorHandler (ErrorStates) const;
  friend class CPincell;
};

#endif
