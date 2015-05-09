/*********************************************
Jan,10
Reactor Geometry Generator
Argonne National Laboratory

CNrgen class definition.
*********************************************/
#include "nrgen.hpp"
/* ==================================================================
   ======================= CNrgen class =============================
   ================================================================== */

CNrgen::CNrgen ()
// ---------------------------------------------------------------------------
// Function: default constructor
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
  err = 0;
  tmpSB = 1;
  m_nPlanar = 0; //default is 3D
  m_nLineNumber = 0;
  geom = NULL;
  root_set= NULL;
  szComment = "!";
  MAXCHARS = 1500;
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
  save_exodus = false;
  have_common = true;
  com_run_count = 0;
  m_nBLAssemblyMat = 0;
}

CNrgen::~CNrgen ()
// ---------------------------------------------------------------------------
// Function: destructor
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
  // deallocate ... deallocate ... deallocate
}
