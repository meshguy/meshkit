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
  m_nPlanar = 0; //default is 3D
  m_nLineNumber = 0;
  geom = NULL;
  root_set= NULL;
  szComment = "!";
  MAXCHARS = 300;
  pi = M_PI;
  m_dRadialSize = -1.0;
  m_dAxialSize = -1.0;
  m_nDimensions = 0;
  m_nMaterialSetId = 1;
  m_nNeumannSetId = 1;
  m_szEngine = "acis";
  m_nDuct = 0;
  m_nDuctNum = 0;
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
