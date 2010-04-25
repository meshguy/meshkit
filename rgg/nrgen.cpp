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
  m_nPlanar = 0; //default is 3D
  geom = NULL;
  root_set= NULL;
  szComment = "!";
  MAXCHARS = 300;
  pi = M_PI;
  m_RadialSize = -1.0;
  m_AxialSize = -1.0;
  m_Centered = false;
  m_nDimensions = 0;
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

void CNrgen::SetSize ()
// ---------------------------------------------------------------------------
// Function: memory allocation for all major arrays in the program
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
}


void CNrgen::ErrorHandler (int nCode) const
// ---------------------------------------------------------------------------
// Function: displays error messages related to frame analysis
// Input:    error code
// Output:   none
// ---------------------------------------------------------------------------
{
    std::cerr << '\n';

    if (nCode == 1) 
        std::cerr << "This is an error";
    else if (nCode == 2) // error in matrix operation
        std::cerr << "Error in.. operation.";	
    else
        std::cerr << "Unknown error ...?";

    std::cerr << std::endl;
    // exit (1);
}
