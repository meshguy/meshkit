/*********************************************
Utility Library Function
*********************************************/
#include <iostream>
#include "arraybase.hpp"

double CArrayBase::m_dAllocated = 0.0;
double CArrayBase::m_dDeAllocated = 0.0;

CArrayBase::CArrayBase ()
// ----------------------------------------------------------------------------
// Function: constructor
// Input:    None 
// Output:   None
// ----------------------------------------------------------------------------
{
}

CArrayBase::~CArrayBase ()
// ----------------------------------------------------------------------------
// Function: destructor
// Input:    None 
// Output:   None
// ----------------------------------------------------------------------------
{
}

void CArrayBase::ShowStatistics ()
// ----------------------------------------------------------------------------
// Function: destructor
// Input:    None 
// Output:   None
// ----------------------------------------------------------------------------
{
    std::cout << "  Allocated : " << m_dAllocated   << " bytes\n"
              << "DeAllocated : " << m_dDeAllocated << " bytes\n";
}
