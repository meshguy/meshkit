/*********************************************
Dec,09
Reactor Geometry Generator
Argonne National Laboratory

Pin Cylinder class definition.
*********************************************/
#include "cylinder.hpp"

CCylinder::CCylinder ()
// ---------------------------------------------------------------------------
// Function: default constructor
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
	m_dVXYPos.SetSize(2);
	m_dVZPos.SetSize(2);

}

CCylinder::CCylinder (const CCylinder& NO)
// ---------------------------------------------------------------------------
// Function: copy constructor
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{

}

CCylinder::~CCylinder ()
// ---------------------------------------------------------------------------
// Function: destructor
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
}

void CCylinder::SetSizes(int nRadii)
// ---------------------------------------------------------------------------
// Function: sets the first line of pin input
// Input:    volume id of the pin, alias and total no. of lines in the pin input
// Output:   none
// ---------------------------------------------------------------------------
{
	m_nRadii = nRadii;
	m_szVMat.SetSize(nRadii);
	m_dVRadii.SetSize(nRadii);
}

void CCylinder::SetRadii(CVector<double> dVRadii)
// ---------------------------------------------------------------------------
// Function: sets the first line of pin input
// Input:    volume id of the pin, alias and total no. of lines in the pin input
// Output:   none
// ---------------------------------------------------------------------------
{
	m_dVRadii = dVRadii;
}

void CCylinder::SetPos(CVector<double> dVXYPos)
// ---------------------------------------------------------------------------
// Function: sets the first line of pin input
// Input:    volume id of the pin, alias and total no. of lines in the pin input
// Output:   none
// ---------------------------------------------------------------------------
{
	m_dVXYPos = dVXYPos;
}


void CCylinder::SetZPos(CVector<double> dVZPos)
// ---------------------------------------------------------------------------
// Function: sets the first line of pin input
// Input:    volume id of the pin, alias and total no. of lines in the pin input
// Output:   none
// ---------------------------------------------------------------------------
{
	m_dVZPos = dVZPos;
}

void CCylinder::SetMat(CVector<std::string> szVMat)
// ---------------------------------------------------------------------------
// Function: sets the first line of pin input
// Input:    volume id of the pin, alias and total no. of lines in the pin input
// Output:   none
// ---------------------------------------------------------------------------
{
	m_szVMat = szVMat;
}

void CCylinder::GetSizes(int &nRadii)
// ---------------------------------------------------------------------------
// Function: sets the first line of pin input
// Input:    volume id of the pin, alias and total no. of lines in the pin input
// Output:   none
// ---------------------------------------------------------------------------
{
	nRadii = m_nRadii;
}

void CCylinder::GetRadii(CVector<double> &dVRadii)
// ---------------------------------------------------------------------------
// Function: sets the first line of pin input
// Input:    volume id of the pin, alias and total no. of lines in the pin input
// Output:   none
// ---------------------------------------------------------------------------
{
	dVRadii = m_dVRadii;
}

void CCylinder::GetMat(CVector<std::string> &szVMat)
// ---------------------------------------------------------------------------
// Function: sets the first line of pin input
// Input:    volume id of the pin, alias and total no. of lines in the pin input
// Output:   none
// ---------------------------------------------------------------------------
{
	szVMat = m_szVMat;
}

void CCylinder::GetPos(CVector<double> &dVCoor)
// ---------------------------------------------------------------------------
// Function: sets the first line of pin input
// Input:    volume id of the pin, alias and total no. of lines in the pin input
// Output:   none
// ---------------------------------------------------------------------------
{
	dVCoor = m_dVXYPos;
}

void CCylinder::GetZPos(CVector<double> &dVZPos)
// ---------------------------------------------------------------------------
// Function: sets the first line of pin input
// Input:    volume id of the pin, alias and total no. of lines in the pin input
// Output:   none
// ---------------------------------------------------------------------------
{
	dVZPos = m_dVZPos;
}
