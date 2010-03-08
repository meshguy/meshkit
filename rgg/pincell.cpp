/*********************************************
Reactor Geometry Generator
Argonne National Laboratory

CPincell class definition.
*********************************************/
#include "pincell.hpp"

CPincell::CPincell ()
// ---------------------------------------------------------------------------
// Function: default constructor
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
}

CPincell::CPincell (const CPincell& NO)
// ---------------------------------------------------------------------------
// Function: copy constructor
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
// TBC
}

CPincell::~CPincell ()
// ---------------------------------------------------------------------------
// Function: destructor
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
}


void CPincell::SetLineOne (std::string szVolId, std::string szVolAlias, int nInputLines)
// ---------------------------------------------------------------------------
// Function: sets the first line of pin input
// Input:    volume id of the pin, alias and total no. of lines in the pin input
// Output:   none
// ---------------------------------------------------------------------------
{
	m_szVolId = szVolId;
	m_szVolAlias = szVolAlias;
	m_nInputLines = nInputLines;
}

void CPincell::SetPitch (double dFlatF, double dZL)
// ---------------------------------------------------------------------------
// Function: sets the first line of pin input
// Input:    volume id of the pin, alias and total no. of lines in the pin input
// Output:   none
// ---------------------------------------------------------------------------
{
	m_dFlatF = dFlatF;
	m_dZL = dZL;
}

void CPincell::SetPitch (double dPX, double dPY, double dPZ)
// ---------------------------------------------------------------------------
// Function: sets the first line of pin input
// Input:    volume id of the pin, alias and total no. of lines in the pin input
// Output:   none
// ---------------------------------------------------------------------------
{
	m_dPX = dPX;
	m_dPY = dPY;
	m_dPZ = dPZ;
}

void CPincell::SetMatArray (int nMaterials)
// ---------------------------------------------------------------------------
// Function: sets the first line of pin input
// Input:    volume id of the pin, alias and total no. of lines in the pin input
// Output:   none
// ---------------------------------------------------------------------------
{
	m_nMaterials = nMaterials;
	m_szVMatName.SetSize(nMaterials);
	m_szVMatAlias.SetSize(nMaterials);
}

void CPincell::SetMat(CVector<std::string> szVMatName, CVector<std::string> szVMatAlias)
// ---------------------------------------------------------------------------
// Function: sets the first line of pin input
// Input:    volume id of the pin, alias and total no. of lines in the pin input
// Output:   none
// ---------------------------------------------------------------------------
{
	m_szVMatName = szVMatName;
	m_szVMatAlias = szVMatAlias;
}

void CPincell::SetNumCyl(const int nCyl)
// ---------------------------------------------------------------------------
// Function: sets the first line of pin input
// Input:    volume id of the pin, alias and total no. of lines in the pin input
// Output:   none
// ---------------------------------------------------------------------------
{
	m_nNumCyl = nCyl;
	m_VCyl.SetSize(nCyl);
}

void CPincell::GetNumCyl(int &nCyl)
// ---------------------------------------------------------------------------
// Function: sets the first line of pin input
// Input:    volume id of the pin, alias and total no. of lines in the pin input
// Output:   none
// ---------------------------------------------------------------------------
{
	nCyl = m_VCyl.GetSize();
}

void CPincell::SetCylSizes(int nCyl,int nRadii)
// ---------------------------------------------------------------------------
// Function: sets the first line of pin input
// Input:    volume id of the pin, alias and total no. of lines in the pin input
// Output:   none
// ---------------------------------------------------------------------------
{
	m_VCyl(nCyl).SetSizes(nRadii);
}

void CPincell::SetCylPos(int nCyl, CVector<double> dVCoor)
// ---------------------------------------------------------------------------
// Function: sets the first line of pin input
// Input:    volume id of the pin, alias and total no. of lines in the pin input
// Output:   none
// ---------------------------------------------------------------------------
{
	m_VCyl(nCyl).SetPos(dVCoor);
}


void CPincell::SetCylRadii(int nCyl, CVector<double> dVRadii)
// ---------------------------------------------------------------------------
// Function: sets the first line of pin input
// Input:    volume id of the pin, alias and total no. of lines in the pin input
// Output:   none
// ---------------------------------------------------------------------------
{
	m_VCyl(nCyl).SetRadii(dVRadii);
}



void CPincell::SetCylZPos(int nCyl, CVector<double> dVZCoor)
// ---------------------------------------------------------------------------
// Function: sets the first line of pin input
// Input:    volume id of the pin, alias and total no. of lines in the pin input
// Output:   none
// ---------------------------------------------------------------------------
{
	m_VCyl(nCyl).SetZPos(dVZCoor);
}

void CPincell::SetCylMat(int nCyl, CVector<std::string> szVMat)
// ---------------------------------------------------------------------------
// Function: sets the first line of pin input
// Input:    volume id of the pin, alias and total no. of lines in the pin input
// Output:   none
// ---------------------------------------------------------------------------
{
	m_VCyl(nCyl).SetMat(szVMat);
}

void CPincell::SetCellMatSize(int nSize)
// ---------------------------------------------------------------------------
// Function: sets the first line of pin input
// Input:    volume id of the pin, alias and total no. of lines in the pin input
// Output:   none
// ---------------------------------------------------------------------------
{
	m_dVZStart.SetSize(nSize);
	m_dVZEnd.SetSize(nSize);
	m_szVCellMat.SetSize(nSize);
}

void CPincell::SetCellMat(CVector<double> dZVStart, CVector<double> dVZEnd, CVector<std::string> szVCellMat) 
// ---------------------------------------------------------------------------
// Function: sets the first line of pin input
// Input:    volume id of the pin, alias and total no. of lines in the pin input
// Output:   none
// ---------------------------------------------------------------------------
{
	m_dVZStart = dZVStart;
	m_dVZEnd = dVZEnd;
	m_szVCellMat = szVCellMat;
}

//////


void CPincell::GetLineOne (std::string &szVolId, std::string &szVolAlias, int &nInputLines)
// ---------------------------------------------------------------------------
// Function: sets the first line of pin input
// Input:    volume id of the pin, alias and total no. of lines in the pin input
// Output:   none
// ---------------------------------------------------------------------------
{
	szVolId = m_szVolId;
	szVolAlias = m_szVolAlias;
	nInputLines = m_nInputLines;
}

void CPincell::GetPitch (double &dFlatF, double &dZL)
// ---------------------------------------------------------------------------
// Function: sets the first line of pin input
// Input:    volume id of the pin, alias and total no. of lines in the pin input
// Output:   none
// ---------------------------------------------------------------------------
{
	dFlatF = m_dFlatF;
	dZL = m_dZL;
}

void CPincell::GetPitch (double &dPX, double &dPY, double  &dPZ)
// ---------------------------------------------------------------------------
// Function: sets the first line of pin input
// Input:    volume id of the pin, alias and total no. of lines in the pin input
// Output:   none
// ---------------------------------------------------------------------------
{
	dPX = m_dPX;
	dPY = m_dPY;
	dPZ = m_dPZ;
}

void CPincell::GetCylSizes(int &nCyl,int &nRadii)
// ---------------------------------------------------------------------------
// Function: sets the first line of pin input
// Input:    volume id of the pin, alias and total no. of lines in the pin input
// Output:   none
// ---------------------------------------------------------------------------
{
	m_VCyl(nCyl).GetSizes(nRadii);
}

void CPincell::GetMat(CVector<std::string> &szVMatName, CVector<std::string> &szVMatAlias)
// ---------------------------------------------------------------------------
// Function: sets the first line of pin input
// Input:    volume id of the pin, alias and total no. of lines in the pin input
// Output:   none
// ---------------------------------------------------------------------------
{
	szVMatName = m_szVMatName;
	szVMatAlias = m_szVMatAlias;
}

void CPincell::GetCylPos(int &nCyl, CVector<double> &dVCoor)
// ---------------------------------------------------------------------------
// Function: sets the first line of pin input
// Input:    volume id of the pin, alias and total no. of lines in the pin input
// Output:   none
// ---------------------------------------------------------------------------
{
	m_VCyl(nCyl).GetPos(dVCoor);
}

void CPincell::GetCylZPos(int &nCyl, CVector<double> &dVCoor)
// ---------------------------------------------------------------------------
// Function: sets the first line of pin input
// Input:    volume id of the pin, alias and total no. of lines in the pin input
// Output:   none
// ---------------------------------------------------------------------------
{
	m_VCyl(nCyl).GetZPos(dVCoor);
}


void CPincell::GetCylRadii(int &nCyl, CVector<double> &dVRadii)
// ---------------------------------------------------------------------------
// Function: sets the first line of pin input
// Input:    volume id of the pin, alias and total no. of lines in the pin input
// Output:   none
// ---------------------------------------------------------------------------
{
	m_VCyl(nCyl).GetRadii(dVRadii);
}

void CPincell::GetCylMat(int &nCyl, CVector<std::string> &szVMat)
// ---------------------------------------------------------------------------
// Function: sets the first line of pin input
// Input:    volume id of the pin, alias and total no. of lines in the pin input
// Output:   none
// ---------------------------------------------------------------------------
{
	m_VCyl(nCyl).GetMat(szVMat);
}



void CPincell::GetCellMat(CVector<double> &dVZStart, CVector<double> &dVZEnd, CVector<std::string> &szVCellMat) 
// ---------------------------------------------------------------------------
// Function: sets the first line of pin input
// Input:    volume id of the pin, alias and total no. of lines in the pin input
// Output:   none
// ---------------------------------------------------------------------------
{
	dVZStart = m_dVZStart;
	dVZEnd = m_dVZEnd;
	szVCellMat = m_szVCellMat;
}

void CPincell::GetCellMatSize(int &nSize)
// ---------------------------------------------------------------------------
// Function: sets the first line of pin input
// Input:    volume id of the pin, alias and total no. of lines in the pin input
// Output:   none
// ---------------------------------------------------------------------------
{
   nSize = m_szVCellMat.GetSize();
}

void CPincell::GetMatArray (int &nMaterials)
// ---------------------------------------------------------------------------
// Function: sets the first line of pin input
// Input:    volume id of the pin, alias and total no. of lines in the pin input
// Output:   none
// ---------------------------------------------------------------------------
{
	nMaterials = m_nMaterials;
	
}
