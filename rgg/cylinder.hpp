/*********************************************
Dec,09
Reactor Geometry Generator
Argonne National Laboratory

Pin cylinder class definition.
*********************************************/
#ifndef __RGG_CYLINDER_H__
#define __RGG_CYLINDER_H__
#include <iostream>
#include "vectortemplate.hpp"


class CCylinder
{
public:
	CCylinder ();             // ctor
	CCylinder (const CCylinder&); // copy ctor
	~CCylinder ();            // dtor

	// accessor functions
	void GetSizes(int &nRadii);
	void GetPos(CVector<double>&);
	void GetRadii(CVector<double>&);
	void GetMat(CVector<std::string>&);
	void GetZPos(CVector<double>&);

	// modifier functions
	void SetSizes(int nRadii);
	void SetPos(CVector<double>);
	void SetRadii(CVector<double>);
	void SetMat(CVector<std::string>);
	void SetZPos(CVector<double>);

private:

	// pin related input
	//line one
	int m_nRadii;
	CVector<double> m_dVXYPos;
	CVector<double> m_dVZPos;
	CVector<std::string> m_szVMat;
	CVector<double> m_dVRadii;

	CVector<CCylinder> m_Cyl;
};

#endif	
