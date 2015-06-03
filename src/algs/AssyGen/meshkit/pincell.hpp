/*********************************************
Dec,09
Reactor Geometry Generator
Argonne National Laboratory

CPincell class definition.
*********************************************/
#ifndef __RGG_PINCELL_H__
#define __RGG_PINCELL_H__
#include <iostream>
#include "meshkit/vectortemplate.hpp"
#include "cylinder.hpp"

class CPincell
{
public:
  CPincell ();             // ctor
  CPincell (const CPincell&); // copy ctor
  ~CPincell ();            // dtor

  // accessor functions
  void GetNumCyl(int &nCyl);
  void GetLineOne (std::string &szVolId, std::string &szVolAlias, int &nInputLines);
  void GetIntersectFlag(int &nIFlag);
  void GetPitch (double &dFlatF, double &dZL);
  void GetPitch (double &dPX, double &dPY, double &dPZ);
  void GetMatArray(int &nMaterials);
  void GetMat(CVector<std::string> &szVMatName, CVector<std::string> &szVMatAlias);
  void GetNumCyl(const int &nCyl);
  void GetCylSizes(int &nCyl,int &nRadii);
  void GetCylRadii(int &nCyl, CVector<double>&);
  void GetCylZPos(int &nCyl, CVector<double>&);
  void GetCylPos(int &nCyl, CVector<double>&);
  void GetCylMat(int &nCyl, CVector<std::string>&);
  void GetCellMatSize(int &nSize);
  void GetCellMat(CVector<double> &dZStart, CVector<double> &dZEnd, CVector<std::string> &szVCellMat);


  // modifier functions

  void SetLineOne (std::string szVolId, std::string szVolAlias, int nInputLines);
  void SetIntersectFlag(int nIFlag);
  void SetPitch (double dFlatF, double dZL);
  void SetPitch (double dPX, double dPY, double dPZ);
  void SetMatArray(int nMaterials);
  void SetMat(CVector<std::string> szVMatName, CVector<std::string> szVMatAlias);
  // cylinder
  void SetNumCyl(const int nCyl);
  void SetCylSizes(int nCyl,int nRadii);
  void SetCylZPos(int nCyl, CVector<double>);
  void SetCylRadii(int nCyl, CVector<double>);
  void SetCylPos(int nCyl, CVector<double>);
  void SetCylMat(int nCyl, CVector<std::string>);
  //cell
  void SetCellMatSize(int nSize);
  void SetCellMat(CVector<double> dZStart, CVector<double> dZEnd, CVector<std::string> szVCellMat);
  void SetCellType(int nCyl, int nType);

private:

  // pin related input
  //line one
  std::string m_szVolId, m_szVolAlias;
  int m_nInputLines, m_nIFlag;
  // pitch
  double m_dPX; double m_dPY; double m_dPZ;
  double m_dFlatF; double m_dZL;
  // material
  int m_nMaterials;
  CVector<std::string> m_szVMatName, m_szVMatAlias;
  //cylinder
  CVector<CCylinder> m_VCyl;
  int m_nNumCyl;
  //cellmaterial
  CVector<double>		m_dVZStart;
  CVector<double> m_dVZEnd;
  CVector<std::string> m_szVCellMat;
};

#endif	
