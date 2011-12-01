//-----------------------------------C++-------------------------------------//
// File: src/algs/EquipotentialSmooth.hpp
// Monday Sep 26 10:50 2011
// Brief: EquipotentialSmooth class definition: do the EquipotentialSmooth smoothing for
//        the surface mesh 
//---------------------------------------------------------------------------//


#ifndef EQUIPOTENTIALSMOOTH_HPP
#define EQUIPOTENTIALSMOOTH_HPP

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string>
#include <iostream>
#include <fstream>
#include <string.h>
#include <vector>
#include <set>
#include <list>
#include <math.h>
#include <map>


using namespace std;


//===========================================================================//
  /*!
   * \class EquipotentialSmooth
   * \brief do the EquipotentialSmooth smoothing for surface mesh (Winslow)
   * This class is only for smoothing structured quad mesh
   * 
   */
//===========================================================================//

class EquipotentialSmooth
{	
public:
    EquipotentialSmooth();
    ~EquipotentialSmooth();


    //Execute function
    void Execute();

    //input the point coordinates, edges, triangles
    void SetupData(std::vector<std::set<int> > AdjElements, std::vector<std::vector<int> > Quads, std::vector<std::vector<double> > coords, std::vector<bool> isBnd, std::vector<std::vector<int> > conn);

    //return the results
    void GetCoords(std::vector<std::vector<double> > &coords);


private://private member functions

	

private://private member variable
    std::vector<std::vector<double> > coordinates;
    std::vector<std::set<int> > adjElements;
    std::vector<std::vector<int> > quads;
    std::vector<bool> isBoundary;
    std::vector<std::vector<int> > connect;

};

#endif
