//-----------------------------------C++-------------------------------------//
// File: src/algs/IsoLaplace.hpp
// Monday Sep 26 10:50 2011
// Brief: IsoLaplace class definition: do the IsoLaplace smoothing for
//        the surface mesh 
//---------------------------------------------------------------------------//


#ifndef ISOLAPLACE_HPP
#define ISOLAPLACE_HPP

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
   * \class IsoLaplace
   * \brief do the IsoLaplace smoothing for surface mesh
   * 
   * 
   */
//===========================================================================//

class IsoLaplace
{	
public:
    IsoLaplace();
    ~IsoLaplace();


    //Execute function
    void Execute();

    //input the point coordinates, edges, triangles
    void SetupData(std::vector<std::set<int> > AdjElements, std::vector<std::vector<int> > Quads, std::vector<std::vector<double> > coords, std::vector<bool> isBnd, std::vector<double> w);

    //return the results
    void GetCoords(std::vector<std::vector<double> > &coords);


private://private member functions

    void UpdateWeight();
	

private://private member variable
    std::vector<std::vector<double> > coordinates;
    std::vector<std::set<int> > adjElements;
    std::vector<std::vector<int> > quads;
    std::vector<bool> isBoundary;
    std::vector<double> weight;

};

#endif
