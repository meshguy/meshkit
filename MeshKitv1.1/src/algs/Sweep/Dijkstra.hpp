//-----------------------------------C++-------------------------------------//
// File: src/algs/Dijkstra.hpp
// Wednesday February 11 10:50 2011
// Brief: Dijkstra typically only provides the lengths between all pairs of 
// vertices
//        
//---------------------------------------------------------------------------//


#ifndef MESHKIT_DIJKSTRAIMPRINTING_HPP
#define MESHKIT_DIJKSTRAIMPRINTING_HPP

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string>
#include <iostream>
#include <fstream>
#include <string.h>
#include <limits.h>

#include <iGeom.h>
#include "meshkit/MKCore.hpp"
#include "meshkit/ModelEnt.hpp"
#include "meshkit/SizingFunction.hpp"
#include "moab/ReadUtilIface.hpp"
#include "Global.hpp"
#include <iMesh.h>
#include <iGeom.h>
#include <set>
#include <iRel.h>
#include <vector>
#include <set>
#include <list>

#include "meshkit/MeshScheme.hpp"

using namespace std;

namespace MeshKit
{
//===========================================================================//
  /*!
   * \class Dijkstra
   * \brief build lengths between each pair of vertices
   * 
   * Dijkstra
   */
//===========================================================================//

class Dijkstra
{	
public:

	Dijkstra(vector<vector<double> > t);
	~Dijkstra();
	void getResults(vector<vector<double> > &d);
	
	void getSurfList(vector<vector<vector<int> > > &l, int src_size, int tgt_size);
	
	

private:
	void algs(int s, vector<double> &f_list);
	int minDistance(vector<double> d, vector<bool> sptSet);
	void addSurfToList(int layer_index, int value, int src_size, vector<vector<int> > &datalist);
	void getTopMostSurf();

private://private member variable
	vector<vector<double> > dist;
	int topmost_target_surf;
	vector<vector<double> > order_dist;
	vector<int> surf_list;
};

}

#endif


