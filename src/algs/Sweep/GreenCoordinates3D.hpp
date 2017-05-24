#ifndef __GREENCOORDINATES3D_HPP
#define __GREENCOORDINATES3D_HPP

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string>
#include <iostream>
#include <fstream>
#include <string.h>
#include <limits.h>

#include "meshkit/MKCore.hpp"
#include "meshkit/iMesh.hpp"
#include "meshkit/MeshOp.hpp"
#include "meshkit/ModelEnt.hpp"
#include <set>
#include <vector>
#include <map>

//Green Coordinates cage-based deformation
//The coordinates respect both the cage vertices position and face orientation
//It has shape-preserving property and closed-form.



using namespace std;

namespace MeshKit {

class GreenCoordinates3D
{	
public:
	//public function
	GreenCoordinates3D(MKCore* core, vector< vector<double> > iNodes);
	void SetupCages(vector<vector<double> > cageNodes, vector< vector<int> > cageFaces, vector<vector<double> > normal);
	~GreenCoordinates3D();
	void Execute();
	void GetDeformedVertices(vector< vector<double> > nodes, vector< vector<double> > norm, vector<vector<double> > &ReturnNodes);

private:
	double GCTriInt(vector<double> p, vector<double> v1, vector<double> v2, vector<double> iNodes);
	
private:
	//member variables
	MKCore* mk_core;
    vector< vector<double> > CageNodes;
	vector< vector<int> > CageFaces;
	vector< vector<double> > InteriorNodes;
	vector<vector<double> > Normals;
	vector<vector<double> > weight_v, weight_t;
	
};

}
#endif

