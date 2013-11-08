#ifndef __HARMONICMAPPER_HPP
#define __HARMONICMAPPER_HPP

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
#include "meshkit/iGeom.hpp"
#include "meshkit/iRel.hpp"
#include "meshkit/MeshOp.hpp"
#include "meshkit/ModelEnt.hpp"
#include "meshkit/SimpleArray.hpp"
#include "Global.hpp"
#include <armadillo>
#include <set>
#include <vector>
#include <map>


using namespace std;
using namespace arma;

namespace MeshKit {

class HarmonicMapper
{	
public:
	//public function
	HarmonicMapper(MKCore* core, vector<Vertex> &v, vector<Face> &t, vector<Edge> &e, vector<set<int> > a);
	~HarmonicMapper();
	void execute();
	void getUV(vector<Vertex> &v);
	
private:
	void _iterative_map(double epsilon);
	void _map();
	
	
private:
	//member variables
	MKCore* mk_core;
	vector<set<int> > adj;
	vector<Vertex> vtx;
	vector<Edge> edges;
	vector<Face> tri;
	
};

}
#endif

