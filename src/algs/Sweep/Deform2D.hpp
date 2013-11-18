#ifndef __DEFORM2D_HPP
#define __DEFORM2D_HPP

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
#include "Global.hpp"
#include <set>
#include <vector>
#include <map>
//#include <Eigen/Core>
//#include <Eigen/Dense>
#include <armadillo>

//using namespace Eigen;
using namespace arma;

using namespace std;
//using namespace Mesquite;

namespace MeshKit {

class Deform2D
{	
public:
	//public function
	Deform2D(vector<vector<double> > undeformed_cage_vertices, vector<vector<double> > deformed_cage_vertices);
	void SetupInteriorNodes(vector<vector<double> > undeformed_in_nodes);
	void GetInteriorNodes(vector<vector<double> > &final_locations);
    	~Deform2D();
    	void Execute();

private:
	double TriHarmonicFun(vector<double> xi, vector<double> xj);
        double InterpolatedFun(vector<double> x, vec coeffs);
private:
    vector<vector<double> > un_cage_vertices;
    vector<vector<double> > cage_vertices;
    vector<vector<double> > un_InnerNodes;
    vector<vector<double> > InnerNodes;
    unsigned int num_vertices;
    unsigned int num_interior;
	
};

}
#endif

