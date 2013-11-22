#include "GreenCoordinates3D.hpp"
#include <iostream>
#include <math.h>
#include <map>




namespace MeshKit {

GreenCoordinates3D::GreenCoordinates3D(MKCore* core, vector< vector<double> > iNodes)
{
	mk_core = core;
	InteriorNodes.resize(iNodes.size());
	for (unsigned int i = 0; i < iNodes.size(); i++){
		InteriorNodes[i].resize(3);
		for (int j = 0; j < 3; j++)
			InteriorNodes[i][j] = iNodes[i][j];
	}
	
}

//Setup the cages for surface mesh
void GreenCoordinates3D::SetupCages(vector<vector<double> > cageNodes, vector< vector<int> > cageFaces, vector<vector<double> > normal)
{
	int num_CageNodes = cageNodes.size(), num_CageFaces = cageFaces.size();	
	CageNodes.resize(num_CageNodes);
	CageFaces.resize(num_CageFaces);
	
	for (int i = 0; i < num_CageNodes; i++){
		int size_j = cageNodes[i].size();
		CageNodes.resize(size_j);
		for (int j = 0; j < size_j; j++)
			CageNodes[i][j] = cageNodes[i][j];
	}

	for (int i = 0; i < num_CageFaces; i++){
		int size_j = CageFaces[i].size();
		CageFaces.resize(size_j);
		for (int j = 0; j < size_j; j++)
			CageFaces[i][j] = cageFaces[i][j];
	}

	Normals.resize(normal.size());
	for (unsigned int i = 0; i < normal.size(); i++){
		unsigned int size_j = normal[i].size();
		Normals[i].resize(size_j);
		for (unsigned int j = 0; j < size_j; j++)		
			Normals[i][j] = normal[i][j];
	}
}

//calculate the green coordinates
void GreenCoordinates3D::Execute()
{
	//Initialization
	double error = 1.0e-5;
	weight_t.resize(InteriorNodes.size());
	weight_v.resize(InteriorNodes.size());
	for (unsigned int i = 0; i < weight_t.size(); i++){
		weight_t[i].resize(CageFaces.size());		
		for (unsigned int j = 0; j < CageFaces.size(); j++)
			weight_t[i][j] = 0.0;
	}

	for (unsigned int i = 0; i < weight_v.size(); i++){
		weight_v[i].resize(CageNodes.size());
		for (unsigned int j = 0; j < CageNodes.size(); j++)		
			weight_v[i][j] = 0.0;
	}

	//Main Loop--Calculate the weight for cage nodes and cage face normal
	//From the appendix of the paper: Green Coordinates, Yaron Lipman, David Levin, Daniel Cohen-Or
	for (unsigned int i = 0; i < InteriorNodes.size(); i++){
		for (unsigned int j = 0; j < CageFaces.size(); j++){
			vector<vector<double> > node(3, vector<double>(3));
			for (int k = 0; k < 3; k++){
				node[0][k] = CageNodes[CageFaces[j][0]][k] - InteriorNodes[i][k];
				node[1][k] = CageNodes[CageFaces[j][1]][k] - InteriorNodes[i][k];
				node[2][k] = CageNodes[CageFaces[j][2]][k] - InteriorNodes[i][k];
			}
			double mag = node[0][0]*Normals[j][0] + node[1][1]*Normals[j][1] + node[2][2]*Normals[j][2];
			vector<double> p(3);
			for (int k = 0; k < 3; k++)
				p[k] = mag*Normals[j][k];

			vector<double> s(3),I(3),II(3);
			vector<vector<double> > q(3, vector<double>(3)), N(3, vector<double>(3));
			for (int k = 0; k < 3; k++){
				vector<double> ivector(3), jvector(3);
				for (int l = 0; l < 3; l++){
					ivector[l] = node[k][l] - p[l];
					jvector[l] = node[(k+1)%3][l]- p[l];
				}
				vector<double> crossproduct(3);
				crossproduct[0] = ivector[1]*jvector[2] - ivector[2]*jvector[1];
				crossproduct[1] = -1.0*(ivector[0]*jvector[2]-ivector[2]*jvector[0]);
				crossproduct[2] = ivector[0]*jvector[1] - ivector[1]*jvector[0];
				double dotproduct = crossproduct[0]*Normals[j][0]+crossproduct[1]*Normals[j][1]+crossproduct[2]*Normals[j][2];
				s[k] = dotproduct>0? 1.0:-1.0;
				vector<double> tmp(3);
				tmp[0] = 0.0; tmp[1] = 0.0; tmp[2] = 0.0;
				I[k] = GCTriInt(p, node[k], node[(k+1)%3], tmp);
				II[k] = GCTriInt(tmp, node[k], node[(k+1)%3], tmp);	

				q[k][0] = node[(k+1)%3][1]*node[k][2]-node[(k+1)%3][2]*node[k][1];
				q[k][1] = -1.0*(node[(k+1)%3][0]*node[k][2]-node[(k+1)%3][2]*node[k][0]);
				q[k][2] = node[(k+1)%3][0]*node[k][1]-node[(k+1)%3][1]*node[k][0];

				mag = sqrt(pow(q[k][0],2) + pow(q[k][1],2) + pow(q[k][2], 2));
				N[k][0] = q[k][0]/mag;
				N[k][1] = q[k][1]/mag;
				N[k][2] = q[k][2]/mag;						
			}

			double I_scalar = -1.0*fabs(s[0]*I[0]+s[1]*I[1]+s[2]*I[2]);
			weight_t[i][j] = -1.0*I_scalar;
			
			vector<double> w(3);
			w[0] = I_scalar*Normals[j][0] + N[0][0]*I_scalar*I[0] + N[1][0]*I_scalar*I[1] + N[2][0]*I_scalar*I[2];
			w[1] = I_scalar*Normals[j][1] + N[0][1]*I_scalar*I[0] + N[1][1]*I_scalar*I[1] + N[2][1]*I_scalar*I[2];
			w[2] = I_scalar*Normals[j][2] + N[0][2]*I_scalar*I[0] + N[1][2]*I_scalar*I[1] + N[2][2]*I_scalar*I[2];

			double w_mag = sqrt(pow(w[0],2)+pow(w[1],2)+pow(w[2],2));
			if (w_mag > error)//this guarantees the local deformation property
				for (int k = 0; k < 3; k++)
					weight_v[i][CageFaces[j][k]] += (N[(k+1)%3][0]*w[0]+N[(k+1)%3][1]*w[1]+N[(k+1)%3][2]*w[2])/(N[(k+1)%3][0]*node[k][0]+N[(k+1)%3][1]*node[k][1]+N[(k+1)%3][2]*node[k][2]);
				
			
		}

	}
}



double GreenCoordinates3D::GCTriInt(vector<double> p, vector<double> v1, vector<double> v2, vector<double> iNodes)
{
	//calculate the angle
	double dotproduct_21p1 = (v2[0]-v1[0])*(p[0]-v1[0]) + (v2[1]-v1[1])*(p[1]-v1[1]) + (v2[2]-v1[2])*(p[2]-v1[2]);
	double dotproduct_1p2p = (v1[0]-p[0])*(v2[0]-p[0]) + (v1[1]-p[1])*(v2[1]-p[1]) + (v1[2]-p[2])*(v2[2]-p[2]);
	double mag_21 = sqrt((v2[0]-v1[0])*(v2[0]-v1[0])+(v2[1]-v1[1])*(v2[1]-v1[1])+(v2[2]-v1[2])*(v2[2]-v1[2]));
	double mag_p1 = sqrt((p[0]-v1[0])*(p[0]-v1[0])+(p[1]-v1[1])*(p[1]-v1[1])+(p[2]-v1[2])*(p[2]-v1[2]));
	double mag_p2 = sqrt((p[0]-v2[0])*(p[0]-v2[0])+(p[1]-v2[1])*(p[1]-v2[1])+(p[2]-v2[2])*(p[2]-v2[2]));

	double alpha = acos(dotproduct_21p1/(mag_21*mag_p1));
	double belta = acos(dotproduct_1p2p/(mag_p1*mag_p2));

	double lambda = pow(mag_p1,2)*pow(sin(alpha), 2);
	double c = pow(p[0]-iNodes[0],2) + pow(p[1]-iNodes[1],2) + pow(p[2]-iNodes[2],3);

	double pi = atan(1.0)*4.0;
	double S = sin(pi-alpha), C = cos(pi-alpha);
	double theta_pi_alpha = -1.0*(S >0? 1.0:-1.0)*0.5*(2.0*sqrt(C)*atan(sqrt(c)*C/sqrt(lambda+S*S*c)) + sqrt(lambda)*log(2*sqrt(lambda)*pow(S,2)*(1-2*c*C/(c*(1+C)+lambda+sqrt(pow(lambda,2)+lambda*c*pow(S,2))))/pow(1-C,2)));
	
	S = sin(pi-alpha-belta);
	C = cos(pi-alpha-belta);
	double theta_pi_alpha_belta = -1.0*(S >0? 1.0:-1.0)*0.5*(2.0*sqrt(C)*atan(sqrt(c)*C/sqrt(lambda+S*S*c)) + sqrt(lambda)*log(2*sqrt(lambda)*pow(S,2)*(1-2*c*C/(c*(1+C)+lambda+sqrt(pow(lambda,2)+lambda*c*pow(S,2))))/pow(1-C,2)));

	return -1.0*fabs(theta_pi_alpha-theta_pi_alpha_belta-sqrt(c)*belta)/(4*pi);	
}


//input new cage and calculate deformed interior points based on the existing weights for cage vertices and cage face normals
void GreenCoordinates3D::GetDeformedVertices(vector< vector<double> > nodes, vector< vector<double> > norm, vector<vector<double> > &ReturnNodes)
{
	ReturnNodes.resize(InteriorNodes.size());
	//loop over the interior nodes
	for (unsigned int i = 0; i < InteriorNodes.size(); i++){
		ReturnNodes[i].resize(3);
		for (int j = 0; j < 3; j++)
			ReturnNodes[i][j] = 0.0;
		//consider the effect of cage nodes		
		for (unsigned int j = 0; j < nodes.size(); j++)
			for (int k = 0; k < 3; k++)
				ReturnNodes[i][k] += weight_v[i][j]*nodes[j][k];
		
		//consider the effect of cage face normals
		for (unsigned int j = 0; j < norm.size(); j++)
			for (int k = 0; k < 3; k++)
				ReturnNodes[i][k] += weight_t[i][j]*norm[j][k];
	}
}

GreenCoordinates3D::~GreenCoordinates3D()
{
	cout << "It is over now in calculating the green coordinates for " << endl;
}

}
