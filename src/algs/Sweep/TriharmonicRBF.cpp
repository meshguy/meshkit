#include "TriharmonicRBF.hpp"
#include <iostream>
#include <math.h>
#include <map>


extern "C" void openblas_set_num_threads(int);



namespace MeshKit {

TriharmonicRBF::TriharmonicRBF(vector<vector<double> > undeformed_cage_vertices, vector<vector<double> > deformed_cage_vertices, vector<vector<int> > t_cage_faces)
{
    num_vertices = undeformed_cage_vertices.size();
    //pass the data
    un_cage_vertices.resize(num_vertices);
    cage_vertices.resize(num_vertices);
    for (unsigned int i = 0; i < num_vertices; i++){
        un_cage_vertices[i].resize(3);
        cage_vertices[i].resize(3);
        for (int j = 0; j < 3; j++){
            un_cage_vertices[i][j] = undeformed_cage_vertices[i][j];
            cage_vertices[i][j] = deformed_cage_vertices[i][j];
        }
    }
    
}

void TriharmonicRBF::SetupInteriorNodes(vector<vector<double> > undeformed_in_nodes)
{
    num_interior = undeformed_in_nodes.size();
    un_InnerNodes.resize(num_interior);
    InnerNodes.resize(num_interior);
    for (unsigned int i = 0; i < undeformed_in_nodes.size(); i++){
        un_InnerNodes[i].resize(3);
        InnerNodes[i].resize(3);
        for (int j = 0; j < 3; j++){
            un_InnerNodes[i][j]=undeformed_in_nodes[i][j];
        }
    }
    
}
void TriharmonicRBF::GetInteriorNodes(vector<vector<double> > &final_locations)
{
   for (unsigned int i = 0; i < InnerNodes.size(); i++){
        for (unsigned int j = 0; j < InnerNodes[i].size(); j++){
            final_locations[i][j] = InnerNodes[i][j];
        }
   }
}

void TriharmonicRBF::Execute()
{
    const unsigned int vertex_size = un_cage_vertices.size();
    const unsigned int tmp_var = vertex_size + 4;
    
    mat A = zeros<mat>(vertex_size+4,vertex_size+4);
    //mat::fixed<tmp_var,tmp_var> A;
    //A.zeros();
    /**************************/
    /*	    |	G	H |	 */
    /*	A = |		  |	 */			       
    /*	    |	H'	0 |      */			       
    /*****************************/
    //Assign H and H'
    for (unsigned int j = 0; j < vertex_size; j++){
	A(j,vertex_size) = 1.0;
        A(j,vertex_size+1) = un_cage_vertices[j][0];
	A(j,vertex_size+2) = un_cage_vertices[j][1];
	A(j,vertex_size+3) = un_cage_vertices[j][2];
	A(vertex_size,j) = 1.0;
	A(vertex_size+1,j) = un_cage_vertices[j][0];
	A(vertex_size+2,j) = un_cage_vertices[j][1];
	A(vertex_size+3,j) = un_cage_vertices[j][2];
    }
	
    //assign G
    for (unsigned int i = 0; i < vertex_size; i++){
	for (unsigned int j = 0; j < vertex_size; j++){
	    A(i,j) = TriHarmonicFun(un_cage_vertices[i], un_cage_vertices[j]);
	}
    }
    //done with the big matrix A
    //set up the displacement in x direction
    mat b=zeros<mat>(vertex_size+4,3);
    //mat::fixed<tmp_var,3> b;
    //b.zeros();
    for (unsigned int i = 0; i < vertex_size; i++)
	b(i,0) = cage_vertices[i][0] - un_cage_vertices[i][0];
    //solve the linear equations: least-squares solving with a SVD decomposition.
    //Eigen provides one as the JacobiSVD class and its solve() is doing least-square solving
    //MatrixXf coeffs_x = A.jacobiSvd(ComputeThinU | ComputeThinV).solve(b);
    for (unsigned int i = 0; i < vertex_size; i++)
	b(i,1) = cage_vertices[i][1] - un_cage_vertices[i][1];
    //solve the linear equations
    //MatrixXf coeffs_y = A.jacobiSvd(ComputeThinU | ComputeThinV).solve(b);

    //set up the displacement in z direction
    for (unsigned int i = 0; i < vertex_size; i++)
	b(i,2) = cage_vertices[i][2] - un_cage_vertices[i][2];
    //solve the linear equations
    //MatrixXf coeffs_z = A.jacobiSvd(ComputeThinU | ComputeThinV).solve(b);
    //MatrixXf coeffs_z = A.householderQr().solve(b);
    std::cout << "Matrix size is " << vertex_size+4 << "*" << vertex_size+4 << "\n";
    openblas_set_num_threads(4);
    mat coeffs = solve(A,b);
    //MatrixXf coeffs = A.jacobiSvd(ComputeThinU | ComputeThinV).solve(b);
    

    //compute the displacement in x direction for interior nodes
    for (unsigned int i = 0; i < num_interior; i++){
	InnerNodes[i][0] = un_InnerNodes[i][0] + InterpolatedFun(un_InnerNodes[i], coeffs.col(0));
	InnerNodes[i][1] = un_InnerNodes[i][1] + InterpolatedFun(un_InnerNodes[i], coeffs.col(1));
	InnerNodes[i][2] = un_InnerNodes[i][2] + InterpolatedFun(un_InnerNodes[i], coeffs.col(2));
    }
    
}

double TriharmonicRBF::InterpolatedFun(vector<double> x, vec coeffs)
{
    double tmp = 0.0;
    double tmp_variable[4] = {1.0, x[0], x[1], x[2]};
    for (unsigned int i = 0; i < num_vertices; i++)
	    tmp += coeffs(i,0)*TriHarmonicFun(x, un_cage_vertices[i]);
    for (int i = 0; i < 4; i++)
	    tmp += coeffs(num_vertices+i,0)*tmp_variable[i];

    return tmp;
}

double TriharmonicRBF::TriHarmonicFun(vector<double> xi, vector<double> xj)
{
    double value = 0.0;
    double delta = pow((xi[0] - xj[0]),2) + pow((xi[1] - xj[1]),2) + pow((xi[2] - xj[2]),2);
    //value = sqrt(log10(delta+1));
    //value = pow(delta, 0.5);
    //double a = 1.0e-3;
    //value = 1.0/sqrt(a*a+delta);
    // value = sqrt(a*a+delta);
    //value = exp(-1.0*delta);
    value = pow(delta, 1.5);
    return value;
}

TriharmonicRBF::~TriharmonicRBF()
{
	cout << "It is over now in calculating the TriHarmonic" << endl;
}

}
