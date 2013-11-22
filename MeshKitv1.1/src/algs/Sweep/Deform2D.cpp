#include "Deform2D.hpp"
#include <iostream>
#include <math.h>
#include <map>


extern "C" void openblas_set_num_threads(int);



namespace MeshKit {

Deform2D::Deform2D(vector<vector<double> > undeformed_cage_vertices, vector<vector<double> > deformed_cage_vertices)
{
    num_vertices = undeformed_cage_vertices.size();
    //pass the data
	un_cage_vertices.insert(un_cage_vertices.begin(), undeformed_cage_vertices.begin(), undeformed_cage_vertices.end());
	cage_vertices.insert(cage_vertices.begin(), deformed_cage_vertices.begin(), deformed_cage_vertices.end());
}

void Deform2D::SetupInteriorNodes(vector<vector<double> > undeformed_in_nodes)
{
    num_interior = undeformed_in_nodes.size();
    
    un_InnerNodes.insert(un_InnerNodes.begin(), undeformed_in_nodes.begin(), undeformed_in_nodes.end());
	InnerNodes.insert(InnerNodes.begin(), undeformed_in_nodes.begin(), undeformed_in_nodes.end());    
}
void Deform2D::GetInteriorNodes(vector<vector<double> > &final_locations)
{
   /*
   for (unsigned int i = 0; i < InnerNodes.size(); i++){
        for (unsigned int j = 0; j < InnerNodes[i].size(); j++){
            final_locations[i][j] = InnerNodes[i][j];
        }
   }
   */
   final_locations.clear();
   final_locations.insert(final_locations.begin(), InnerNodes.begin(), InnerNodes.end());
}

void Deform2D::Execute()
{
    const unsigned int vertex_size = un_cage_vertices.size();
    const unsigned int tmp_var = vertex_size + 3;
    
    mat A = zeros<mat>(vertex_size+3,vertex_size+3);
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
		A(vertex_size,j) = 1.0;
		A(vertex_size+1,j) = un_cage_vertices[j][0];
		A(vertex_size+2,j) = un_cage_vertices[j][1];
    }
	
    //assign G
    for (unsigned int i = 0; i < vertex_size; i++){
		for (unsigned int j = 0; j < vertex_size; j++){
			A(i,j) = TriHarmonicFun(un_cage_vertices[i], un_cage_vertices[j]);
		}
    }
    //done with the big matrix A
    //set up the displacement in x direction
    mat b=zeros<mat>(vertex_size+3,2);
    
    for (unsigned int i = 0; i < vertex_size; i++)
		b(i,0) = cage_vertices[i][0] - un_cage_vertices[i][0];
    //set up the displacement in y direction
    for (unsigned int i = 0; i < vertex_size; i++)
		b(i,1) = cage_vertices[i][1] - un_cage_vertices[i][1];

    std::cout << "Matrix size is " << vertex_size+3 << "*" << vertex_size+3 << "\n";
    openblas_set_num_threads(1);
    mat coeffs = solve(A,b);

    //compute the displacement in x and y direction for interior nodes
    for (unsigned int i = 0; i < num_interior; i++){
		InnerNodes[i][0] = un_InnerNodes[i][0] + InterpolatedFun(un_InnerNodes[i], coeffs.col(0));
		InnerNodes[i][1] = un_InnerNodes[i][1] + InterpolatedFun(un_InnerNodes[i], coeffs.col(1));
    }
    
}

double Deform2D::InterpolatedFun(vector<double> x, vec coeffs)
{
    double tmp = 0.0;
    double tmp_variable[3] = {1.0, x[0], x[1]};
    for (unsigned int i = 0; i < num_vertices; i++)
	    tmp += coeffs(i,0)*TriHarmonicFun(x, un_cage_vertices[i]);
    for (int i = 0; i < 3; i++)
	    tmp += coeffs(num_vertices+i,0)*tmp_variable[i];

    return tmp;
}

double Deform2D::TriHarmonicFun(vector<double> xi, vector<double> xj)
{
    double value = 0.0;
    double delta = sqrt(pow((xi[0] - xj[0]),2) + pow((xi[1] - xj[1]),2));
    value = delta*delta*log10(delta);
    return value;
}

Deform2D::~Deform2D()
{
	cout << "It is over now in computing 2D deformation" << endl;
}

}
