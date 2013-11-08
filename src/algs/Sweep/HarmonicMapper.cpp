#include "HarmonicMapper.hpp"
#include <iostream>
#include <math.h>
#include <map>

namespace MeshKit {

HarmonicMapper::HarmonicMapper(MKCore* core, vector<Vertex> &v, vector<Face> &t, vector<Edge> &e, vector<set<int> > a)
{
  mk_core = core;
  vtx.insert(vtx.begin(), v.begin(), v.end());
  tri.insert(tri.begin(), t.begin(), t.end());
  edges.insert(edges.begin(), e.begin(), e.end());
  adj.insert(adj.begin(), a.begin(), a.end());
}

void HarmonicMapper::execute()
{
	_iterative_map(0.001);

}

void HarmonicMapper::getUV(vector<Vertex> &v)
{
	int count = -1;
	for (vector<Vertex>::iterator it = vtx.begin(); it != vtx.end(); it++){
		count++;
		if (it->onBoundary)  continue;
		v[count].uv[0] = it->uv[0];
		v[count].uv[1] = it->uv[1];
	}	

}

//matrix method
void HarmonicMapper::_map()
{
	int n_interior = 0, n_boundary = 0;
	for (vector<Vertex>::iterator it = vtx.begin(); it != vtx.end(); it++){
		if (it->onBoundary){
			it->id = n_boundary; 			
			n_boundary++;
		}
		else{
			it->id = n_interior;
			n_interior++;
		}
	}
	
	//using Armadillo to solve linear equations
	mat A = zeros<mat>(n_interior, n_interior);
	mat B = zeros<mat>(n_interior, n_boundary);
	//setting up the matrix A

	for (vector<Vertex>::iterator it = vtx.begin(); it != vtx.end(); it++){
		if (it->onBoundary) continue;
		int index_i = it->id;
		double sw = 0.0;
		for (set<int>::iterator it_adj = adj[it->index].begin(); it_adj != adj[it->index].end(); it_adj++){
			int index_vtx = -1;
			if (edges[*it_adj].connect[0]->index == it->index)
				index_vtx = edges[*it_adj].connect[1]->index;
			else
				index_vtx = edges[*it_adj].connect[0]->index;
			int index_j = vtx[index_vtx].id;
			double w = edges[*it_adj].e;
			if (vtx[index_vtx].onBoundary)
				B(index_i, index_j) = w;
			else
				A(index_i, index_j) = -1.0*w;
			
			sw += w;
		}
		A(index_i, index_i) = sw;
	}
	
	mat b = zeros<mat>(n_interior, 2);
	
	for (vector<Vertex>::iterator it = vtx.begin(); it != vtx.end(); it++){
		if (!it->onBoundary)  continue;
		b(it->id,0) = it->uv[0];
		b(it->id,1) = it->uv[1];
	}
	mat c = B * b;//n_interior * 2
    mat uv_coord = solve(A, c);
	
	int count = -1;
	for (vector<Vertex>::iterator it = vtx.begin(); it != vtx.end(); it++){
		if (it->onBoundary)  continue;
		count++;
		it->uv[0] = uv_coord(count, 0);
		it->uv[1] = uv_coord(count, 1);
	}	
	

}

//iterative method
void HarmonicMapper::_iterative_map(double epsilon)
{
	//move interior vertices to its center of neighbors
	//boundary nodes are set to (0,0)
	vector<int> interior;
	for (std::vector<Vertex>::iterator it = vtx.begin(); it != vtx.end(); it++)
		if (((*it).uv[0] == 0) && ((*it).uv[1] == 0))
			interior.push_back((*it).index);	

	while(true){
		double error = -1.0e+10;		
		
		
		for (std::vector<int>::iterator it = interior.begin(); it != interior.end(); it++){
			Vector2D uv;
			uv[0] = 0.0;
			uv[1] = 0.0;
			double weight = 0.0;
			for (set<int>::iterator it_e = adj[*it].begin(); it_e != adj[*it].end(); it_e++){
				int adj_v = -1;
				if (edges[*it_e].connect[0]->index == (*it))
					adj_v = edges[*it_e].connect[1]->index;
				else
					adj_v = edges[*it_e].connect[0]->index;
				uv[0] += edges[*it_e].e*vtx[adj_v].uv[0];
				uv[1] += edges[*it_e].e*vtx[adj_v].uv[1];
				weight += edges[*it_e].e;
			}
			Vector2D pre_uv;
			pre_uv[0] = vtx[*it].uv[0];
			pre_uv[1] = vtx[*it].uv[1];
			vtx[*it].uv[0] = uv[0]/weight;
			vtx[*it].uv[1] = uv[1]/weight;
			double v_err = sqrt(pow(pre_uv[0]-vtx[*it].uv[0],2)+pow(pre_uv[1]-vtx[*it].uv[1],2));
			error = (v_err > error)? v_err : error;
		}
		
		if (error < epsilon) break;
	}
}


HarmonicMapper::~HarmonicMapper()
{
  std::cout << "It is over now in smoothing" << endl;
}

}
