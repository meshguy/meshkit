
#include "meshkit/MKCore.hpp"
#include "Dijkstra.hpp"
#include <iostream>
#include <algorithm>
#include <math.h>
#include <map>

namespace MeshKit
{


//---------------------------------------------------------------------------//
// construction function for Dijkstra class
Dijkstra::Dijkstra(vector<vector<double> > t)
{
	dist.insert(dist.begin(), t.begin(), t.end());
}

void Dijkstra::getResults(vector<vector<double> > &d)
{
	d.resize(dist.size());
	order_dist.resize(dist.size());
	for (unsigned int i = 0; i < dist.size(); i++){
		algs(i, d[i]);
		order_dist[i].insert(order_dist[i].begin(), d[i].begin(), d[i].end());
	}
}

void Dijkstra::getTopMostSurf()
{

	for (unsigned int i = 0; i < order_dist.size(); i++){
		bool is_negative = true;
		for (unsigned int j = 0; j < order_dist[i].size(); j++)
			is_negative = is_negative && (order_dist[i][j] <= 0);
		if (is_negative)
			topmost_target_surf = i;
	}

	std::cout << "topmost target surf index = " << topmost_target_surf << std::endl;
}

void Dijkstra::getSurfList(vector<vector<vector<int> > > &l, int src_size, int tgt_size)
{
	std::cout << "====================\n";
	for (unsigned int i = 0; i < order_dist.size(); i++){

		for (unsigned int j = 0; j < order_dist[i].size(); j++)
			std::cout << order_dist[i][j] << "\t";
		std::cout << std::endl;
	}

	std::cout << "====================\n";

	std::vector<double> dist_array;
	getTopMostSurf();

	for (unsigned int i = 0; i < order_dist[topmost_target_surf].size(); i++)
		std::cout << order_dist[topmost_target_surf][i] << "\t";
	std::cout << std::endl;

	dist_array.insert(dist_array.begin(), order_dist[topmost_target_surf].begin(), order_dist[topmost_target_surf].end());
	std::vector<double> bak_dist_array;
	bak_dist_array.insert(bak_dist_array.begin(), dist_array.begin(), dist_array.end());
	//sort the distances
	std::sort(dist_array.begin(), dist_array.end());
	//get the index
	int layer_index = -1;
	double previous_dist = 1.0e10;

	std::cout << "sorted array\n";
	for (unsigned int i = 0; i < dist_array.size(); i++)
		std::cout << dist_array[i] << "\t";
	std::cout << std::endl;


	for (int i = dist_array.size() - 1; i >= 0; i--){
		std::vector<double>::iterator it;
		it = std::find(bak_dist_array.begin(), bak_dist_array.end(), dist_array[i]);
		int index = std::distance(bak_dist_array.begin(), it);
		surf_list.push_back(index);
		//l.push_back(index);
		if (i==( (int)dist_array.size()-1)){
			layer_index++;
			l.resize(layer_index+1);
			addSurfToList(layer_index, index, src_size, l[layer_index]);
		}
		else{
			if (fabs(dist_array[i]-previous_dist) < 1.0e-5){//on the same layer
				addSurfToList(layer_index, index, src_size, l[layer_index]);
			}
			else{
				layer_index++;
				l.resize(layer_index+1);
				addSurfToList(layer_index, index, src_size, l[layer_index]);
			}

		}
		previous_dist = dist_array[i];
		bak_dist_array[index] = 1.0e10;
	}

	//
	for (unsigned int i = 0; i < l.size(); i++){
		for (unsigned int j = 0; j < l[i].size(); j++){
			std::cout << "layer = " << i << "\tsurf index = " << j << "\tsurf_id = " << l[i][j][0] << std::endl;
		}
	}
}

void Dijkstra::addSurfToList(int layer_index, int value, int src_size, vector<vector<int> > &datalist){
	int surf_index = datalist.size();
	datalist.resize(surf_index+1);
	if (value < src_size){
		datalist[surf_index].push_back(value);
		datalist[surf_index].push_back(0);
	}
	else{
		datalist[surf_index].push_back(value-src_size);
		datalist[surf_index].push_back(1);
	}
}

void Dijkstra::algs(int s, vector<double> &f_list)
{
	//Initializaton: set every distance to INFINITY until we discover a path
	vector<double> d(dist.size(), 1.0e10);
	//sptSet[i] will true if vertex i is included in the shortest path tree or
	//shortest distance from s to i is finalized.
	vector<bool> sptSet(dist.size(), false);
	//the distance from the source to the source is defined to be zero.
	d[s] = 0.0;
	/*
	This loop corresponds to sending out the explorers walking the paths,
	where the step of picking "the vertex, v, with the shortest path to s"
	corresponds to an explorer arriving at an unexplored vertex.
	*/
	int V = (int)dist.size();
	// Find shortest path for all vertices
	for (int count = 0; count < V-1; count++){
		//pick the minimum distance vertex from the set of vertices not yet
		//processed. u is always equal to s in the first iteration
		int u = minDistance(d, sptSet);
		//Mark the picked vertex as processed.
		sptSet[u] = true;
		//Update d value of the adjacent vertices of the picked vertex
		for (int v = 0; v < V; v++){
			//update d[v] only if is not in sptSet, there is an edge from u to
			//v, and total weight of apth from s to v through u is smaller than
			//current value of d[v]
			if (!sptSet[v] && dist[u][v] && d[u] != 1.0e10 && d[u]+dist[u][v]<d[v])
				d[v] = d[u]+dist[u][v];
		}
	}

	//pass results to the list
	for (unsigned int i = 0; i < d.size(); i++)
		f_list.push_back(d[i]);

}
//A utility function to find the vertex with minimum distance value, from the
//set of vertices not yet included in shortest path tree
int Dijkstra::minDistance(vector<double> d, vector<bool> sptSet)
{
	//Initialize min value
	double minvalue = 1.0e10;
	int min_index = 0;
	for (int v = 0; v < (int) dist.size(); v++){
		if (sptSet[v] == false && d[v] <= minvalue){
			minvalue = d[v];
			min_index = v;
		}
	}
	return min_index;
}

//---------------------------------------------------------------------------//
// deconstruction function for Dijkstra class
Dijkstra::~Dijkstra()
{
	std::cout << "Dijkstra algorithm is over\n";
}




}

