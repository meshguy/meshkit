#include <string>
#include <iostream>
#include <cassert>

#include "Tri2Quad.h"
#include "QuadCleanUp.h"

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/max_cardinality_matching.hpp>

using namespace boost;
using namespace Jaal;

typedef adjacency_list<vecS, vecS, undirectedS> Graph; 

int nx = 100;
int ny = 100;

////////////////////////////////////////////////////////////////////////

void edmond_graph_matching( const string &fname)
{
   Jaal::Mesh *inmesh  = new Mesh;

// inmesh = Jaal::struct_tri_grid( nx, ny);
  inmesh->readFromFile( fname );

  QuadCleanUp qClean(inmesh);

  inmesh->get_topological_statistics();

  inmesh->saveAs( "inmesh");
  cout << "Input Mesh: " << endl;
  cout << "      #Nodes     : " << inmesh->getSize(0) << endl;
  cout << "      #Triangles : " << inmesh->getSize(2) << endl;
  cout << "Building Dual Graph ... " << endl;

  ticks start_tick = getticks();
  Jaal::DualGraph dgraph;
  dgraph.build( inmesh );
  dgraph.setAdjTable(1,0);

  int numnodes = dgraph.getSize(0);
  int numedges = dgraph.getSize(1);

  Graph graph(numnodes);
  for( int i = 0; i < numedges; i++) {
       PEdge edge = dgraph.getEdge(i);
       int v0 = edge->getNodeAt(0)->getID();
       int v1 = edge->getNodeAt(1)->getID();
       add_edge(v0,v1,graph);
  }

  std::vector<graph_traits<Graph>::vertex_descriptor> mate(numnodes);

  cout << "Edmonds Graph Matching ... " << endl;
  bool success = checked_edmonds_maximum_cardinality_matching(graph, &mate[0]);
  assert(success);

  cout << "Boost:: maximum matching size " << matching_size(graph, &mate[0]) << endl;

  vector<FacePair>  boostpairs;
  FacePair newpair;
  graph_traits<Graph>::vertex_iterator vi, vi_end;
  for(tie(vi,vi_end) = vertices(graph); vi != vi_end; ++vi) {
    if (mate[*vi] != graph_traits<Graph>::null_vertex() && *vi < mate[*vi]) {
      newpair.first  = *vi;
      newpair.second = mate[*vi];
      boostpairs.push_back( newpair );
    }
  }
  sort( boostpairs.begin(), boostpairs.end() );

  Tri2Quads::collapse_matched_triangles(inmesh, boostpairs, 1 );

  ticks end_tick = getticks();

  double elapsed_ticks = elapsed(end_tick, start_tick);

  cout << "Elapsed Performance Counter " << elapsed_ticks << endl;
  cout << "Info: Storing Quad Mesh in boost.dat" << endl;
  cout << "Output Mesh: " << endl;
  cout << "      #Nodes   : " << inmesh->getSize(0) << endl;
  cout << "      #Quads   : " << inmesh->getSize(2) << endl;

  inmesh->get_topological_statistics();
  inmesh->saveAs( "boost.dat");
}

////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
  if( argc == 2) {
      edmond_graph_matching( argv[1] );
      return 1;
  } else {
      string nofile;
      edmond_graph_matching( nofile );
  }
  
  return 0;
}

