#include "meshkit/DijkstraShortestPath.hpp"
#include "meshkit/StopWatch.hpp"
#include <limits>
///////////////////////////////////////////////////////////////////////////////

int DijkstraShortestPath::atomicOp(LVertex &currnode)
{
     Vertex *vi = currnode.vertex;

     if (vi == vdst) return 0;

     if( filter ) {
          if( vi != vsrc  ) {
               if( !filter->pass( vi ) ) {
                    vdst = vi;
                    return 0;
               }
          }
     }

     NodeSequence vneighs;
     vi->getRelations( vneighs );
     int nSize = vneighs.size();

//   int offset = Math::random_value((int)0, nSize-1);
     int offset = 0;

     LVertex lv;
     for (int i = 0; i < nSize; i++) {
          Vertex *vj = vneighs[(i+offset)%nSize];
          assert( !vj->isRemoved() );
          miter = vmap.find( vj );
          if( miter == vmap.end()) {
               lv.distance = std::numeric_limits< double >:: max();
               lv.vertex   = vj;
               lv.previous = vi;
               vmap.insert(make_pair(vj,lv));
          } else {
               lv = miter->second;
          }
          assert( lv.vertex == vj );

          double vcost = getCost(currnode, lv);
          if (vcost < lv.distance) {
               lv.previous = vi;
               lv.distance = vcost;
               vmap[vj]    = lv;
               vertexQ.push(lv);
          }
     }

     return 1;
}

///////////////////////////////////////////////////////////////////////////////

void DijkstraShortestPath::traceback()
{
     nodepath.clear();
     if( vdst == NULL ) return;

     miter = vmap.find(vdst);
     if( miter == vmap.end() ) return;

     LVertex currnode = miter->second;
     while(1) {
          nodepath.push_back( currnode.vertex );
          Vertex *v = currnode.previous;
          if (v == NULL) break;
          currnode = vmap[v];
     }

     assert( nodepath.front() == vdst );
     assert( nodepath.back()  == vsrc );
}
///////////////////////////////////////////////////////////////////////////////

void DijkstraShortestPath::fastmarching()
{
     assert( mesh->getAdjTable(0,0) );

     while (!vertexQ.empty()) vertexQ.pop();

     vmap.clear();

     LVertex lv;
     lv.vertex   = vsrc;
     lv.distance = 0.0;
     lv.previous = NULL;
     vmap[vsrc]  = lv;
     vertexQ.push( lv );

     assert( !vsrc->isRemoved() );

     int progress;
     while (!vertexQ.empty()) {
          LVertex currVertex = vertexQ.top();
          vertexQ.pop();
          progress = atomicOp(currVertex);
          if (!progress) break;
     }
}
///////////////////////////////////////////////////////////////////////////////

int Jaal::dijkstra_shortest_path_test()
{
     double origin[]   = { 0.0, 0.0, 0.0};
     double length[]   = { 1.0, 1.0, 1.0};
     int    gridim[]   = { 5, 5, 5};

     Jaal::Mesh *mesh = Jaal::create_structured_mesh(origin, length, gridim, 2 );

     struct MyFilter: public MeshFilter {
          MyFilter( int i ) {
               id = i;
          }
          size_t id;
          bool pass( const Vertex *v) const{
               return v->getID() != id;
          }
     };

     MeshFilter *filter = new MyFilter(18);

     DijkstraShortestPath djk(mesh);

     NodeSequence sq = djk.getPath(mesh->getNodeAt(0), filter);

     for( size_t i = 0; i < sq.size(); i++)
          cout << sq[i]->getID() << endl;

     delete filter;

     return 0;
}

///////////////////////////////////////////////////////////////////////////////

