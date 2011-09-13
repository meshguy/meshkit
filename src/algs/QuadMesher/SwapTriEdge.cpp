#include "SwapTriEdge.hpp"

#include <iostream>
#include <iomanip>

//#############################################################################

int SwapTriEdge:: one_sweep(int entity, int rule)
{
     size_t total_count = 0;
     size_t numnodes = mesh->getSize(0);
     size_t numfaces = mesh->getSize(2);

     // Sweep over nodes in the mesh...
     if( entity == 0 ) {
          while(1) {
               size_t curr_count = 0;
               for (size_t i = 0; i < numnodes; i++) {
                    int err = atomicOp( mesh->getNodeAt(i), rule);
                    if (!err) curr_count++;
               }
               if( curr_count == 0) break;
               total_count += curr_count;
          }
     }

     // Sweep over faces in the mesh...
     if( entity == 2 ) {
          while(1) {
               size_t curr_count = 0;
               for (size_t i = 0; i < numfaces; i++) {
                    int err = atomicOp( mesh->getFaceAt(i), rule);
                    if (!err) curr_count++;
               }
               if( curr_count == 0) break;
               total_count += curr_count;
          }
     }
 
     assert( numnodes ==  mesh->getSize(0) );
     assert( numfaces ==  mesh->getSize(2) );

     mesh->make_consistently_oriented();

     return total_count;
}

//#############################################################################

int SwapTriEdge::apply_rule(int rule )
{
     int relexist2 = mesh->build_relations(0, 2);

     mesh->search_boundary();

     assert(mesh->getAdjTable(0, 2));

     num_edges_flipped = 0;
     Jaal::MeshOptimization mopt;

     int curr_count;
     while(1) {
          curr_count = one_sweep( 2, rule);
          if( curr_count == 0) {
               mopt.shape_optimize(mesh);
               curr_count = one_sweep( 2, rule );
               if( curr_count == 0) break;
               num_edges_flipped+= curr_count;
          } else
               num_edges_flipped+= curr_count;
     }

     if (!relexist2) mesh->clear_relations(0, 2);

     cout << "# of edges swapped  " << num_edges_flipped << endl;

     return num_edges_flipped;
}

//#############################################################################

int SwapTriEdge::atomicOp(const Face *face, int rule)
{
     if( !face->isActive() )  return 1;

     for (int i = 0; i < 3; i++) {
          Vertex *v1 = face->getNodeAt(i + 1);
          Vertex *v2 = face->getNodeAt(i + 2);
          FlipEdge edge(v1, v2);
          if (is_edge_flip_allowed(edge, rule)) {
               int err = commit(edge);
               if( !err) return 0;
          }
     }
     return 2;
}

//#############################################################################

void SwapTriEdge::FlipEdge::process(Vertex *v1, Vertex *v2)
{
     assert(v1 != v2);

     this->setNodes(v1, v2);

     faces[0] = NULL;
     faces[1] = NULL;
     opposite_nodes[0] = NULL;
     opposite_nodes[1] = NULL;

     if (v1->isBoundary() && v2->isBoundary()) return;

     FaceSequence neighs;
     Mesh::getRelations112(v1, v2, neighs);

     assert(neighs.size() == 2);

     Vertex *ov1 = Face::opposite_node(neighs[0], v1, v2);
     Vertex *ov2 = Face::opposite_node(neighs[1], v1, v2);
     assert(ov1 && ov2);

     faces[0] = neighs[0];
     faces[1] = neighs[1];
     opposite_nodes[0] = ov1;
     opposite_nodes[1] = ov2;
}

//#############################################################################

bool SwapTriEdge::FlipEdge::isSharp(double featureAngle) const
{
     Vertex *v1 = getNodeAt(0);
     Vertex *v2 = getNodeAt(1);
     Vertex *ov1 = opposite_nodes[0];
     Vertex *ov2 = opposite_nodes[1];

     Vec3D A, B;
     A = Face::normal(ov1, v1, v2);
     B = Face::normal(ov2, v2, v1);

     double angle = Math::getVectorAngle(A, B, ANGLE_IN_DEGREES);
     if (angle > featureAngle && fabs(180 - angle) > 1.0E-06) return 1;

     return 0;
}

//#############################################################################

bool SwapTriEdge::FlipEdge::isConcave() const
{
     Vertex *v1 = getNodeAt(0);
     Vertex *v2 = getNodeAt(1);
     Vertex *ov1 = opposite_nodes[0];
     Vertex *ov2 = opposite_nodes[1];

     assert( v1 && v2  );
     assert( ov1 && ov2  );

     bool convex = Face::is_convex_quad(v1->getXYZCoords(),
                                        ov1->getXYZCoords(),
                                        v2->getXYZCoords(),
                                        ov2->getXYZCoords());
     if (!convex) return 1;

     return 0;
}


//*****************************************************************************

bool SwapTriEdge::is_edge_flip_allowed(const FlipEdge &edge, int rule) const
{
     if (!edge.isValid())   return 0;
     if (edge.isConcave()) return 0;

     Vertex *v1  = edge.getNodeAt(0);
     Vertex *v2  = edge.getNodeAt(1);
     Vertex *ov1 = edge.opposite_nodes[0];
     Vertex *ov2 = edge.opposite_nodes[1];

     if (v1->isBoundary()  && v2->isBoundary()  ) return 0;
     if (ov1->isBoundary() || ov2->isBoundary() ) return 0;

     int d1 = v1->getNumRelations(2);
     int d2 = v2->getNumRelations(2);
     int d3 = ov1->getNumRelations(2);
     int d4 = ov2->getNumRelations(2);

     if( rule == DELAUNAY_RULE) {
          double len1 = Vertex::length2(v1, v2);
          double len2 = Vertex::length2(ov1, ov2);
          if (len1 <= len2) return 0;

          if (edge.isSharp(creaseAngle)) return 0;
          return 1;
     }

     if( rule == ADVANCE_FRONT_RULE ) {
          int ideal_v1 =  v1->get_ideal_face_degree(3);
          int ideal_v2 =  v2->get_ideal_face_degree(3);
          int ideal_v3 =  ov1->get_ideal_face_degree(3);
          int ideal_v4 =  ov2->get_ideal_face_degree(3);

          int l1 = 0;
          v1->getAttribute("Layer",  l1);

          int l2 = 0;
          v2->getAttribute("Layer",  l2);

          int l3 = 0;
          ov1->getAttribute("Layer", l3);

          int l4 = 0;
          ov2->getAttribute("Layer", l4);

          // Decrease the vertex degree
          if( (d1  > ideal_v1) && (l2 > l1) && (l3 > l1) && (l4 > l1) ) return 1;
          if( (d2  > ideal_v2) && (l1 > l2) && (l3 > l2) && (l4 > l2) ) return 1;

          // Increase the vertex degree ...
          if( (d3  < ideal_v3) && (l1 > l3) && (l2 > l3) && (l4 > l3) ) return 1;
          if( (d4  < ideal_v4) && (l1 > l4) && (l2 > l4) && (l3 > l4) ) return 1;

          return 0;
     }

     if( rule == DEGREE_REDUCTION_RULE) {
          int ideal_v1 =  v1->get_ideal_face_degree(3);
          int ideal_v2 =  v2->get_ideal_face_degree(3);

          if( v1->isBoundary() &&  (d1 > ideal_v1) && (d2 >3) )  return 1;
          if( v2->isBoundary() &&  (d2 > ideal_v2) && (d1 >3) )  return 1;

          int relaxation_index = d1 + d2 - d3 - d4;
          if (relaxation_index < 3) return 0;
          return 1;
     }

     return 0;
}

//#############################################################################

int SwapTriEdge::commit(const FlipEdge &edge)
{
     Face *t1 = edge.faces[0];
     Face *t2 = edge.faces[1];
     Vertex *v1 = edge.getNodeAt(0);
     Vertex *v2 = edge.getNodeAt(1);
     Vertex *ov1 = edge.opposite_nodes[0];
     Vertex *ov2 = edge.opposite_nodes[1];

     int pos1 = t1->getPosOf(v1);
     int pos2 = t2->getPosOf(v1);

     v1->removeRelation(t1);
     v1->removeRelation(t2);
     v2->removeRelation(t1);
     v2->removeRelation(t2);
     v1->removeRelation(v2);
     v2->removeRelation(v1);

     NodeSequence vconn(3);

     if( t1->getNodeAt(pos1+1) == ov1 && t2->getNodeAt(pos2+2) == ov2 ) {
          vconn[0] = v1;
          vconn[1] = ov1;
          vconn[2] = ov2;
          t1->setNodes(vconn);

          vconn[0] = v2;
          vconn[1] = ov2;
          vconn[2] = ov1;
          t2->setNodes(vconn);

          mesh->reactivate(t1);
          mesh->reactivate(t2);
          return 0;
     }

     if( t1->getNodeAt(pos1+2) == ov1 && t2->getNodeAt(pos2+1) == ov2 ) {
          vconn[0] = v1;
          vconn[1] = ov2;
          vconn[2] = ov1;
          t1->setNodes(vconn);

          vconn[0] = v2;
          vconn[1] = ov1;
          vconn[2] = ov2;
          t2->setNodes(vconn);

          mesh->reactivate(t1);
          mesh->reactivate(t2);
          return 0;
     }

     cout << "Fatal Error: Orientation of faces unknown " << endl;
     exit(0);

     return 0;
}

//#############################################################################

int SwapTriEdge ::atomicOp(Vertex *apexVertex, int rule)
{
     FaceSequence vneighs;
     apexVertex->getRelations( vneighs );
     int numNeighs = vneighs.size();

     for( int i = 0; i < numNeighs; i++) {
          if( unchecked(vneighs[i] ) ) {
               int err = atomicOp( vneighs[i], rule);
               if( err  == 0) return 0;
          }
     }

     return 1;
}
///////////////////////////////////////////////////////////////////////////
int SwapTriEdge :: apply_advance_front_rule()
{
     int relexist2 = mesh->build_relations(0, 2);
     int relexist0 = mesh->build_relations(0, 0);

     mesh->search_boundary();

     size_t numNodes = mesh->getSize(0);
     NodeSequence currlayer, vneighs;
     NodeSet   nextlayer;

     for(size_t i = 0; i < numNodes; i++) {
          Vertex *v = mesh->getNodeAt(i);
          if( v->isBoundary() ) {
               v->setAttribute("Layer", 0);
               currlayer.push_back(v);
          } else
               v->setAttribute("Layer",INT_MAX);
     }

     Jaal::MeshOptimization mopt;
     size_t ncount, nSize, num_edges_flipped = 0;
     mesh->make_consistently_oriented();
     mopt.shape_optimize(mesh);

     LaplaceNoWeight lw;
     LaplaceSmoothing lapsmooth(mesh);
     lapsmooth.setWeight(&lw);
     lapsmooth.setNumIterations(100);

     int curr_layer_id = 0;
     int nirregular0, nirregular1;

     while(1) {
          nSize = currlayer.size();
          nirregular0 = 0;
          for( size_t i = 0; i < nSize; i++) {
               Vertex *v = currlayer[i];
               if( !v->isRemoved() && !isIdeal(v) ) nirregular0++;
          }

          for( int k = 0; k < 3; k++) {  // Three attempts for single layer
               while(1) {
                    nSize = currlayer.size();
                    ncount = 0;
                    for( size_t i = 0; i < nSize; i++) {
                         Vertex *v = currlayer[i];
                         if( !v->isRemoved() &&  !isIdeal(v)) {
                              int err  = atomicOp( v, ADVANCE_FRONT_RULE );
                              if( !err ) ncount++;
                         }
                    }
                    if( ncount == 0) break;
                    num_edges_flipped += ncount;
               }

               nSize = currlayer.size();
               nirregular1 = 0;
               for( size_t i = 0; i < nSize; i++) {
                    Vertex *v = currlayer[i];
                    if( !v->isRemoved() &&  !isIdeal(v) ) nirregular1++;
               }
               assert( nirregular1 <= nirregular0);
               if( nirregular1 == 0) break;
          }
          lapsmooth.execute();

//        mopt.shape_optimize(mesh);

          cout << "Layer : " << curr_layer_id << endl;
          cout << "# of Irregular nodes before swapping : " << nirregular0 << endl;
          cout << "# of Irregular nodes after swapping  : " << nirregular1 << endl;

          int lid;

          nextlayer.clear();
          nSize = currlayer.size();
          for( size_t i = 0; i < nSize; i++) {
               Vertex *v = currlayer[i];
               v->getRelations( vneighs );
               for( size_t k = 0; k < vneighs.size(); k++) {
                    vneighs[k]->getAttribute( "Layer", lid);
                    if( lid > curr_layer_id ) {
                         vneighs[k]->setAttribute("Layer", curr_layer_id+1 );
                         nextlayer.insert( vneighs[k] );
                    }
               }
          }
          if( nextlayer.empty() ) break;

          NodeSet::const_iterator it;
          currlayer.resize(nextlayer.size() );
          int index = 0;
          for( it = nextlayer.begin(); it != nextlayer.end(); ++it)
               currlayer[index++] = *it;
          curr_layer_id++;
     }
     cout << "# of Edges Swapped " << num_edges_flipped << endl;

     vector<int>  less_than_ideal, more_than_ideal, total_ideal;

     int numLayers = curr_layer_id;
     less_than_ideal.resize( numLayers );
     more_than_ideal.resize( numLayers );
     total_ideal.resize( numLayers );

     for( int i = 0; i < numLayers; i++) {
          less_than_ideal[i] = 0;
          more_than_ideal[i] = 0;
          total_ideal[i] = 0;
     }
     int lid;

     numNodes = mesh->getSize(0);
     int final_irregular = 0;
     for( size_t i = 0; i < numNodes; i++) {
          Vertex *v = mesh->getNodeAt(i);
          if( !v->isRemoved()) {
               v->getAttribute("Layer", lid );
               int curr_degree  = v->getNumRelations(2);
               int ideal_degree = v->get_ideal_face_degree(3);
               if( curr_degree != ideal_degree ) {
                    final_irregular++;
                    if( curr_degree < ideal_degree) less_than_ideal[lid]++;
                    if( curr_degree > ideal_degree) more_than_ideal[lid]++;
               } else
                    total_ideal[lid]++;
          }
     }

     cout << " Layer   Less   More  Ideal " << endl;
     for( int i = 0; i < numLayers; i++)
          cout << i << setw(10) <<  less_than_ideal[i]
               << setw(10) <<  more_than_ideal[i]
               << setw(10) <<  total_ideal[i] << endl;
     cout << " Final # of irregular nodes : " << final_irregular << endl;

     mopt.shape_optimize(mesh);

     if (!relexist2) mesh->clear_relations(0, 2);
     if (!relexist0) mesh->clear_relations(0, 0);

     return num_edges_flipped;
}

