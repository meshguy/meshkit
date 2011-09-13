#ifndef EDGEFLIP_H
#define EDGEFLIP_H

#include "Mesh.hpp"

using namespace Jaal;

class SwapTriEdge : public MeshOptimization {

public:
     static const int  DELAUNAY_RULE = 0;
     static const int  DEGREE_REDUCTION_RULE = 1;
     static const int  ADVANCE_FRONT_RULE = 2;

     //!  Constructor ...

     SwapTriEdge(Mesh *m, double angle = 10.0) {
          mesh = m;
          creaseAngle = angle;
     }

     ~SwapTriEdge() {
     }

     void setCreaseAngle(double a) {
          creaseAngle = a;
     }

     void setConstraintEdges(vector<Edge*> &edges) {
          //	  constraint_edges.add(emesh);
     }

     size_t get_number_of_edges_flipped() const {
          return num_edges_flipped;
     }

     int apply_rule(int r = DELAUNAY_RULE);

private:
     Mesh *mesh;
     double creaseAngle;
     size_t num_edges_flipped;

     int apply_advance_front_rule();

     bool  isIdeal( const Vertex *v)  const {
          int ideal_degree = v->get_ideal_face_degree(3);
          int curr_degree  = v->getNumRelations(2);
          if( curr_degree == ideal_degree ) return 1;
          return 0;
     }

     bool unchecked( const Face *f ) const {
          int lid = 0;
          for( int i = 0; i < 3; i++) {
               Vertex *v = f->getNodeAt(i);
               v->getAttribute("Layer", lid );
               if( lid == INT_MAX) return 1;
          }
          return 0;
     }

     struct FlipEdge : public Edge {

          FlipEdge() {
          }

          FlipEdge(Vertex *v1, Vertex * v2) {
               process(v1, v2);
          }

          ~FlipEdge() {
          }

          bool isValid() const {
               if( faces[0] == NULL || faces[1] == NULL ) return 0;
               if( opposite_nodes[0] == NULL || opposite_nodes[1] == NULL ) return 0;
               return 1;
          }

          bool isSharp(double creseAngle) const;
          bool isConcave() const;

          Face * faces[2];
          Vertex * opposite_nodes[2];
private:
          void process(Vertex *v1, Vertex * v2);
     };
     int one_sweep( int entity, int r);
     int atomicOp(const Face *face, int r);
     int atomicOp(Vertex *v, int r);
     virtual int commit(const FlipEdge &edge);
     virtual bool is_edge_flip_allowed(const FlipEdge &edge, int r ) const;
};

#endif
