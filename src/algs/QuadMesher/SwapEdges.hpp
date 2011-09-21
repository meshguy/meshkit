#ifndef EDGEFLIP_H
#define EDGEFLIP_H

#include "Mesh.hpp"

namespace Jaal
{
struct SwapEdges : public MeshOptimization
{
     static const int  DELAUNAY_RULE = 0;
     static const int  DEGREE_REDUCTION_RULE = 1;
     static const int  ADVANCE_FRONT_RULE = 2;

     void setCreaseAngle(double a) {
          creaseAngle = a;
     }

     void setConstraintEdges(vector<Edge*> &edges) {
          //	  constraint_edges.add(emesh);
     }

     size_t get_number_of_edges_flipped() const {
          return num_edges_flipped;
     }

     Mesh *mesh;
     double creaseAngle;
     size_t num_edges_flipped;

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

          FaceSequence  oldFaces, newFaces;
          NodeSequence  oldNodes, newNodes;
       private:
          void process(Vertex *v1, Vertex * v2);
     };
};

//////////////////////////////////////////////////////////////////////////////

class SwapTriEdges : public SwapEdges {

public:

     //!  Constructor ...

     SwapTriEdges(Mesh *m, double angle = 10.0) {
          mesh = m;
          creaseAngle = angle;
     }

     ~SwapTriEdges() { }

     int apply_rule(int r = DELAUNAY_RULE);

private:

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
     int one_sweep( int entity, int r);
     int atomicOp(const Face *face, int r);
     int atomicOp(Vertex *v, int r);

     virtual int commit(const FlipEdge &edge);
     virtual bool is_edge_flip_allowed(const FlipEdge &edge, int r ) const;
};


struct QuadEdge {

     QuadEdge() {
          mesh = NULL;
          connect[0] = NULL;
          connect[1] = NULL;
     }

     ~QuadEdge() {
          for( size_t i = 0; i < newNodes.size(); i++)
               if( newNodes[i] ) delete newNodes[i];

          for( size_t i = 0; i < newFaces.size(); i++)
               if( newFaces[i] ) delete newFaces[i];
     }

     bool isBoundary() const {
          if (adjFaces[0] == NULL || adjFaces[1] == NULL) return 1;
          return 0;
     }

     Vertex*  getNodeAt( int i ) const {
          if( i == 0) return connect[0];
          if( i == 1) return connect[1];
          return NULL;
     }

     Vertex * connect[2];
protected:
     Mesh *mesh;
     FaceSequence  adjFaces;
     NodeSequence  newNodes;
     FaceSequence  newFaces;
     int commit();
};


class SwapQuadEdge : public QuadEdge {
public:
     static bool is_topologically_valid_swap(int d1, int d2, int d3, int d4);

     SwapQuadEdge(Mesh *m, Vertex *v0, Vertex *v1, Face *firstface = NULL) {
          mesh = m;
          edge = NULL;
          connect[0] = v0;
          connect[1] = v1;
          firstFace = firstface;
          assert(mesh);
          check_fronts = 0;
     }

     SwapQuadEdge(Mesh *m, Edge *e, Face *firstface = NULL) {
          mesh = m;
          edge = e;
          connect[0] = e->getNodeAt(0);
          connect[1] = e->getNodeAt(1);
          firstFace = firstface;
          assert(mesh);
     }

     void clear() {
          if (newFaces[0] != NULL) delete newFaces[0];
          if (newFaces[1] != NULL) delete newFaces[1];
     }

     void modify_fronts( int v ) {
          check_fronts = v;
     }

     int apply_reduce_degree_rule();
     int apply_concave_rule();    // Swap some concave edge
     int apply_bound_rule();     // Swap some boundary edge
     int apply_singlet_rule(Vertex *singlet); // Force creating diagonal at singlet..
     int apply_deficient_rule(Vertex *v); // Force creating diagonal at deficient vertex..
     int apply_advance_front_rule();

private:
     bool  check_fronts;
     Face *firstFace;           // Which one of the two faces is the first one. It is needed
     Edge *edge;                // Swapping edge.
     NodeSequence bound_nodes;  // It is always going to be six nodes...

     struct BackUpData {
          Vertex *diagonalConnect[2]; // Diagonal of the two quads.
          NodeSequence face1Connect, face2Connect;
          map<Vertex*, Point3D>  pCoords; // Six Boundary node's coordinates.
     };
     BackUpData bkp_data;

     void backup();
     void rollback();

     int getPosOf(const Vertex *v) const;
     int build_boundary();
     int get_boundary_nodes_chain();
     int make_new_diagonal_at(int pos, bool bound_check = 1);

     // For front cleaning ops;
     int hasLessNodes(Vertex *vertex, int layerid);
     int hasExcessNodes(Vertex *vertex, int layerid);
     void update_front();
};
}

/////////////////////////////////////////////////////////////////////////////////////

/*
struct RemeshTemplate {
     NodeSequence newnodes;
     FaceSequence newfaces;

     void addNewElements( const NodeSequence &nnodes, const FaceSequence &nfaces) {
          size_t nSize;
          nSize = nnodes.size();
          for (size_t i = 0; i < nSize; i++)
               newnodes.push_back(nnodes[i]);

          nSize = nfaces.size();
          for (size_t i = 0; i < nSize; i++)
               newfaces.push_back(nfaces[i]);
     }

     void discard() {
          size_t nSize = newfaces.size();
          for (size_t i = 0; i < nSize; i++)
               mesh->remove(newfaces[i]);

          nSize = newnodes.size();
          for (size_t i = 0; i < nSize; i++)
               mesh->remove(newnodes[i]);
     }

     Mesh   *mesh;
};
*/



#endif
