#ifndef QUADCLEAN_H
#define QUADCLEAN_H

////////////////////////////////////////////////////////////////////////////////
//                      Quad-Cleanup
//
// Developed by:  Chaman Singh Verma
//                Department of Computer Sciences.
//                The University of Wisconsin, Madison
//
// Work Supported by:
//                 Dr. Tim Tautges
//                 Argonne National Lab, Chicago
//
//
// Objective:  Given a quadrilateral mesh, this class implements various strategies
// to improve the quadrilateral mesh both geometrically and topologically. The
// Laplacian ( local and global ) is used for geometric quality improvement, and for
// topological improvements various operations are used.
// The two basis operations for topological improvements are
//  1)   Face close
//  2)   doublet insertion and removal.
//

// Reference Papers:
//  1) Topological Improvement Procedures for Quadrilateral Finite Element Meshes
//     S.A. Canann,  S.N. Muthikrishnan and R.K. Phillips

//  2) Automated All Quadrilateral Mesh Adaptation Through Refinment and Coarsening
//     Bret Dallas Anderson
//     Master Thesis, Brigham Young University.
//
//  3) Non-Local Topological Clean-Up ( The idea of yring  is from this paper)
//     Guy Bunin.
//
// For suggestios, bugs and criticisms, please send e-mail to
//                      csverma@cs.wisc.edu
//
// Last Date update:  16th Feb 2010.
//
///////////////////////////////////////////////////////////////////////////////////

#include "Mesh.hpp"
#include "Tri2Quad.hpp"
#include "basic_math.hpp"
#include "StopWatch.hpp"
#include "tfiblend.hpp"

#include "DijkstraShortestPath.hpp"
#include "QuadPatchRemesh.hpp"

extern double area_of_poly3d(int n, double *x, double *y, double *z);

namespace Jaal {

struct FirstIrregularNode : public MeshFilter {
     bool pass( const Vertex *vertex ) const {
          if( vertex->isBoundary() ) return 1;
          if( vertex->getNumRelations(2) != 4 ) return 0;
          return 1;
     }
};

///////////////////////////////////////////////////////////////////////////////////
// Diamond:  An element whose at least one of the opposite vertex is surrounded by
//           three faces. In many cases, diamonds are essential in the quadrilateral
// mesh and they can not be removed, Finding the minimum number of diamonds is hard,
// and we are working towards it.
///////////////////////////////////////////////////////////////////////////////////

class FaceClose {
public:

     FaceClose(Mesh *m, Face *f, Vertex *v0, Vertex *v1) {
          mesh = m;
          face = f;
          vertex0 = v0;
          vertex2 = v1;
          replacedNode = NULL;
     }

     ~FaceClose() {
          if( replacedNode ) delete replacedNode;
     }

     int remove();

     int build();
     int commit();

     Mesh *mesh;
     Face *face;
     Vertex *vertex0, *vertex2;
     Vertex *replacedNode;
private:
     bool isSafe() const;
     int  backup();
     int  rollback();
};

///////////////////////////////////////////////////////////////////////////////////

struct Diamond {

     Diamond(Mesh *m, Face *f, int p) {
          mesh = m;
          face = f;
          position = p;
          faceclose = NULL;
          vertex0   = NULL;
          vertex2   = NULL;

          if (position == 0 || position == 2) {
               vertex0 = face->getNodeAt(0);
               vertex2 = face->getNodeAt(2);
          }

          if (position == 1 || position == 3) {
               vertex0 = face->getNodeAt(1);
               vertex2 = face->getNodeAt(3);
          }
     }
     ~Diamond() {
          if( faceclose) delete faceclose;
     }

     int remove();
     int commit();
     int isSafe();
     int makeShield();

     bool operator<(const Diamond & rhs) const {
          return face->getArea() < rhs.face->getArea();
     }

     Vertex * getNewNode() const {
          if (faceclose) return faceclose->replacedNode;
          return NULL;
     }

     double getDiagonalRatio() const {
          Vertex *v0 = face->getNodeAt((position + 0) % 4);
          Vertex *v1 = face->getNodeAt((position + 1) % 4);
          Vertex *v2 = face->getNodeAt((position + 2) % 4);
          Vertex *v3 = face->getNodeAt((position + 3) % 4);
          double len0 = Vertex::length(v0, v2);
          double len1 = Vertex::length(v1, v3);
          return len0 / len1;
     }
     int build();

     Face *face;
private:
     Vertex *vertex0, *vertex2;
     Mesh *mesh;
     int position;
     FaceClose *faceclose;
};

///////////////////////////////////////////////////////////////////////////////////

struct Doublet {

     Doublet(Mesh *m, Vertex * v) {
          mesh = m;
          vertex = v;
          replacedFace = NULL;
     }

     bool isSafe() const;
     void makeShield();
     int remove();

     Mesh *mesh;
     Vertex *vertex;
     Face *replacedFace;
     Face * shield[2];
};

///////////////////////////////////////////////////////////////////////////////////

struct Singlet {
     Singlet(Mesh *m, Vertex * v) {
          mesh = m;
          vertex = v;
          active = 1;
     }

     int remove();

private:
     Mesh *mesh;
     Vertex *vertex;
     bool active;
     NodeSequence oldNodes, newNodes;
     FaceSequence oldFaces, newFaces;

     int remove_by_refinement();
     int remove_by_swapping();

     int commit();
     void clear();
};

/////////////////////////////////////////////////////////////////////////////////////

class OneDefectPatch {
public:
     static size_t MAX_FACES_ALLOWED;

     static size_t num_boundaries;
     static double exec_time;
     static size_t num_3_patches;
     static size_t num_4_patches;
     static size_t num_5_patches;
     static size_t disk_remeshable;

     OneDefectPatch( Mesh *m ) {
          mesh = m;
          apex = NULL;
          quad_splitting_node = NULL;
          new_defective_node  = NULL;
     }

     OneDefectPatch( Mesh *m, Vertex *v) {
          mesh = m;
          apex = v;
          quad_splitting_node = NULL;
          new_defective_node  = NULL;
     }

     void set_initial_path( const NodeSequence &sq) {
          nodepath = sq;
     }

     size_t getSize(int e) const {
          if( e == 0) return inner_nodes.size() + bound_nodes.size();
          if( e == 2) return faces.size();
          return 0;
     }

     const NodeSequence &get_irregular_nodes_removed() {
          return irregular_nodes_removed;
     }

     bool isBoundaryEven() const {
          if( bound_nodes.size()%2 == 0 ) return 1;
          return 1;
     }

     size_t count_irregular_nodes(int where);

     int  build_remeshable_boundary();

     double get_isoperimetic_quotient() const {
          // This definiation is taken from Wikipedia..
          double A = getArea();
          double L = getPerimeter();
          double q = 4*M_PI*A/(L*L);
          return q;
     }

     void getFaces( FaceSequence &result) const {
          size_t nSize = faces.size();

          result.clear();
          if( nSize == 0 ) return;

          result.resize( nSize );

          int index = 0;
          FaceSet::const_iterator it;
          for( it = faces.begin(); it != faces.end(); ++it)
               result[index++] = *it;
     }

     const NodeSequence &getBoundaryNodes() const {
          return bound_nodes;
     }

     bool operator < ( const OneDefectPatch &rhs) const {
          return getSize(2) < rhs.getSize(2);
     }

     Vertex *get_new_defective_node() {
          return new_defective_node;
     }

     int  remesh();

     void setAttributes();

     void clear() {
          nodepath.clear();
          new_defective_node  = NULL;
          quad_splitting_node = NULL;
          corners.clear();
          inner_nodes.clear();
          bound_nodes.clear();
          faces.clear();
          inner_faces.clear();
          irregular_nodes_removed.clear();
          boundary.clear();
          cornerPos.clear();
          segSize.clear();
          newnodes.clear();
          newfaces.clear();
     }

private:
     bool isSafe();

     // Input data.
     Mesh   *mesh;
     Vertex *apex;               // Seed: Irregular vertex to start from.
     NodeSequence  nodepath;     // Initial joining two irregular nodes..

     FaceSet faces;
     FaceSet inner_faces;

#ifdef USE_HASHMAP
     std::tr1::unordered_map<Vertex*, FaceSet> relations02;
     std::tr1::unordered_map<Vertex*, FaceSet>::iterator miter;
//   std::tr1::unordered_set<Face*> inner_faces;
#else
     std::map<Vertex*, FaceSet> relations02;
     std::map<Vertex*, FaceSet>::iterator miter;
#endif

     Vertex *new_defective_node;

     // Local data ...
     Vertex *quad_splitting_node;     // One special node that splits a quad loop
     int quad_splitting_node_degree;  // Valence of the splitting node.

     NodeSet corners;                 // Corners of the blob
     NodeSet nodes;                   // All the nodes (inner + boundary)
     NodeSequence inner_nodes;        // Inner nodes (not on the boundary ) of the blob
     NodeSequence bound_nodes;        // Boundary nodes
     NodeSequence irregular_nodes_removed;

     vector<Edge> boundary;          // boundary of the blob.
     vector<int>  cornerPos;         // Positions of the corners in the bound_nodes.
     vector<int>  segSize;
     int   partSegments[10];

     TriRemeshTemplate    template3;
     QuadRemeshTemplate   template4;
     PentaRemeshTemplate  template5;

     // Variable used in 3-5 sided patch...
     NodeSequence anodes, bnodes, cnodes, dnodes, enodes;  // Nodes on each segment.

     // Variables used in 4 sided patch..
     NodeSequence a1nodes, a2nodes, b1nodes, c0nodes, c1nodes, c2nodes,
     abnodes, canodes, bcnodes, d1nodes;

     // New nodes and faces in the patch...
     NodeSequence  newnodes, nnodes;
     FaceSequence  newfaces, nfaces;

     // Get the position on the boundary ...
     int getPosOf( const Vertex *v) const {
          size_t nSize = bound_nodes.size();
          for (size_t i = 0; i <  nSize; i++)
               if (bound_nodes[i] == v) return i;
          return -1;
     }
     // Return nodes within the range (src, dst)
     void get_bound_nodes( const Vertex *src, const Vertex *dst, NodeSequence &s);

     // randomly select one irregular node
     bool  has_irregular_node_on_first_segment() const;

     // re-orient boundary nodes so that it starts from a given vertex.
     void start_boundary_loop_from (Vertex *v);

     // re-orient loops ...
     int reorient_4_sided_loop();

     // Patch creation functions...
     int   init_blob();
     int   update_boundary();
     int   finalize_boundary();
     int   expand_blob(Vertex *v);
     int   expand_blob();
     int   get_topological_outer_angle( Vertex *v);
     bool  is_simply_connected();

     bool  is_quad_breakable_at( const Vertex *v);
     // Query for the validity of 3-4-5 sided patches.
     bool  is_4_sided_convex_loop_quad_meshable();

     // Set the boundary pattern string.
     void set_boundary_segments();

     // If the resulting mesh is invalid for some reasons, revert back to
     // original and restore all information.
     void rollback();

     void pre_remesh();  // Before we start remeshing, do some clean-up
     int  remesh_3_sided_patch();
     int  remesh_4_sided_patch();
     int  remesh_5_sided_patch();
     void local_smoothing();
     void post_remesh(); // After successful remeshing, do some clean-up

     double getArea() const {
          double a = 0.0;
          FaceSet::const_iterator it;
          assert( faces.size() );
          for( it = faces.begin(); it != faces.end(); ++it) {
               a += fabs( (*it)->getArea() );
          }
          return a;
     }

     double getPerimeter() const {
          double l = 0.0;
          int nSize = bound_nodes.size();
          for( int i = 0; i < nSize; i++)
               l += Vertex::length( bound_nodes[i], bound_nodes[(i+1)%nSize] );
          return l;
     }

};

/////////////////////////////////////////////////////////////////////////////////////

class QuadCleanUp {
public:
     static bool isDoublet(const Vertex *v);
     static bool isSinglet(const Vertex *v);
     static bool isRegular( const Vertex *v);
     static bool hasSinglet(const Face *f);
     static bool isTunnel(const Edge *e);
     static bool isEdge33(const Edge *e);
     static bool isEdge35(const Edge *e);
     static bool isDiamond(Face *f, int &pos, int type = 33);

     QuadCleanUp() {
          djkpath = NULL;
          defective_patch = NULL;
     }

     QuadCleanUp(Mesh *m) {
          setMesh(m);
          djkpath = NULL;
          defective_patch = NULL;
     }

     ~QuadCleanUp() {
          /*
                    if( lapweight ) delete lapweight;
                    if( lapsmooth ) delete lapsmooth;
          */
          if( djkpath   ) delete djkpath;
          if( defective_patch ) delete defective_patch;
     }

     void setMesh( Mesh *m ) {
          mesh = m;
          /*
                    lapsmooth = new LaplaceSmoothing(mesh);
                    lapweight = new LaplaceLengthWeight;
                    lapsmooth->setWeight( lapweight );
          */
     }

     // Query methods ...
     NodeSequence  search_restricted_nodes();
     FaceSequence  search_restricted_faces();
     FaceSequence  search_flat_quads();

     vector<Diamond> search_diamonds(int type = 33 );
     vector<Singlet> search_boundary_singlets();
     vector<Doublet> search_interior_doublets();
     vector<Edge>    search_tunnels();
     vector<OneDefectPatch> search_one_defect_patches();

     OneDefectPatch* build_one_defect_patch(Vertex *vertex = NULL );

     int  degree_5_dominated();

     // Global Cleanup methods ..
     int remesh_defective_patches();

     // Local Cleanup methods ..
     int reduce_degree( Vertex *v );
     int vertex_degree_reduction();

     int swap_concave_faces();

     // Removal Methods ...
     int  remove_diamonds();
     int  remove_tunnels();
     int  remove_interior_doublets();
     int  remove_boundary_singlets();
     int  remove_bridges();
     int  shift_irregular_nodes();

//  int irregular_nodes_clustering();

     //  void remove_ynodes();
     int  clean_layer(int id);
     void cleanup_boundary(double cutOffAngle = 100.0);
     void advancing_front_cleanup();
     void advancing_front_edges_swap();

     int  automatic();
     void report();

     int atomic_op_swap_edge( Vertex *v0, Vertex *v1);
     int atomic_op_face_close( Face *f);

     // Some Feature that may be obsolete in the next version...
     Vertex* insert_doublet(Face *face);
     Vertex* insert_boundary_doublet(Face *face);
     Vertex* insert_doublet(Face *face, Vertex *v0, Vertex *v2);

//   vector<Edge33>  search_edges33();
//   int remove_edges35();

     int refine_restricted_node(Vertex *resnode, Vertex *bndnode);
     int refine_degree3_faces();
     int refine_bridges_face();
     // Utility functions ...
     void get_strips(Face *face, FaceSequence &strip1, FaceSequence strip2);

     int  reduce_internal_vertex_degree(Vertex *v);
     int  reduce_boundary_vertex_degree(Vertex *v);

private:
     Mesh *mesh;

     MeshOptimization mopt;
     /*
          LaplaceSmoothing *lapsmooth;
          LaplaceWeight *lapweight;
     */

     DijkstraShortestPath *djkpath; // Used in one defect remeshing ....
     OneDefectPatch* defective_patch;

     int  has_interior_nodes_degree_345();

     NodeSequence  irregular_nodes;

     vector<OneDefectPatch>  vDefectPatches;
     vector<Doublet> vDoublets;
     vector<Singlet> vSinglets;
     vector<Diamond> vDiamonds; // Diamonds in the mesh;
     vector<Diamond> search_diamonds_in_layer(int l);

     int clean_layer_once(int id);
     int face_close(Face *face, Vertex *v0, Vertex *v2);
     int diamond_collapse(FaceClose &d);
     int remove_interior_doublet(Doublet &d);
     int remove_boundary_singlet_type1(const Singlet &s);
     int remove_boundary_singlet_type2(const Singlet &s);
     int remove_boundary_singlets_once();
     int remove_bridges_in_layer( int l);
     int remove_bridges_once();
     int remove_diamonds_once();
     int remove_diamonds_in_layer( int l);
     int advance_front_edges_swap_once(int layerid);

     int apply_advance_front_bridge_rule( Vertex *v0, Vertex *v1);
     int apply_advance_front_excess_rule( Vertex *v);
     int apply_advance_front_triplet_rule( Vertex *v);
     int apply_advance_front_singlet_rule( Vertex *v);

     int remove_doublets_once();
     int remove_interior_doublets_once();

     int boundary_vertex_degree_reduction_once();
     int internal_vertex_degree_reduction_once();

     // High level utility function composed of basic functions...
     void cleanup_internal_boundary_face();

     // May become obsolere
     int refine_3434_pattern( Face *face, int pos);
     int refine_3454_pattern( Face *face, int pos);
     int refine_3444_pattern( Face *face, int pos);
     int apply_shift_node3_rule( Vertex *vertex);
};

////////////////////////////////////////////////////////////////////////////////

inline bool
QuadCleanUp::isRegular (const Vertex *v)
{
     assert(v);
     // Any interior vertex having four nodes( or faces ) is a regular node.
     if (!v->isBoundary () && (v->getNumRelations(2) == 4)) return 1;
     return 0;
}

////////////////////////////////////////////////////////////////////////////////

inline bool
QuadCleanUp::isDoublet (const Vertex *v)
{
     assert(v);
     // Any interior node having two neighboring face is a  doublet node.
     if (!v->isBoundary () && (v->getNumRelations(2) == 2)) return 1;
     return 0;
}

////////////////////////////////////////////////////////////////////////////////

inline bool
QuadCleanUp::isSinglet (const Vertex *v)
{
     assert( v );
     // Any boundary node having only one neigbour cell is a singlet node ...
     int numfaces = v->getNumRelations(2);
     assert (numfaces >= 0);
     if (v->isBoundary () && (numfaces == 1)) return 1;
     return 0;
}

////////////////////////////////////////////////////////////////////////////////

inline bool
QuadCleanUp::hasSinglet (const Face *face)
{
     assert( face );

     if( face->isRemoved() ) return 0;

     for (int i = 0; i < face->getSize (0); i++) {
          if (isSinglet (face->getNodeAt (i))) return 1;
     }
     return 0;
}

////////////////////////////////////////////////////////////////////////////////

void set_singlet_tag(Mesh *m, const string &s = "Singlet" );
void set_doublet_tag(Mesh *m, const string &s = "Doublet" );
void set_diamond_tag(Mesh *mesh, const string &s = "Diamond" );

} // namespace Jaal

#endif

///////////////////////////////////////////////////////////////////////////////
