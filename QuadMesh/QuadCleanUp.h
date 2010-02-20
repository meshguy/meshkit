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
// For suggestios, bugs and criticisms, please send e-mail to
//                      csverma@cs.wisc.edu
//
// Last Date update:  16th Feb 2010.
//
///////////////////////////////////////////////////////////////////////////////////

#include "Mesh.h"

BEGIN_JAAL_NAMESPACE

///////////////////////////////////////////////////////////////////////////////////
// Diamond:  An element whose at least one of the opposite vertex is surrounded by
//           three faces. In many cases, diamonds are essential in the quadrilateral
// mesh and they can not be removed, Finding the minimum number of diamonds is hard,
// and we are working towards it.
///////////////////////////////////////////////////////////////////////////////////

struct Diamond
{
  Face *face;
  Vertex *vertex0, *vertex1;
};

/////////////////////////////////////////////////////////////////////////////////////
//Bridge:  An Edge whose two end vertices are surrounded by three faces. By removing
//         a bridge, we essentially remove two diamonds. But the removal method is
//         different. In the bridge we use element removal followed by edge swapping.
/////////////////////////////////////////////////////////////////////////////////////

struct Bridge
{
  Vertex *vertex0, *vertex1;
};


class QuadCleanUp
{
public:
  QuadCleanUp(Mesh *m)
  {
    mesh = m;
  }

  // Query methods ...
  vector<Face*> search_diamonds(bool check_both_sides = 1,
      bool allow_boundary_faces = 1);
  vector<Face*> search_flat_quads();
  vector<Vertex*> search_interior_doublets();
  vector<Vertex*> search_boundary_singlets();
  vector<Bridge> search_bridges();

  // Removal Methods ...
  void remove_diamonds(bool recursive = 1, bool check_both_sides = 1,
      bool allow_boundary_faces = 1);
  void remove_doublets(bool recursive = 1, bool allow_boundary_nodes = 0);
  void remove_bridges();
  void cleanup_boundary(double cutOffAngle = 100.0);

  // Insert methods ...
  Vertex* insert_doublet(Face *face);
  Vertex* insert_boundary_doublet(Face *face);
  Vertex* insert_doublet(Face *face, Vertex *v0, Vertex *v2);

  // Utility functions ...
  void get_strips(Face *face, vector<Face*> &strip1, vector<Face*> strip2);

  // Topological Quality method ...
  vector<int> getVertexFaceDegrees();

private:
  // Input-output instance. Input mesh is modified...
  Mesh *mesh;

  vector<Diamond> vDiamonds; // Diamonds in the mesh;
  vector<Bridge> vBridges; // Bridges in the mesh.

  // Basic Operations ...
  int face_close(Face *face, Vertex *v0, Vertex *v2);
  int remove_interior_doublet(Vertex *vertex);
  int remove_boundary_singlet(Vertex *vertex);
  int diamond_collapse(Diamond &d);
  int remove_bridge( const Bridge &b);
  int remove_diamonds_once(bool check_both_sides = 1, bool allow_boundary_faces = 1);
  int remove_doublets_once( bool allow_boundary_nodes = 0);

  // High level utility function composed of basic functions...
  void cleanup_internal_boundary_face();

  // Create wavefront of nodes/faces ...
  void initialize_wavefront();
  vector<Vertex*> next_front_nodes() const;

  // Quality: Set the tag for regular (= 0)/irregular node (> 0) value
  void set_regular_node_tag();

  // Get the histogram of vertex-face topological information. (ideal is 4)
  Vertex* get_VertexOf_FaceDegree(int n);
};

END_JAAL_NAMESPACE

#endif

///////////////////////////////////////////////////////////////////////////////
