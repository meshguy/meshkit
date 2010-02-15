#ifndef QUADCLEAN_H
#define QUADCLEAN_H

#include "Mesh.h"

BEGIN_JAAL_NAMESPACE

struct Diamond
{
  Face *face;
  Vertex *vertex0, *vertex1;
};

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
  void remove_doublets(bool recursive = 1);
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
