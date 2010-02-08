#ifndef QUADCLEAN_H
#define QUADCLEAN_H

#include "Mesh.h"

BEGIN_JAAL_NAMESPACE

class QuadCleanUp
{
  public:
     QuadCleanUp(Mesh *m) { mesh = m; }

     void get_strips(Face *face, vector<Face*> &strip1, vector<Face*> strip2);

     vector<int> getVertexFaceDegrees();

     vector<Face*> search_diamonds();
     vector<Vertex*> search_doublets();

     void remove_diamonds();
     void remove_doublets();
     void cleanup_boundary( double cutOffAngle = 100.0);

     Vertex* insert_doublet(Face *face);
     Vertex* insert_boundary_doublet(Face *face);
     Vertex* insert_doublet(Face *face, Vertex *v0, Vertex *v2);

  private:
     Mesh *mesh;
     Vertex* get_VertexOf_FaceDegree(int n);

     int face_close(Face *face, Vertex *v0, Vertex *v2);

     void cleanup_internal_boundary_face();

     int remove_unconstrained_doublet(Vertex *vertex);
     int diamond_collapse(Face *face);

     void initialize_wavefront();
     vector<Vertex*> next_front_nodes() const;
};

END_JAAL_NAMESPACE

#endif
