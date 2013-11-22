#ifndef QPATCHREMESH_H
#define QPATCHREMESH_H

#include "Mesh.hpp"

namespace Jaal {
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

/////////////////////////////////////////////////////////////////////////////////////

struct QuadRemeshTemplate : public RemeshTemplate {
     int remesh(Mesh *mesh,
                NodeSequence &anodes,
                NodeSequence &bnodes,
                NodeSequence &cnodes,
                NodeSequence &dnodes);
private:
     int remesh();
     int nx, ny;
     NodeSequence boundnodes, qnodes;
};

/////////////////////////////////////////////////////////////////////////////////////

struct TriRemeshTemplate : public RemeshTemplate {
     int remesh(Mesh *mesh, NodeSequence &anodes,
                NodeSequence &bnodes,
                NodeSequence &cnodes, int *partition);
private:
     QuadRemeshTemplate quadtemplate;
     int segments[3], partSegments[6];

     NodeSequence a1nodes, b1nodes;
     NodeSequence a2nodes, b2nodes;
     NodeSequence a3nodes, b3nodes;
     NodeSequence oa_nodes, ob_nodes, oc_nodes;
};

/////////////////////////////////////////////////////////////////////////////////////

struct PentaRemeshTemplate : public RemeshTemplate {
     int  remesh(Mesh *mesh,
                 NodeSequence &anodes,
                 NodeSequence &bnodes,
                 NodeSequence &cnodes,
                 NodeSequence &dnodes,
                 NodeSequence &enodes,
                 int *partition);
private:
     QuadRemeshTemplate quadtemplate;

     int segments[5], partSegments[10];
     NodeSequence a1nodes, b1nodes;
     NodeSequence a2nodes, b2nodes;
     NodeSequence a3nodes, b3nodes;
     NodeSequence a4nodes, b4nodes;
     NodeSequence a5nodes, b5nodes;
     NodeSequence oa_nodes, ob_nodes, oc_nodes, od_nodes, oe_nodes;
};

/////////////////////////////////////////////////////////////////////////////////////

struct CircleRemeshTemplate : public RemeshTemplate {
     int  remesh(Mesh *mesh, NodeSequence &anodes);

};
}

#endif

