#ifndef EDGEFLIP_H
#define EDGEFLIP_H

#include <meshkit/Mesh.hpp>

using namespace Jaal;

class SwapTriEdge : public MeshOptimization {
public:
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

    virtual int execute();

    void setConstraintEdges(vector<Edge*> &edges) {
        //	  constraint_edges.add(emesh);
    }

    size_t get_number_of_edges_flipped() const {
        return num_edges_flipped;
    }

    int apply_reduce_degree_rule();
    int apply_advance_front_rule();
    int apply_delaunay_rule();

protected:
    Mesh *mesh;
    double creaseAngle;
    size_t num_edges_flipped;
    //   EdgeMap constraint_edges;

    struct FlipEdge : public Edge {

        FlipEdge() {
        }

        FlipEdge(Vertex *v1, Vertex * v2) {
            process(v1, v2);
        }

        ~FlipEdge() {
        }

        bool isSharp(double creseAngle) const;
        bool isConcave() const;

        Face * faces[2];
        Vertex * opposite_nodes[2];

    private:
        void process(Vertex *v1, Vertex * v2);
    };

    int atomicOp(const Face *face);
    virtual int commit(const FlipEdge &edge);
    virtual bool is_edge_flip_allowed(const FlipEdge &edge) const;
};

class VertexDegreeReduction : public SwapTriEdge {
public:

    VertexDegreeReduction(Mesh *m) : SwapTriEdge(m) {
    }

    size_t getMinDegree() const;
    size_t getMaxDegree() const;

    int execute();

private:
    int atomicOp(const Vertex *apexVertex);
    bool is_edge_flip_allowed(const FlipEdge &edge) const;
    int getVertexDegree(const Vertex *v) const;
};

#endif
