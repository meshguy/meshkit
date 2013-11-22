#ifndef ITAP_EDGE_H
#define ITAP_EDGE_H

#include <stdlib.h>

#include <string.h>

#include <gmsh/GEdge.h>
#include <gmsh/GFace.h>

#include <iGeom.h>

#include "ITAPVertex.h"

#include "SimpleArray.hpp"
#include "SearchUV.h"

class ITAPEdge : public GEdge
{
  public:
    ITAPEdge(GModel *model, iGeom_Instance &g, iBase_EntityHandle *eHandle, 
             GVertex *v1, GVertex *v2);

    virtual ~ITAPEdge() { }

    virtual Range<double> parBounds(int i) const;
    virtual GPoint point(double p) const;
    virtual SVector3 firstDer(double par) const;

    void * getNativePtr() const { return edgeHandle; }

    virtual int minimumMeshSegments() const;

    virtual SPoint2 reparamOnFace(const GFace *face, double epar, int dir) const;

    double  get_approx_length () const;

    bool  isSeam( const GFace *f) const;

 private:
    void discretize( int numSegments );

    static iBase_TagHandle geom_id_tag;
    static bool tag_available;

    void setGeomTag();
    bool  close_curve;

    iBase_EntityHandle *edgeHandle;
    iGeom_Instance geometry;
    double umin, umax;

    void discretize_open_edge( int numEdges );
    void discretize_close_edge( int numEdges );
};
void save( GEdge *e);

#endif





