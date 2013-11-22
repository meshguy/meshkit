#ifndef ITAP_VERTEX_H
#define ITAP_VERTEX_H

#include <string.h>
#include <iostream>

using namespace std;

#include <gmsh/GVertex.h>
#include <gmsh/GEdge.h>
#include <gmsh/GFace.h>
#include <gmsh/MElement.h>

#include <iGeom.h>

class ITAPVertex : public GVertex
{
public:
    ITAPVertex(GModel *m, iGeom_Instance &g, iBase_EntityHandle *h);

    virtual ~ITAPVertex() { }

    virtual GPoint point() const
    {
        return GPoint(xyzCoord[0], xyzCoord[1], xyzCoord[2]);
    }

    virtual double x() const { return xyzCoord[0]; }
    virtual double y() const { return xyzCoord[1]; }
    virtual double z() const { return xyzCoord[2]; }

    virtual void setPosition(GPoint &p);

    void *getNativePtr() const { return vertexHandle; }

    SPoint2 reparamOnFace(const GFace *gf, int dir) const;

protected:
    static iBase_TagHandle geom_id_tag;
    static bool tag_available;

    void setGeomTag();

    iBase_EntityHandle *vertexHandle;
    double xyzCoord[3];
    iGeom_Instance geometry;
};

#endif
