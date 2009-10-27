#include "ITAPVertex.h"

iBase_TagHandle ITAPVertex:: geom_id_tag = 0; 
bool ITAPVertex:: tag_available = 0; 

///////////////////////////////////////////////////////////////////////////////

void ITAPVertex::setGeomTag()
{
    if( !tag_available )  {
        int err;
        const char *tag = "GLOBAL_ID";
        int namelen = strlen(tag);
        iGeom_getTagHandle(geometry, tag, &geom_id_tag, &err, namelen);
        tag_available = 1;
    }
}

///////////////////////////////////////////////////////////////////////////////

ITAPVertex::ITAPVertex(GModel *model, iGeom_Instance &g, iBase_EntityHandle *v)
: GVertex(model, 0), geometry(g), vertexHandle(v)
{
    int err;
    int vertexID;

    setGeomTag();
    iGeom_getIntData(geometry, *vertexHandle, geom_id_tag, &vertexID, &err);
    setTag( vertexID );
   
    double x, y, z;
    iGeom_getVtxCoord(geometry, *vertexHandle, &x, &y, &z, &err);

    xyzCoord[0] = x;
    xyzCoord[1] = y;
    xyzCoord[2] = z;

    mesh_vertices.push_back(new MVertex(x, y, z, this));
    points.push_back(new MPoint(mesh_vertices.back()));
}

///////////////////////////////////////////////////////////////////////////////

void ITAPVertex::setPosition(GPoint &p)
{
    xyzCoord[0] = p.x();
    xyzCoord[1] = p.y();
    xyzCoord[2] = p.z();

    if (mesh_vertices.size())
    {
        mesh_vertices[0]->x() = p.x();
        mesh_vertices[0]->y() = p.y();
        mesh_vertices[0]->z() = p.z();
    }
}

///////////////////////////////////////////////////////////////////////////////
SPoint2 ITAPVertex::reparamOnFace(const GFace *gface, int dir) const
{
  std::list<GEdge*>::const_iterator it = l_edges.begin();

  int err;
  double s1,s0;

  std::list<GEdge*> l = gface->edges();

  while(it != l_edges.end())
  {
    GEdge *curredge = *it;
    if(std::find(l.begin(), l.end(), curredge) != l.end()){
       if((curredge)->isSeam(gface))
       {
           iBase_EntityHandle *edgeHandle = (iBase_EntityHandle*)(curredge)->getNativePtr();
           iGeom_getEntURange( geometry, *edgeHandle, &s0, &s1, &err);
           if(curredge->getBeginVertex() == this)
              return curredge->reparamOnFace(gface, s0, dir);
           else if(curredge->getEndVertex() == this)
              return curredge->reparamOnFace(gface, s1, dir);
      }
    }
    ++it;
  }  

  it = l_edges.begin();
  while(it != l_edges.end()){
    GEdge *curredge = *it;
    if(std::find(l.begin(), l.end(), curredge) != l.end())
    {
      iBase_EntityHandle *edgeHandle = (iBase_EntityHandle*)(*it)->getNativePtr();
      iGeom_getEntURange( geometry, *edgeHandle, &s0, &s1, &err);
      if(curredge->getBeginVertex() == this)
        return curredge->reparamOnFace(gface, s0, dir);
      else if( curredge->getEndVertex() == this)
        return curredge->reparamOnFace(gface, s1, dir);
    }
    ++it;
  }

  // normally never here
  return GVertex::reparamOnFace(gface, dir);
}

