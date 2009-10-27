#include "ITAPEdge.h"
#include <sstream>
#include <fstream>
#include <iomanip>

///////////////////////////////////////////////////////////////////////////////

iBase_TagHandle ITAPEdge:: geom_id_tag = 0; 
bool ITAPEdge:: tag_available = 0; 

///////////////////////////////////////////////////////////////////////////////

void ITAPEdge::setGeomTag()
{
   if( !tag_available )  {
      int err;
      const char *tag = "GLOBAL_ID";
      int namelen = strlen(tag);
      iGeom_getTagHandle(geometry, tag, &geom_id_tag, &err, namelen);
      assert( !err );
      tag_available = 1;
    }
}

///////////////////////////////////////////////////////////////////////////////

ITAPEdge::ITAPEdge( GModel *model, iGeom_Instance &g, 
                    iBase_EntityHandle *eHandle, GVertex *v1, GVertex *v2)
  : GEdge(model, 0, v1, v2), geometry(g), edgeHandle(eHandle)
{
    int err;
    int edgeID;
    setGeomTag();
    iGeom_getIntData(geometry, *edgeHandle, geom_id_tag, &edgeID, &err);
    assert( !err );

    setTag( edgeID );

    iGeom_getEntURange(geometry, *edgeHandle, &umin, &umax, &err);

    assert( !err );
}

///////////////////////////////////////////////////////////////////////////////

Range<double> ITAPEdge::parBounds(int i) const
{ 
  return Range<double>(umin, umax);
}

///////////////////////////////////////////////////////////////////////////////

GPoint ITAPEdge::point(double u) const
{
    if( u < umin || u > umax ) { 
       cout << "Warning: U value " << u << " Out of Range : " << umin << " " << umax << endl;
    }

    int err;
    double x, y, z;
    iGeom_getEntUtoXYZ(geometry, *edgeHandle, u, &x, &y, &z, &err);

    assert( !err );
    return GPoint(x,y,z);
}

///////////////////////////////////////////////////////////////////////////////

SPoint2 ITAPEdge::reparamOnFace(const GFace *face, double epar, int dir) const
{
  double u, v, u1, v1, uC, vC;

  u = epar;

  GPoint ePoint = point(u);
  SPoint3 xyz1(ePoint.x(), ePoint.y(), ePoint.z());
  SPoint2 uv = face->parFromPoint(xyz1);

  u = uv.x();
  v = uv.y();

  if( !periodic(0) ) {
      uC = u;
      vC = v;
      if( face->periodic(0) && dir == -1 ) {
          Range<double> ur = face->parBounds(0);
          u1 = uC;
          if( u1 == ur.low()  ) uC = ur.high();
          if( u1 == ur.high() ) uC = ur.low();
      }

      if( face->periodic(1) && dir == -1 ) {
          Range<double> vr = face->parBounds(1);
          v1 = vC;
          if( v1 == vr.low()  ) vC = vr.high();
          if( v1 == vr.high() ) vC = vr.low();
     }

     return SPoint2(uC,vC);
  }

  if( epar > umin &&  epar < umax) return SPoint2(u,v);

  uC = u;
  vC = v;

  double du, dv;

  if( face->periodic(0) ) {
      Range<double> ur = face->parBounds(0);
  if( epar == umin) {
      u1 = uC;
      v1 = vC;
      if( uC == ur.low() )  u1 = ur.high();
      if( uC == ur.high() ) u1 = ur.low();
      double delu = 1.0E-03*(umax - umin);

      SPoint2 uvR = reparamOnFace( face, umin + delu, dir);
      du = uvR.x() - uC;
      dv = uvR.y() - vC;
      double dR   =  sqrt(du*du + dv*dv );

      du = uvR.x() - u1;
      dv = uvR.y() - v1;
      double dL   =  sqrt(du*du + dv*dv );
      if( dL < dR )  uC = u1;
  }

  if( epar == umax) {
      SPoint2 uvR = reparamOnFace( face, umin, dir);

      if( uvR.x() == ur.low()  ) uC = ur.high();
      if( uvR.x() == ur.high() ) uC = ur.low();
  }

  if( dir == -1 ) {
      u1 = uC;
      if( u1 == ur.low()  ) uC = ur.high();
      if( u1 == ur.high() ) uC = ur.low();
  }
  }

  if( face->periodic(1) ) {
      Range<double> vr = face->parBounds(1);

  if( epar == umin ) {
      u1 = uC;
      v1 = vC;
      if( vC == vr.low() )  v1 = vr.high();
      if( vC == vr.high() ) v1 = vr.low();

      double delu = 1.0E-03*(umax - umin);

      SPoint2 uvR = reparamOnFace( face, umin + delu, dir);
      du = uvR.x() - uC;
      dv = uvR.y() - vC;
      double dR   =  sqrt(du*du + dv*dv );

      du = uvR.x() - u1;
      dv = uvR.y() - v1;
      double dL   =  sqrt(du*du + dv*dv );
      if( dL < dR )  vC = v1;
  }

  if( epar == umax) {
      SPoint2 uvR = reparamOnFace( face, umin, dir);

      if( uvR.x() == vr.low()  ) vC = vr.high();
      if( uvR.x() == vr.high() ) vC = vr.low();
  }

  if( dir == -1 ) {
      v1 = vC;
      if( v1 == vr.low()  ) vC = vr.high();
      if( v1 == vr.high() ) vC = vr.low();
  }
  }

  return SPoint2(uC,vC);
}
///////////////////////////////////////////////////////////////////////////////

SVector3 ITAPEdge::firstDer(double u) const
{  
 if( u < umin || u > umax ) 
    cout << "Warning: U value " << u << " Out of Range : " << umin << " " << umax << endl;

  int err;
  double nx, ny, nz;
  double du = 1.0E-06*(umax-umin);
  double uc, xc, yc, zc;  // Center
  double uf, xf, yf, zf;  // Forward 
  double ub, xb, yb, zb;  // backward

  uc = u;
  iGeom_getEntUtoXYZ(geometry, *edgeHandle, uc, &xc, &yc, &zc, &err);

  if( u == umin ) {
      uf = u + du;
      iGeom_getEntUtoXYZ(geometry, *edgeHandle, uf, &xf, &yf, &zf, &err);
      nx = (xf-xc)/(uf-uc);
      ny = (yf-yc)/(uf-uc);
      nz = (zf-zc)/(uf-uc);
      return SVector3( nx, ny, nz);
  }

  if( u == umax ) {
      ub = u - du;
      iGeom_getEntUtoXYZ(geometry, *edgeHandle, ub, &xb, &yb, &zb, &err);
      nx = (xb-xc)/(ub-uc);
      ny = (yb-yc)/(ub-uc);
      nz = (zb-zc)/(ub-uc);
      return SVector3( nx, ny, nz);
  }

  uf = u + du;
  iGeom_getEntUtoXYZ(geometry, *edgeHandle, uf, &xf, &yf, &zf, &err);

  ub = u - du;
  iGeom_getEntUtoXYZ(geometry, *edgeHandle, ub, &xb, &yb, &zb, &err);

  nx = (xf-xb)/(ub-uf);
  ny = (yf-yb)/(ub-uf);
  nz = (zf-zb)/(ub-uf);

  return SVector3( nx, ny, nz);
}

///////////////////////////////////////////////////////////////////////////////
//
bool ITAPEdge::isSeam (const GFace *face) const
{
  int err, sense_out;
   
  const iBase_EntityHandle *faceHandle = 
                   (iBase_EntityHandle *)face->getNativePtr();

  iGeom_getEgFcSense( geometry, *edgeHandle, *faceHandle, &sense_out, &err);

  if( sense_out == 0) return 1;
  
  return 0;
}
///////////////////////////////////////////////////////////////////////////////

int ITAPEdge::minimumMeshSegments() const
{
  int np = GEdge::minimumMeshSegments();

/*
   int np = GEdge::minimumMeshSegments();
  if(geomType() == Line)
    np = GEdge::minimumMeshSegments();
  else 
    np = CTX::instance()->mesh.minCurvPoints - 1;
*/

  // if the edge is closed, ensure that at least 3 points are
  // generated in the 1D mesh (4 segments, one of which is
  // degenerated)
  if (getBeginVertex() == getEndVertex()) np = std::max(4, np);

  return std::max(np, meshAttributes.minimumMeshSegments);
}

///////////////////////////////////////////////////////////////////////////////
double ITAPEdge::get_approx_length() const
{
    int err;
    SimpleArray<double> measure;
    iGeom_measure( geometry, edgeHandle, 1, ARRAY_INOUT( measure), &err);
    return measure[0];
}
///////////////////////////////////////////////////////////////////////////////
void ITAPEdge :: discretize_close_edge(int numEdges)
{
    close_curve = 1;

    int err, status;
    double x, y, z;

    int numNodes = numEdges;

    double u, du = (umax - umin) / (double) numEdges;

    mesh_vertices.resize(numNodes -1);
    vector<MVertex*>  nodes( numNodes );

    nodes[0] =  this->getBeginVertex()->mesh_vertices[0];

    for (int i = 1; i < numNodes; i++)
    {
        u = umin + i*du;
        iGeom_getEntUtoXYZ(geometry, *edgeHandle, u, &x, &y, &z, &err);
        MVertex *v = new MEdgeVertex(x, y, z, this, u);
        mesh_vertices[i-1] = v;
        nodes[i] = v;
    }

    lines.resize( numEdges);
    for (int i = 0; i < numEdges; i++)
    {
        MVertex *v0  = nodes[i];
        MVertex *v1  = nodes[(i + 1) % numNodes];
        lines[i] = new MLine(v0,v1);
    }
}

////////////////////////////////////////////////////////////////////////////////

void ITAPEdge :: discretize_open_edge( int numEdges )
{
    close_curve = 0;
    int err, status, numNodes, index;
    double x, y, z, u, du;

    du = (umax - umin) / (double) numEdges;

    numNodes = numEdges + 1;

    vector<MVertex*>  nodes( numNodes );

    mesh_vertices.resize( numNodes - 2);
     
    for (int i = 1; i < numNodes - 1; i++)
    {
        u = umin + i*du;
        iGeom_getEntUtoXYZ(geometry, *edgeHandle, u, &x, &y, &z, &err);
        MVertex *v = new MEdgeVertex( x, y, z, this, u);
        mesh_vertices[i-1] = v;
        nodes[i] = v;
    }

    nodes[0] =  this->getBeginVertex()->mesh_vertices[0];
    assert( nodes[0] );
    nodes[numNodes-1] =  this->getEndVertex()->mesh_vertices[0];
    assert( nodes[numNodes-1] );

    lines.resize(numEdges);
    for (int i = 0; i < numEdges; i++)
    {
        MVertex *v0 = nodes[i]; assert(v0);
        MVertex *v1 = nodes[i+1]; assert(v1);
        lines[i] = new MLine(v0,v1);
    }
}

////////////////////////////////////////////////////////////////////////////////
void ITAPEdge :: discretize(int numEdges)
{
    int err;

    SimpleArray<iBase_EntityHandle> gNodes;
    iGeom_getEntAdj(geometry, *edgeHandle, iBase_VERTEX, ARRAY_INOUT(gNodes), &err);

    switch (gNodes.size())
    {
    case 1:
        discretize_close_edge(numEdges);
        break;
    case 2:
        discretize_open_edge(numEdges);
        break;
    default:
        cout << "Fatal Error: Invalid Geometric edge " << endl;
        exit(0);
    }
}
////////////////////////////////////////////////////////////////////////////////
void save( GEdge *edge)
{
    ostringstream oss;
    oss << "./ModelMesh/itapedge" << edge->tag() << ".off";
    ofstream ofile( oss.str().c_str(), ios::out);
    if( ofile.fail() ) {
        cout << " Warning: Cann't open file " << oss.str() << endl;
        return;
    }

    int numSegments = edge->lines.size();

    std::map<int,MVertex*> vmap;
    std::map<int,MVertex*>::iterator it; 

    for( int i = 0; i < numSegments; i++) {
         MVertex *v0 = edge->lines[i]->getVertex(0); assert(v0);
         MVertex *v1 = edge->lines[i]->getVertex(1); assert(v1);
         vmap[v0->getNum()] = v0;
         vmap[v1->getNum()] = v1;
   }

   ofile << "OFF" << endl;
   ofile << vmap.size() << " 0 " << numSegments << endl;

   map<int, int> idmap;
   int index = 0;
   for( it = vmap.begin(); it != vmap.end(); ++it) 
   {
        MVertex *v = it->second;
        ofile << fixed << v->x() << " " << v->y() << " " << v->z() << endl;
        idmap[it->first] = index++;
  }

  for( int i = 0; i < numSegments; i++) 
  {
         MVertex *v0 = edge->lines[i]->getVertex(0);
         MVertex *v1 = edge->lines[i]->getVertex(1);
         ofile << idmap[v0->getNum()] << "  " << idmap[v1->getNum()] << endl;
  }
}
////////////////////////////////////////////////////////////////////////////////
