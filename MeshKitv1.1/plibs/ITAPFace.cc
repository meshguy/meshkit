#include "ITAPFace.h"
#include <sstream>
#include <fstream>

iBase_TagHandle ITAPFace::geom_id_tag = 0;
bool ITAPFace::tag_available = 0;

//////////////////////////////////////////////////////////////////////////////////

void ITAPFace::setGeomTag()
{
    if (!tag_available)
    {
        int err;
        const char *tag = "GLOBAL_ID";
        int namelen = strlen(tag);
        iGeom_getTagHandle(geometry, tag, &geom_id_tag, &err, namelen);
        tag_available = 1;
    }
}

//////////////////////////////////////////////////////////////////////////////////

ITAPFace::ITAPFace(GModel *model, iGeom_Instance &g, iBase_EntityHandle *fHandle)
: GFace(model, 0), geometry(g), faceHandle(fHandle)
{
    int err;

    int faceID;
    setGeomTag();
    iGeom_getIntData(geometry, *faceHandle, geom_id_tag, &faceID, &err);
    assert(!err);

    setTag(faceID);

    iGeom_getEntUVRange(geometry, *faceHandle, &umin, &vmin, &umax, &vmax, &err);
    assert( !err );

    iGeom_getEntBoundBox( geometry, *faceHandle, &xmin, &ymin, &zmin, 
                          &xmax, &ymax, &zmax, &err);

    _periodic[0] = 0;
    _periodic[1] = 0;
    int in_u, in_v;
    iGeom_isEntPeriodic( geometry, *faceHandle, &in_u, &in_v, &err);
    if( in_u ) _periodic[0] = 1;
    if( in_v ) _periodic[1] = 1;

    xlength = fabs(xmax-xmin);
    ylength = fabs(ymax-ymin);
    zlength = fabs(zmax-zmin);
    maxlength = max( xlength, max( ylength, zlength) );

    SimpleArray<iBase_EntityHandle> edgeHandles;
    iGeom_getEntAdj(geometry, *faceHandle, iBase_EDGE, ARRAY_INOUT(edgeHandles), &err);
    assert(!err);

    int edgeID;
    std::list<GEdge*> l_wire, closed_loop;

    for (int i = 0; i < edgeHandles.size(); i++)
    {
        iGeom_getIntData(geometry, edgeHandles[i], geom_id_tag, &edgeID, &err);
        GEdge *edge = model->getEdgeByTag(edgeID); assert(edge);
        l_wire.push_back(edge);
        edge->addFace(this);
    }

    vector<list<GEdge*> > icclist = get_icc_list( l_wire );
    for( int i = 0; i < icclist.size(); i++) 
         addEdgeLoop( icclist[i] );

   kdtree   = NULL;
   uvCoords = NULL;
   numNeighs = 1;
   annIdx.resize(numNeighs);
   anndist.resize(numNeighs);
}

///////////////////////////////////////////////////////////////////////////////

bool ITAPFace :: hasSeam() const 
{
    int err, sense_out;

    SimpleArray<iBase_EntityHandle> edgeHandles;
    iGeom_getEntAdj(geometry, *faceHandle, iBase_EDGE, ARRAY_INOUT(edgeHandles), &err);
    assert(!err);

    for( int i = 0; i < edgeHandles.size(); i++) {
         iGeom_getEgFcSense( geometry, edgeHandles[i], *faceHandle, &sense_out, &err);
         if( sense_out == 0) return 1;
    }

    return 0;
}

///////////////////////////////////////////////////////////////////////////////

void  ITAPFace :: create_kdtree()
{
    delete_kdtree(); // If Allocated earlier

    int err;

    int N = 1000;
    double x, y, z;

    double du = (umax-umin)/(double)N;
    double dv = (vmax-vmin)/(double)N;

    int numnodes = (N + 1)*(N+1);
    kdnodes = annAllocPts( numnodes, 3);

    uvCoords = new double[2*numnodes];
 
    int index = 0;
    for( int j = 0; j < N+1; j++) {
         double v = vmin + j*dv; if( v > vmax ) v = vmax;
         for( int i = 0; i < N+1; i++) {
              double u = umin + i*du; if( u > umax ) u = umax;
              iGeom_getEntUVtoXYZ(geometry, *faceHandle, u, v, &x, &y, &z, &err);
              kdnodes[index][0]   = x;
              kdnodes[index][1]   = y;
              kdnodes[index][2]   = z;
              uvCoords[2*index+0] = u;
              uvCoords[2*index+1] = v;
              index++;
          }
    }

    kdtree = new ANNkd_tree( kdnodes, numnodes,3);
    assert( kdtree );
}

///////////////////////////////////////////////////////////////////////////////

void  ITAPFace :: delete_kdtree()
{
   if( kdtree == NULL ) return;

   annDeallocPts( kdnodes );
   delete [] uvCoords;
   delete kdtree; kdtree = NULL;
   annClose();
}

///////////////////////////////////////////////////////////////////////////////

vector<std::list<GEdge*> > ITAPFace :: get_icc_list( std::list<GEdge*> &l)
{
   GEdge *edge, *start_edge;

   std::list<GEdge*> newlist;
   vector<std::list<GEdge*> > icclist;

   std::set<GEdge*> eset;
   BOOST_FOREACH( edge, l) eset.insert( edge );

   std::set<GVertex*> vset;
   while( !eset.empty() ) 
   {
       start_edge = *eset.begin();  
       BOOST_FOREACH( edge, eset) {
           if( edge->isSeam( this) ) start_edge = edge;
       }
       eset.erase( start_edge );
       newlist.clear();
       newlist.push_back(start_edge);
       if(start_edge->isSeam( this )) newlist.push_back(start_edge);
       vset.clear();
       vset.insert( start_edge->getBeginVertex() );
       vset.insert( start_edge->getEndVertex() );
       while(1) {
            int success = 0;
            BOOST_FOREACH(edge, eset) 
            {
                GVertex *v1 = edge->getBeginVertex();
                GVertex *v2 = edge->getEndVertex();
                if( vset.find(v1) != vset.end() ) {
                    vset.insert(v2);
                    eset.erase(edge);
                    newlist.push_back(edge);
                    if(edge->isSeam( this )) newlist.push_back(edge);
                    success = 1;
                    break;
                }
                if( vset.find(v2) != vset.end() ) {
                    vset.insert(v1);
                    eset.erase(edge);
                    newlist.push_back(edge);
                    if(edge->isSeam( this )) newlist.push_back(edge);
                    success = 1;
                    break;
                }
            }
            if( !success ) break;
        }
        icclist.push_back( newlist );
   }
   return icclist;
}

///////////////////////////////////////////////////////////////////////////////

void ITAPFace :: addEdgeLoop( std::list<GEdge*> &wires)
{
    if( wires.empty()  ) return;

    GEdgeLoop el(wires);
    for (GEdgeLoop::citer it = el.begin(); it != el.end(); ++it)
    {
        l_edges.push_back(it->ge);
        l_dirs.push_back(it->_sign);
        if (el.count() == 2)
        {
            it->ge->meshAttributes.minimumMeshSegments =
                    std::max(it->ge->meshAttributes.minimumMeshSegments, 2);
        }
        if (el.count() == 1)
        {
            it->ge->meshAttributes.minimumMeshSegments =
                    std::max(it->ge->meshAttributes.minimumMeshSegments, 3);
        }
    }
    edgeLoops.push_back(el);
}

///////////////////////////////////////////////////////////////////////////////

Range<double> ITAPFace::parBounds(int i) const
{
    if (i == 0) return Range<double>(umin, umax);

    return Range<double>(vmin, vmax);
}

///////////////////////////////////////////////////////////////////////////////

SVector3 ITAPFace::normal(const SPoint2 &param) const
{
    int err;
    double nx, ny, nz;
    iGeom_getEntNrmlUV(geometry, *faceHandle, param.x(), param.y(),
                       &nx, &ny, &nz, &err); assert( !err );

    SVector3 n(nx, ny, nz);
    n.normalize();

    return n;
}
///////////////////////////////////////////////////////////////////////////////

Pair<SVector3, SVector3> ITAPFace::firstDer(const SPoint2 &param) const
{
    int err;
    static SimpleArray<double> uDeriv, vDeriv;
    iGeom_getEnt1stDrvt(geometry, *faceHandle, param.x(), param.y(),
                        ARRAY_INOUT(uDeriv), ARRAY_INOUT(vDeriv),
                        &err); assert( !err );

    SVector3 uVec(uDeriv[0], uDeriv[1], uDeriv[2]);
    SVector3 vVec(vDeriv[0], vDeriv[1], vDeriv[2]);

    return Pair<SVector3, SVector3 > (uVec, vVec);
}

///////////////////////////////////////////////////////////////////////////////

GPoint ITAPFace::point(double u, double v) const
{
    int err;
    double pp[2] = {u, v};

    double x, y, z;
    iGeom_getEntUVtoXYZ(geometry, *faceHandle, u, v, &x, &y, &z, &err);
    assert( !err );

    return GPoint(x, y, z, this, pp);
}

///////////////////////////////////////////////////////////////////////////////
SPoint2  ITAPFace :: parFromPoint( const SPoint3 &p) const
{
    double tol = 1.0E-06;
    int err;
    double x, y, z, u, v, dx, dy, dz, derr;

    x = p.x();
    y = p.y();
    z = p.z();

    if( x < xmin || x > xmax ) 
        cout << "Warning: Query point outside X Range " << endl;

    if( y < ymin || y > ymax ) 
        cout << "Warning: Query point outside Y Range " << endl;

    if( z < zmin || z > zmax ) 
        cout << "Warning: Query point outside Z Range " << endl;

    iGeom_getEntXYZtoUV(geometry, *faceHandle, x, y, z, &u, &v, &err);
    assert( !err );

    iGeom_getEntUVtoXYZ(geometry, *faceHandle, u, v, &x, &y, &z, &err);
    assert( !err );

    dx = fabs( x - p.x() );
    dy = fabs( y - p.y() );
    dz = fabs( z - p.z() );
    derr = dx*dx + dy*dy + dz*dz;

    if( derr < tol*tol) return SPoint2(u,v);

    double queryPoint[3], eps = 0.0;
    double xon, yon, zon;
    double dist1, dist2; 
 
    queryPoint[0] = p.x();
    queryPoint[1] = p.y();
    queryPoint[2] = p.z();

    kdtree->annkSearch(queryPoint, numNeighs, &annIdx[0], &anndist[0], eps);

    int index = annIdx[0];
    dist1  = anndist[0];

    x = kdnodes[index][0];
    y = kdnodes[index][1];
    z = kdnodes[index][2];
    u = uvCoords[2*index+0];
    v = uvCoords[2*index+1];

    /*
    double uguess = u;
    double vguess = v;
    
    iGeom_getEntXYZtoUVHint( geometry, *faceHandle, x, y, z, &uguess, &vguess, &err);
    iGeom_getEntUVtoXYZ(geometry, *faceHandle, uguess, vguess, &xon, &yon, &zon, &err);
    dx = queryPoint[0] - xon;
    dy = queryPoint[1] - yon;
    dz = queryPoint[2] - zon;
    dist2 =  dx*dx + dy*dy + dz*dz;

    if( dist2 < dist1 ) {
        u = uguess;
        v = vguess;
    }
    */

    return SPoint2(u,v);
}
///////////////////////////////////////////////////////////////////////////////

GPoint ITAPFace::closestPoint(const SPoint3 & qp, const double initialGuess[2]) const
{
    int err;
    double queryPoint[3], eps = 0.0;
    double xon, yon, zon;
    double dx, dy, dz, dist1, dist2; 
 
    queryPoint[0] = qp.x();
    queryPoint[1] = qp.y();
    queryPoint[2] = qp.z();

    iGeom_getEntClosestPt( geometry, *faceHandle, qp.x(), qp.y(), qp.z(),
                           &xon, &yon, &zon, &err);
    assert( !err );
    dx = queryPoint[0] - xon;
    dy = queryPoint[1] - yon;
    dz = queryPoint[2] - zon;
    dist2 =  dx*dx + dy*dy + dz*dz;

    if( kdtree ) 
    {
        kdtree->annkSearch(queryPoint, numNeighs, &annIdx[0], &anndist[0], eps);
        dist1 = anndist[0];
    }

    if( dist2 > dist1) 
        cout << " IGeom Failed : " << dist2 << "  KDTree " << dist1 << endl;
    else
        cout << " IGeom Pass   : " << dist2 << "  KDTree " << dist1 << endl;

    int index = annIdx[0];

    double u = uvCoords[2*index+0];
    double v = uvCoords[2*index+1];

    iGeom_getEntUVtoXYZ(geometry, *faceHandle, u, v, &xon, &yon, &zon, &err);

    dx = queryPoint[0] - xon;
    dy = queryPoint[1] - yon;
    dz = queryPoint[2] - zon;
    dist2 =  dx*dx + dy*dy + dz*dz;

    cout << " Minimum Distance " <<  fixed << dist1 << "  " << dist2 << endl;
   
    return GPoint(xon, yon, zon);
}

///////////////////////////////////////////////////////////////////////////////

void save( GFace *face) 
{
    ostringstream oss;
    oss << "./ModelMesh/itapface" << face->tag() << ".off";
    ofstream ofile(oss.str().c_str(), ios::out);
    if (ofile.fail()) {
        cout << " Warning: Cann't open file " << oss.str() << endl;
        return;
    }

    std::set<MVertex*> vset;
    std::set<MVertex*>::iterator it;

    unsigned elemtype[2];
    elemtype[0] = 0;
    elemtype[1] = 0;
    face->getNumMeshElements( elemtype );
   
    int numTriangles = elemtype[0];
    int numQuads     = elemtype[1];

    int numFaces = face->getNumMeshElements();
    assert( numFaces == numTriangles + numQuads );

    list<GEdge*> ledges = face->edges();
    list<GEdge*>::const_iterator eit;

    int total_edges = 0;
    for (eit = ledges.begin(); eit != ledges.end(); ++eit)
    {
        GEdge *edge = *eit;
        int numSegments = edge->lines.size();
        for (int i = 0; i < numSegments; i++)
        {
            MVertex *v0 = edge->lines[i]->getVertex(0);
            assert(v0);
            MVertex *v1 = edge->lines[i]->getVertex(1);
            assert(v1);
            vset.insert(v0);
            vset.insert(v1);
        }
        total_edges += numSegments;
    }

    ofile << "OFF" << endl;
    ofile << vset.size() << " 0 " << total_edges << endl;

    map<MVertex*, int> idmap;
    int index = 0;
    for (it = vset.begin(); it != vset.end(); ++it)
    {
        MVertex *vertex = *it;
        ofile << vertex->x() << " " << vertex->y() << " " << vertex->z() << endl;
        if( idmap.find(vertex) == idmap.end() ) 
            idmap[vertex] = index++;
    }

    for (eit = ledges.begin(); eit != ledges.end(); ++eit)
    {
        GEdge *edge = *eit;
        int numSegments = edge->lines.size();
        for (int i = 0; i < numSegments; i++)
        {
            MVertex *v0 = edge->lines[i]->getVertex(0);
            MVertex *v1 = edge->lines[i]->getVertex(1);
            ofile << idmap[v0] << "  " << idmap[v1] << endl;
        }
        ofile << endl;
    }

}
////////////////////////////////////////////////////////////////////////////////

