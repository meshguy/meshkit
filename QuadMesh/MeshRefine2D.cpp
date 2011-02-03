#include "MeshRefine2D.h"

using namespace Jaal;

////////////////////////////////////////////////////////////////////////////////

int MeshRefine2D :: initialize()
{
  insertedNodes.clear();
  insertedFaces.clear();
  return 0;
}

////////////////////////////////////////////////////////////////////////////////

int MeshRefine2D :: finalize() 
{
  mesh->prune();
  return 0;
}

///////////////////////////////////////////////////////////////////////////////

void 
MeshRefine2D::RefinedEdgeMap::clear()
{
  std::map<Vertex*, vector<RefinedEdge> >::const_iterator it;

  for( it = refined_edges.begin(); it != refined_edges.end(); ++it) 
  {
   const vector<RefinedEdge> &refedges = it->second;
   for( size_t i = 0; i < refedges.size(); i++) 
        delete refedges[i].edge;
  }
  refined_edges.clear();
}

///////////////////////////////////////////////////////////////////////////////
bool
MeshRefine2D::RefinedEdgeMap::hasEdge( Vertex *v1, Vertex *v2) const
{
   Vertex *vmin = std::min(v1,v2);

   std::map<Vertex*, vector<RefinedEdge> >::const_iterator it;
   it = refined_edges.find(vmin);

   if( it == refined_edges.end() ) return 0;

   const vector<RefinedEdge> &refedges = it->second;
   for( int i = 0; i < refedges.size(); i++) {
        Vertex *ev1 = refedges[i].edge->getNodeAt(0);
        Vertex *ev2 = refedges[i].edge->getNodeAt(1);
	if( ev1 == v1 || ev2 == v2) return 1;
	if( ev1 == v2 || ev2 == v1) return 1;
   }
   return 0;
}

///////////////////////////////////////////////////////////////////////////////

bool 
MeshRefine2D::RefinedEdgeMap::allow_edge_refinement( const Edge *edge) const
{
   if( edge->isBoundary() || edge->isConstrained() ) 
        if( boundary_split_flag == 0) return 0;

   return 1;
}

///////////////////////////////////////////////////////////////////////////////

int 
MeshRefine2D::RefinedEdgeMap::setVertexOnEdge( Vertex *v1, Vertex *v2) 
{
  if( hasEdge(v1,v2) ) return 0;

  RefinedEdge  refedge;
  Edge *edge = new Edge(v1,v2);
  if( allow_edge_refinement(edge) ) 
  {
      Point3D p3d = Vertex::mid_point(v1,v2);
      Vertex *v   = Vertex::newObject();
      v->setXYZCoords( p3d );
      refedge.edge = edge;
      refedge.midVertex = v;
      Vertex *vmin = std::min(v1,v2);
      refined_edges[vmin].push_back( refedge );
      cout << refined_edges.size() << endl;
      return 0;
   }

   return 1;
}

///////////////////////////////////////////////////////////////////////////////

Vertex* MeshRefine2D::RefinedEdgeMap::getVertexOnEdge( Vertex *v1, Vertex *v2 ) const
{
   Vertex *vmin = std::min(v1,v2);

   std::map<Vertex*, vector<RefinedEdge> >::const_iterator it;
   it = refined_edges.find(vmin);

   if( it == refined_edges.end() ) return NULL;
   const vector<RefinedEdge> &refedges = it->second;
   for( int i = 0; i < refedges.size(); i++) {
        Vertex *ev1 = refedges[i].edge->getNodeAt(0);
        Vertex *ev2 = refedges[i].edge->getNodeAt(1);
	if( ev1 == v1 || ev2 == v2) return refedges[i].midVertex;
	if( ev1 == v2 || ev2 == v1) return refedges[i].midVertex;
   }

   return NULL;
}

///////////////////////////////////////////////////////////////////////////////

void MeshRefine2D::append_new_node( Vertex *vertex) 
{
  mesh->addNode( vertex );
  insertedNodes.push_back( vertex );
}

///////////////////////////////////////////////////////////////////////////////

Face* MeshRefine2D::append_new_triangle( Vertex *v0, Vertex *v1 , Vertex *v2) 
{
  NodeSequence pnodes(3);
  pnodes[0] = v0;
  pnodes[1] = v1;
  pnodes[2] = v2;

  Face *face = Face::newObject();
  face->setNodes( pnodes );
  mesh->addFace(face);
  insertedFaces.push_back( face );
  return face;
}

///////////////////////////////////////////////////////////////////////////////

Face* MeshRefine2D::append_new_quad( Vertex *v0, Vertex *v1 , Vertex *v2, Vertex *v3) 
{
  NodeSequence pnodes(4);
  pnodes[0] = v0;
  pnodes[1] = v1;
  pnodes[2] = v2;
  pnodes[3] = v3;

  Face *face = Face::newObject();
  face->setNodes( pnodes );
  mesh->addFace(face);
  insertedFaces.push_back( face );
  return face;
}

///////////////////////////////////////////////////////////////////////////////

int Sqrt3Refine2D :: execute()
{
   CentroidRefine2D refine(mesh);
   SwapTriEdge eflip(mesh);

   for ( int itime = 0; itime < numIterations; itime++) {
       refine.execute();
       eflip.execute();
   }

   return 0;
}

///////////////////////////////////////////////////////////////////////////////

int LongestEdgeRefine2D :: atomicOp( const Face *oldface)
{
  vector<double> elen(3);

  double maxlen = 0.0;
  for( int i = 0; i < 3; i++) {
      Vertex *v1  = oldface->getNodeAt((i+1)%3);
      Vertex *v2  = oldface->getNodeAt((i+2)%3);
      elen[i] = Vertex::length(v1, v2);
      maxlen  = std::max(elen[i], maxlen);
  }

  for( int i = 0; i < 3; i++) {
      if( elen[i]/maxlen > 0.90 ) 
       {
           Vertex *v1 = oldface->getNodeAt( (i+1)%3);
           Vertex *v2 = oldface->getNodeAt( (i+2)%3);
           edgemap->setVertexOnEdge(v1,v2);
        }
   }

   return 0;
}

///////////////////////////////////////////////////////////////////////////////

int LongestEdgeRefine2D::execute() 
{
  edgemap = new RefinedEdgeMap;

  size_t numfaces = mesh->getSize(2);

  size_t ncount = 0;
  for( size_t i = 0; i < numfaces; i++) {
       Face *face = mesh->getFaceAt(i);
	double ratio  = face->getAspectRatio();
        if( ratio < cutOffAspectRatio) {
	     int err = atomicOp( face );
	     if( !err ) ncount++;
        }
  }

  if( ncount ) {
      ConsistencyRefine2D consistency(mesh, edgemap);
      consistency.execute();
      finalize();
  } else
    cout << "Warning: No Edge was refined " << endl;

  edgemap->clear();

  delete edgemap;


  return 0;
}

///////////////////////////////////////////////////////////////////////////////

int  ConsistencyRefine2D :: execute()
{
  hangingVertex.resize(3);
  edge0.set(0);
  edge1.set(1);
  edge2.set(2);

  makeConsistent();

  finalize();

  return 0;
}

//#############################################################################

void ConsistencyRefine2D :: subDivideQuad2Tri( const NodeSequence &connect)
{
  int err, status;
  assert( connect.size() == 4 );
  //********************************************************************
  // Subdivide a quadrilateral cell into two triangles. We can choose
  // either quadrilateral diagonal for spliiting, but we choose to
  // select quadrilateral which gives, maximum of minimum aspect
  // ratio of the resulting two triangles ....
  //*******************************************************************
  double diagonal[2];

  Vertex *v0 = connect[0]; 
  Vertex *v1 = connect[1]; 
  Vertex *v2 = connect[2]; 
  Vertex *v3 = connect[3]; 

  diagonal[0]  = Vertex::length( v0, v3 );
  diagonal[1]  = Vertex::length( v1, v2 );

  if( diagonal[0] < diagonal[1] ) {
    append_new_triangle( v0, v1, v3);
    append_new_triangle( v0, v3, v2);
  } else {
    append_new_triangle( v0, v1, v2);
    append_new_triangle( v1, v3, v2);
  }
}

//#############################################################################

void ConsistencyRefine2D :: makeConsistent1( Face *oldface)
{
  //------------------------------------------------------------------
  // When only one edge is inconsistent, we can direcly join the
  // hanging node to the opposite node of the triangle. Therefore
  // one additional triangle is generated with this case.
  //------------------------------------------------------------------
  if( oldface->getSize(0) != 3 ) return;

  Vertex *n1 = oldface->getNodeAt(0);
  Vertex *n2 = oldface->getNodeAt(1);
  Vertex *n3 = oldface->getNodeAt(2);

  remove_it( oldface );

  // If the only hanging vertex lies of the edge0. i.e. opposite to the
  // vertex 0 of the existing triangle.
  if( bitvec == edge0) {
    append_new_triangle(n1, n2, hangingVertex[0] );

    // Create a new triangle ....
    append_new_triangle(n3, n1, hangingVertex[0] );
    return;
  }


  // If the only hanging vertex lies of the edge1. i.e. opposite to the
  // vertex 1 of the existing triangle.
  if( bitvec == edge1) 
  {
    // Replace the existing triangle ....
    append_new_triangle(n2, n3, hangingVertex[1] );

    // Create a new triangle ....
    append_new_triangle(n2, hangingVertex[1], n1 );
    return;
  }

  // If the only hanging vertex lies of the edge2. i.e. opposite to the
  // vertex 2 of the existing triangle.
  if( bitvec == edge2) {
    // Replace the existing triangle ....
    append_new_triangle(n3, n1, hangingVertex[2] );

    // Create a new triangle ....
    append_new_triangle(n3, hangingVertex[2], n2 );
    return;
  }

}

//#############################################################################

void ConsistencyRefine2D :: refineEdge0(const Face *oldface)
{
  Vertex *n1 = oldface->getNodeAt(0);
  Vertex *n2 = oldface->getNodeAt(1);
  Vertex *n3 = oldface->getNodeAt(2);

  append_new_triangle(n1, hangingVertex[2], hangingVertex[1] );

  // One Quadrilateral is created, divide into 2 triangles.
  NodeSequence qnodes(4);
  qnodes[0] =  hangingVertex[2];
  qnodes[1] =  n2;
  qnodes[2] =  hangingVertex[1];
  qnodes[3] =  n3;

  subDivideQuad2Tri( qnodes );
}

//#############################################################################

void ConsistencyRefine2D :: refineEdge1(const Face *oldface)
{
  Vertex *n1 = oldface->getNodeAt( 0 );
  Vertex *n2 = oldface->getNodeAt( 1 );
  Vertex *n3 = oldface->getNodeAt( 2 );

  append_new_triangle(n2, hangingVertex[0], hangingVertex[2] );

  // One Quadrilateral is created, divide into 2 triangles.
  NodeSequence qnodes(4);
  qnodes[0] =  n1;
  qnodes[1] =  hangingVertex[2];
  qnodes[2] =  n3;
  qnodes[3] =  hangingVertex[0];

  subDivideQuad2Tri( qnodes );
}


//#############################################################################

void ConsistencyRefine2D :: refineEdge2(const Face *oldface)
{
  Vertex *n1 = oldface->getNodeAt(0);
  Vertex *n2 = oldface->getNodeAt(1);
  Vertex *n3 = oldface->getNodeAt(2);

  append_new_triangle(n3, hangingVertex[1], hangingVertex[0] );

  NodeSequence qnodes(4);
  qnodes[0] =  n1;
  qnodes[1] =  n2;
  qnodes[2] =  hangingVertex[1];
  qnodes[3] =  hangingVertex[0];

  subDivideQuad2Tri( qnodes );
}

//#############################################################################

void ConsistencyRefine2D :: makeConsistent2( Face *oldface)
{
  //--------------------------------------------------------------------
  // When there are two edges which are inconsistent, then we create
  // one triangle and one quadrilateral. This quadrilateral is further
  // divided into 2 triangle, which produces better aspect ratio.
  // Therefore, three triangles are generated in this procedure.
  //--------------------------------------------------------------------
  // Find out which edge is consistent ...
  bitvec.flip();

  if( bitvec == edge0) refineEdge0(oldface); 
  if( bitvec == edge1) refineEdge1(oldface); 
  if( bitvec == edge2) refineEdge2(oldface); 
  
  remove_it( oldface );

}

//#############################################################################

void ConsistencyRefine2D :: makeConsistent3( Face *oldface)
{
  Vertex *n1 = oldface->getNodeAt(0);
  Vertex *n2 = oldface->getNodeAt(1);
  Vertex *n3 = oldface->getNodeAt(2);

  // First Triangle  Using 0(old), 1,2(new): Replace existing triangle
  append_new_triangle( n1, hangingVertex[2], hangingVertex[1] );

  // Second Triangle  Using 1(old), 2,0(new): Create new triangle
  append_new_triangle( n2, hangingVertex[0], hangingVertex[2] );

  // Second Triangle  Using 2(old), 1,0(new): Create new triangle
  append_new_triangle( n3, hangingVertex[1], hangingVertex[0] );

  //  All new only : Create new triangle
  append_new_triangle( hangingVertex[0], hangingVertex[1], hangingVertex[2] );

  remove_it(oldface);

}

//#############################################################################

void ConsistencyRefine2D :: checkFaceConsistency( Face *oldface ) 
{
  bitvec.reset();

  for( int i = 0; i < 3; i++) {
    Vertex *v1  = oldface->getNodeAt( (i+1)%3 );
    Vertex *v2  = oldface->getNodeAt( (i+2)%3 );
    hangingVertex[i] = edgemap->getVertexOnEdge( v1, v2 );
    if( hangingVertex[i] ) bitvec.set(i);
  }

}

//#############################################################################

int ConsistencyRefine2D :: atomicOp(Face *oldface)
{
    checkFaceConsistency( oldface );

    switch( bitvec.count() )
    {
	case 1:
	  numfacesRefined++;
	  makeConsistent1( oldface );
	  break;
	case 2:
	  numfacesRefined++;
	  makeConsistent2( oldface );
	  break;
	case 3:
	  numfacesRefined++;
	  makeConsistent3( oldface );
	  break;
     }


    remove_it( oldface );

    return 0;
}

//#############################################################################
void ConsistencyRefine2D :: makeConsistent()
{
  //**********************************************************************
  // The previous step, will leave some hanging nodes. So adjacent cells 
  // will be forced to refined. All the hanging nodes are attached 
  // directly to the opposite node of the triangle, thus creating two 
  // triangle. This may produce some bad triangle, which could be 
  // improved using edge swapping algorithm. In this process, only new 
  // edges are created and  no new vertices are introduced.
  //**********************************************************************

  size_t numfaces = mesh->getSize(2);

  for( size_t i = 0; i < numfaces ; i++) 
        atomicOp( mesh->getFaceAt(i) );
}

//#############################################################################

int CentroidRefine2D::refine_tri(Face *oldface)
{
  Vertex *vcenter = Vertex::newObject();
  vcenter->setXYZCoords( oldface->getCentroid() );

  Vertex *v1 = oldface->getNodeAt( 0 );
  Vertex *v2 = oldface->getNodeAt( 1 );
  Vertex *v3 = oldface->getNodeAt( 2 );

  append_new_node( vcenter );
  append_new_triangle( vcenter, v1, v2);
  append_new_triangle( vcenter, v2, v3);
  append_new_triangle( vcenter, v3, v1);

  remove_it( oldface );
  return 0;
}

///////////////////////////////////////////////////////////////////////////////

int  CentroidRefine2D::refine_quad(Face *oldface)
{
  Vertex *vcenter = Vertex::newObject();
  vcenter->setXYZCoords( oldface->getCentroid() );

  Vertex *v0 = oldface->getNodeAt( 0 );
  Vertex *v1 = oldface->getNodeAt( 1 );
  Vertex *v2 = oldface->getNodeAt( 2 );
  Vertex *v3 = oldface->getNodeAt( 3 );

  append_new_triangle( vcenter, v0, v1);
  append_new_triangle( vcenter, v1, v3);
  append_new_triangle( vcenter, v3, v2);
  append_new_triangle( vcenter, v2, v0);

  remove_it( oldface );

  return 0;
}

///////////////////////////////////////////////////////////////////////////

int CentroidRefine2D::atomicOp(Face *oldface)
{
   if( oldface->isRemoved() ) return 1;

   if( oldface->getSize(0) == 3 ) return refine_tri( oldface );
   if( oldface->getSize(0) == 4 ) return refine_quad( oldface );

   cout << "Warning: Element not supported for refinement " << endl;

   return 1;
}

///////////////////////////////////////////////////////////////////////////

int CentroidRefine2D::execute()
{

  for( int it = 0; it < numIterations; it++) {
       size_t numfaces = mesh->getSize(2);
       for( size_t i = 0; i < numfaces; i++)
            atomicOp( mesh->getFaceAt(i)  );
  }

  mesh->prune();
  mesh->enumerate(0);
  mesh->enumerate(2);

  assert( mesh->isSimple() );

  return 0;
}

///////////////////////////////////////////////////////////////////////////

int ObtuseRefine2D :: atomicOp( const Face *oldface)
{
/*
    int err;
    SimpleArray<iBase_EntityHandle> facenodes;
    iMesh_getEntAdj(mesh, facehandle, iBase_VERTEX, ARRAY_INOUT(facenodes), &err);

    Point3D pv0, pv1, pv2;
    Point3D vec1, vec2;

    iBase_EntityHandle vmid;
    double angle;
    for( int i = 0; i < 3; i++) {
         getVertexCoords( facenodes[(i+0)%3], pv0);
         getVertexCoords( facenodes[(i+1)%3], pv1);
         getVertexCoords( facenodes[(i+2)%3], pv2);
	 vec1 = Math::create_vector( pv2, pv0);
	 vec2 = Math::create_vector( pv1, pv0);
         angle = Math::getVectorAngle(vec1, vec2);
         if( angle > cutoffAngle) {
	    setVertexOnEdge(facenodes[(i+1)%2], facenodes[(i+2)%2], vmid );
	    return 0;
         }
    }
*/

    return 1;
}

///////////////////////////////////////////////////////////////////////////
int ObtuseRefine2D :: execute()
{
  edgemap = new RefinedEdgeMap;

  initialize();

  size_t numfaces = mesh->getSize(2);

  size_t ncount = 0;
  for( size_t i = 0; i < numfaces; i++) {
       int err = atomicOp( mesh->getFaceAt(i) );
       if( !err ) ncount++;
  }

  if( ncount ) {
      ConsistencyRefine2D refine(mesh, edgemap);
      refine.execute();
      finalize();
  } else 
      cout << "Warning: No triangle was refined " << endl;

  edgemap->clear();
  delete edgemap;

  return 0;
}

//////////////////////////////////////////////////////////////////////////////

int Refine2D14::refine_quad(Face *oldface)
{
  Vertex *v0 = oldface->getNodeAt( 0 );
  Vertex *v1 = oldface->getNodeAt( 1 );
  Vertex *v2 = oldface->getNodeAt( 2 );
  Vertex *v3 = oldface->getNodeAt( 3 );

  edgemap->setVertexOnEdge( v0,v1 );
  edgemap->setVertexOnEdge( v2,v3 );
  edgemap->setVertexOnEdge( v0,v2 );
  edgemap->setVertexOnEdge( v1,v3 );

  Vertex *v01  = edgemap->getVertexOnEdge( v0,v1 );
  Vertex *v23  = edgemap->getVertexOnEdge( v2,v3 );
  Vertex *v02  = edgemap->getVertexOnEdge( v0,v2 );
  Vertex *v13  = edgemap->getVertexOnEdge( v1,v3 );

  Vertex *vcenter = Vertex::newObject();
  vcenter->setXYZCoords( oldface->getCentroid() );
  
  append_new_quad( v0, v01, v02, vcenter);
  append_new_quad( v01, v1, vcenter, v13);
  append_new_quad( v02, vcenter, v2, v23);
  append_new_quad(vcenter, v13, v23, v3 );

  remove_it( oldface );

  return 0;
}

///////////////////////////////////////////////////////////////////////////////

int Refine2D14:: refine_tri( Face *oldface)
{
  Vertex *v0 = oldface->getNodeAt( 0 );
  Vertex *v1 = oldface->getNodeAt( 1 );
  Vertex *v2 = oldface->getNodeAt( 2 );

  edgemap->setVertexOnEdge(v0,v1);    
  edgemap->setVertexOnEdge(v1,v2);    
  edgemap->setVertexOnEdge(v2,v0);  

  Vertex *v01 = edgemap->getVertexOnEdge( v0, v1 );
  Vertex *v12 = edgemap->getVertexOnEdge( v1, v2 );
  Vertex *v20 = edgemap->getVertexOnEdge( v2, v0 );

  append_new_triangle( v0, v01, v20);
  append_new_triangle( v01, v1, v12);
  append_new_triangle( v12, v2, v20);
  append_new_triangle( v01, v12, v20);

  remove_it( oldface );

  return 0;

}

///////////////////////////////////////////////////////////////////////////////

int Refine2D14::atomicOp(Face *oldface)
{
   if( oldface->getSize(0) == 3 ) return refine_tri( oldface );
   if( oldface->getSize(0) == 4 ) return refine_quad( oldface );

   cout << "Warning: Element not supported for refinement " << endl;
   return 1;
}

///////////////////////////////////////////////////////////////////////////////

int Refine2D14 ::execute() 
{
  edgemap = new RefinedEdgeMap;
  initialize();

  size_t numfaces = mesh->getSize(2);

  size_t ncount = 0;
  for( size_t i = 0; i < numfaces; i++)  {
       int err = atomicOp( mesh->getFaceAt(i)  );
       if( !err ) ncount++;
  }

  if( ncount ) {
      ConsistencyRefine2D refine(mesh, edgemap);
      refine.execute();
      MeshRefine2D::finalize();
  }

  edgemap->clear();

  delete edgemap;

  return 0;
}

///////////////////////////////////////////////////////////////////////////////

int GradeRefine2D :: initialize()
{
   return 0;
}

///////////////////////////////////////////////////////////////////////////////

int GradeRefine2D :: atomicOp( const Vertex *apexVertex)
{
/*
    SimpleArray<iBase_EntityHandle> vneighs;
    iMesh_getEntAdj(mesh, apexVertex, iBase_VERTEX, ARRAY_INOUT(vneighs), &err);

    size_t numNeighs = vneighs.size();

    if( numNeighs == 0 ) return 0;

    vector<double> elen;
    elen.resize( numNeighs);

    for( int i = 0; i < numNeighs; i++) 
         elen[i] = length( vertex, vneighs[i] );

    sort( elen.begin(), elen.end() );

    double median_value = elen[numNeighs/2];

    for( int i = 0; i < numNeighs; i++) {
        if( elen[i] > 0.90*median_value) 
            setVertexOnEdge( vertex, vneighs[i]);
    }
*/
    return 1;
}
////////////////////////////////////////////////////////////////////////////////

int GradeRefine2D :: finalize()
{  
  return 0; 
}

////////////////////////////////////////////////////////////////////////////////

int GradeRefine2D :: execute()
{
    edgemap = new RefinedEdgeMap;

    initialize();

    size_t numnodes = mesh->getSize(0);

    size_t ncount = 0;
    for( size_t i = 0; i < numnodes; i++) {
         int err = atomicOp( mesh->getNodeAt(i) );
	 if( !err ) ncount++;
    }

    if( ncount ) {
        ConsistencyRefine2D refine(mesh, edgemap);
        refine.execute();
        finalize();
    }

    edgemap->clear();
    delete edgemap;

    return ncount;
}

///////////////////////////////////////////////////////////////////////////////

#ifdef TEST_MESHREFINE
int main(int argc, char **argv)
{
    iMesh_Instance mesh = read_off_file( "model.off");

//  LongestEdgeRefine2D  meshrefine;
    Refine2D14  meshrefine;
    meshrefine.setBoundarySplitFlag(0);
    meshrefine.setMesh( mesh );
    meshrefine.execute();

    write_off_file( mesh, "refine.off");
    
    return 0;
}

#endif

