#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <time.h>
#include <fcntl.h>
#include <sstream>
#include <queue>
#include <algorithm>
#include <string.h>

#include "ProjectShell.hpp"
#include "vec_utils.hpp"
 
#define ERROR(a) {if (iBase_SUCCESS != err) std::cerr << a << std::endl;}
#define ERRORR(a,b) {if (iBase_SUCCESS != err) {std::cerr << a << std::endl; return b;}}

const bool debug = false;
static bool equal_to(double d1, double d2) { return abs(d1 - d2) < 10e-7; }

int dbg = 0;
std::ofstream mout;// some debug file
double epsilon = 1.e-5; // cm, for coincident points in P, the intersection area


int getEdge(PSNode & v1, PSNode & v2, int & edgeId, PSEdge * edgesArr, int sizeArr)
{
  // not found return 0 
  std::list<int> & L = v1.edges;
  std::list<int>::iterator i;
  for (i=L.begin(); i!=L.end(); i++)
    {
      int edge = *i;
      int ae = abs(edge);
      if (ae>sizeArr || ae ==0 )
	{
	  std::cerr<<" bad index in getEdge\n";
	  exit(0);
	}

      PSEdge & pse = edgesArr[ae-1];
      //      int sig = (edge<0) ? -1 : 1;
      if (pse.v[0] == v1.id && pse.v[1] == v2.id)
	{
	  edgeId = edge;
	  return 1;
	}
      if (pse.v[0] == v2.id && pse.v[1] == v1.id)
	{
	  edgeId = edge;
	  return 1;
	}
    }
  return 0;// did not find any edge between v1 and v2 
  
}
ProjectShell::ProjectShell(iMesh_Instance mesh,  iBase_EntitySetHandle root_set, double direction[3])
  : m_mesh(mesh), m_hRootSet(root_set), 
    m_numNodes(0),
    m_numTriangles(0),
    m_xyz(NULL),  // original coordinates
    m_triangles(NULL), // original triangles
    ps_edges(NULL),
    ps_nodes(NULL),
    ps_triangles(NULL),
    m_xy(NULL), // 2d coordinates after projection
    m_redMesh(NULL),
    m_blueMesh(NULL),
    m_num2dPoints(0)
{
  m_direction[0] = direction[0];
  m_direction[1] = direction[1];
  m_direction[2] = direction[2];
}

ProjectShell::~ProjectShell()
{
  delete [] ps_edges;
  delete [] ps_nodes;
  delete [] ps_triangles;
  delete [] m_xyz;
  delete [] m_triangles;
  delete [] m_xy;
  delete [] m_redMesh;
  delete [] m_blueMesh;
}

int ProjectShell::project()
{
  double dist1= dist(m_direction);
  if (dist1==0)
    {
      std::cerr << " null direction \n";
      return 1;
    }
  // normalize direction
  for (int k=0; k<3; k++)
    {
      m_direction[k]/=dist1;
    }
  int ret = getMeshData();
  if (ret)
    {
      std::cout<< " bad mesh\n" ;
      return 1;
    }
  // we have now the 3D mesh
  // verify orientation of the triangles; is it manifold?
  ret =  checkMeshValidity();
  if (ret)
    {
      std::cout<< " bad orientation\n" ;
      return 1;
    }
  
  ret =  projectIn2D();
  if (ret)
    {
      std::cout<< " cannot project in 2d\n" ;
      return 1;
    }
  ret =  computeIntersections();
  if (ret)
    {
      std::cout<< " error in computing intersections\n" ;
      return 1;
    }
  return 0;
}

int ProjectShell::getMeshData()
{
 
  // get original coordinates of mesh
  iBase_EntityHandle *verts = NULL;
  int verts_alloc = 0;
  
  int vert_coords_alloc = 0;
  int vertex_coord_size;

  /* check storage order */
  int result;
  int this_order;
  iMesh_getDfltStorage(m_mesh, &this_order, &result);
  if (iBase_SUCCESS != result) {
    printf("failed to get preferred storage order in getMesh.\n");
    return 1;
  }
  

  /* now get the vertex coordinates from a vertex array */
  /* need to get the vertices in the model */
  verts = NULL;
  verts_alloc = 0;
  iMesh_getEntities(m_mesh, m_hRootSet, iBase_VERTEX, 
                    iMesh_POINT, &verts, &verts_alloc, &m_numNodes, &result);
  if (iBase_SUCCESS != result) {
    printf("failed to get vertices.\n");
    return 1;
  }
  /*
    iMesh_getNumOfType  	(  	iMesh_Instance   	 instance,
    const iBase_EntitySetHandle  	entity_set_handle,
    const int  	entity_type,
    int *  	num_type,
    int *  	err	 
    ) */
  int numEdges=0;
  iMesh_getNumOfType 	(m_mesh, m_hRootSet, iBase_EDGE, 
			 &numEdges, &result);
  		
  // get the triangles and the vertices in one shot
  iBase_EntityHandle *triangles = NULL;
  int triangles_alloc = 0;
  iBase_EntityHandle *vert_adj = NULL;
  int  vert_adj_alloc = 0, vert_adj_size;
  int * offsets = NULL, offsets_alloc = 0, indices_size;
  int * indices = NULL, indices_alloc = 0, offsets_size;
  iMesh_getAdjEntIndices( m_mesh, m_hRootSet, 
			  iBase_FACE, iMesh_TRIANGLE, iBase_VERTEX,
			  &triangles, &triangles_alloc, &m_numTriangles,
			  &vert_adj, &vert_adj_alloc, &vert_adj_size,
			  &indices, &indices_alloc, &indices_size,
			  &offsets, &offsets_alloc, &offsets_size,
			  &result );
  if (iBase_SUCCESS != result) {
    printf("failed to get triangles and vertices.\n");
    return 1;
  }
  

  /* get the coordinates in one array */
   
  vert_coords_alloc = 0;

  iMesh_getVtxArrCoords(m_mesh, vert_adj, vert_adj_size, iBase_INTERLEAVED, 
                        &m_xyz, &vert_coords_alloc, &vertex_coord_size, &result);
  if (iBase_SUCCESS != result || 3*m_numNodes!=vertex_coord_size) {
    printf("failed to get vertex coordinates of entities in getMeshData.\n");
    return 1;
  }
 

  // build list of edges; we cannot rely on iMesh to give them to us
  // can we force imesh to give the edges? It does not seem like that

  // the vertices are identified as the index in vert_adj 
  // indices are the triangles
  // i is index in triangles
  //  offsets[i], offsets[i+1] are offsets in indices
  // vertices of triangle [i] have indices  indices[offsets[i]+j], j=0:(offsets[i+1]-offsets[i])-1

  // first, determine the edges of triangle i
  if (offsets_size - 1 != m_numTriangles)
    {
      std::cerr<<"bad indexing\n";
      return 1;
    }
  int i=0;// index used everywhere
  for (i=0; i<m_numTriangles; i++)
    {
      if (offsets[i+1]-offsets[i] !=3)
	{
	  std::cerr<< "not a triangle\n";
	  return 1;
	}
    }
  // build edges from triangle information
  ps_edges = new PSEdge [m_numTriangles*3];// overestimate; probably only half will be created
  ps_nodes = new PSNode [m_numNodes];
  ps_triangles = new PS3DTriangle [m_numTriangles];
  //
  for (i=0; i<m_numNodes; i++)
    {
      PSNode & psn = ps_nodes[i];
      psn.id = i;// this can be found by address in the array, but why bother
      for (int j=0; j<3; j++)
	{
	  psn.xyz[j] = m_xyz[3*i+j];
	}
    }
  // index in edge: 
  int edgeIndex = 0;
  int j=0;
  for (j=0; j<m_numTriangles; j++)
    {
      PS3DTriangle & curTria = ps_triangles[j];
      int ii=0;
      for (ii=0; ii<3; ii++)
        curTria.v[ii]=indices[offsets[j]+ii];
      for (ii=0; ii<3; ii++)
	{
	  PSNode & v1= ps_nodes[curTria.v[ii]];
          int nextii = ii+1;
	  if (ii==2)
	    nextii=0;
          PSNode & v2= ps_nodes[curTria.v[nextii]];
	  // is there an edge from v1 to v2, or from v2 to v1
          int edgeId;
	  int exi = getEdge(v1, v2, edgeId, ps_edges, edgeIndex);// will be created if not existing already
          if (exi)
	    {
	      curTria.e[ii] = edgeId;
              PSEdge & foundEdge = ps_edges[abs(edgeId)-1];
	      foundEdge.used++;
	      if(foundEdge.used>2)
		{
		  std::cerr<< " bad indexing in ps_edge; exit\n";
		  exit(1);
		}
	      foundEdge.t[1]=j; // mark the triangle it belongs to
	    }
	  else
	    {
	      int cId = curTria.e[ii]=edgeIndex+1;
	      PSEdge & newEdge = ps_edges[edgeIndex];
	      newEdge.v[0] = curTria.v[ii];
	      newEdge.v[1]=curTria.v[nextii];
	      v1.edges.push_back(cId); // positive means starts at vertex
	      v2.edges.push_back(-cId);// negative means incident to the vertex
	      edgeIndex++;
	      newEdge.used=1;
	      newEdge.t[0] = j; // index for triangles
	    }
	}
    }
  m_numEdges = edgeIndex;
  // fill up now the triangle neighbor information
  for (j=0; j<m_numTriangles; j++)
    {
      PS3DTriangle & tri = ps_triangles[j];
      for (int k=0; k<3; k++)
	{
	  int edgeId = tri.e[k];
	  int ae = abs(edgeId);
	  PSEdge & edge = ps_edges[ae-1];
	  int t0 = edge.t[0], t1 = edge.t[1];
	  if (t0==j)
	    tri.t[k] = t1;
	  else if (t1==j)
	    tri.t[k] = t0;
	}
    }

  free(verts);
  free(triangles);
  free(vert_adj);
  free (indices);
  free (offsets);

  //free(vert_coords);
  

  return 0;
}
int ProjectShell::checkMeshValidity()
{
  // basically, see if the edges are all used 2 times, and each is is used in each direction
  int i=0;
  for (i=0; i<m_numEdges; i++)
    {
      PSEdge & edge = ps_edges[i];
      if (edge.used!=2)
	{
	  std::cerr << " edge not used exactly 2 times\n";
	  return 1;
	}
      int edgeId = i+1;
      int sum=0;
      for (int k=0; k<2; k++)
	{
	  PS3DTriangle & tr1=ps_triangles[edge.t[k]];
	  for (int kk=0; kk<3; kk++)
	    {
	      int ec= tr1.e[kk];
	      if (abs(ec)==edgeId)
		{
		  sum+=ec;
		  break;
		}
	
	    }
	}
      if (sum!=0)
	{
	  std::cerr<< " edge " << edgeId << " not used exactly 2 times\n";
	  return 1;
	}
     
    }
  // check that all triangles have 3 neighbors
  for (i=0; i<m_numTriangles; i++)
    {
      PS3DTriangle & tr = ps_triangles[i];
      for (int k=0; k<3; k++)
	{
	  int t = tr.t[k];
	  if (t==i || t < 0 || t >= m_numTriangles)
	    {
	      std::cerr<< " triangle " << i << " does not have 3 neighbors \n";
	      return 1;
	    }
	}
    }
 
  return 0; // everything is OK
}
int ProjectShell::projectIn2D()
{
  // considering the direction, classify each projected triangle as red or blue
  // the red ones are positive (inbound), blue ones are negative
  // first decide the other 2 vectors
  double vv[3]={0.,-1., 0.};
  cross(m_dirX, m_direction, vv);
  
  double d1 = dist(m_dirX);
  if (d1<1.e-5)// consider another direction
    {
      vv[1]=0.; vv[2]=-1.;
      cross(m_dirX, m_direction, vv);
      d1 = dist(m_dirX);
      if (d1<1.e-5)
	{
	  std::cerr << "cannot find a suitable direction; abort\n";
	  return 1;
	}
    }
  int k=0;
  for (k=0; k<3; k++)
    m_dirX[k]/=d1;
  cross(m_dirY, m_direction, m_dirX);
  // dirY must be already normalized, but why worry
  d1 = dist(m_dirY);
  if (d1==0.)
    {
      std::cerr<< " get out of here, it is hopeless\n";
      return 1;
    }
  for (k=0; k<3; k++)
    m_dirY[k]/=d1;
  
  // now do the projection in 2d
  //   we have 3 vectors, m_dirX, m_dirY, m_direction
  // for every point, A, the projection in xy plane will be just
  // xA = A . u; yA = A . v  (dirX and dirY)

  // also, it is important to compute the normal 
  // some triangles will be reverted and some may be projecting to a line

  // coordinates of the projected nodes on 2D
  // what could be a good estimate for total number of points in 2D? 
  //  we start with 2d capacity 3* m_numNodes 
  m_2dcapacity = m_numNodes*3;
  m_num2dPoints = m_numNodes; // directly project 3d nodes in 2d, keep the ids
  m_xy = new double [3*m_2dcapacity];
  // this array will be parallel with the original mesh
  // new nodes will appear after intersection computation
  // triangles are characterized by their orientation: negative, positive or 0

  // the neighbors will be preserved
  // some edges will be collapsed, and also some triangles
  // in those cases, what will be the neighbors?
  // flag them with a high number ()

  //m_finalNodes.resize(3*m_numNodes); // the actual size is n_numNodes first
  for ( k=0; k<m_numNodes; k++)
    {
      // double xx= dot( ps_nodes[k].xyz, m_dirX);
      // double yy= dot( ps_nodes[k].xyz, m_dirY);
      // _finalNodes[k].x =  dot( ps_nodes[k].xyz, m_dirX);
      // _finalNodes[k].y =  dot( ps_nodes[k].xyz, m_dirY);
      m_xy[3*k] = dot( ps_nodes[k].xyz, m_dirX);
      m_xy[3*k+1] = dot( ps_nodes[k].xyz, m_dirY);
      //  m_finalNodes[k].x =  m_xy[2*k];
      //  m_finalNodes[k].y =  m_xy[2*k+1];
    }
  for (k=0; k<m_2dcapacity; k++)
    m_xy[3*k+2] = 0;// the z ccordinate is always 0 !!
  // if any edges are collapsed, the corresponding triangles should be collapsed too
  // we will collapse the triangles first, then the edges
  // when a triangle is collapsed on an edge, it should break in 2 other triangles
  // we should really form another triangle array, in 2D 

  // first, loop over edges and see which are eliminated;
  
  double *edgeLen2D=new double [m_numEdges];
  for ( k=0; k<m_numEdges; k++)
    {
      
      int v0 = ps_edges[k].v[0];
      int v1 = ps_edges[k].v[1];
      edgeLen2D[k] = dist2(&m_xy[3*v0], &m_xy[3*v1]); 
    }
  double *triArea = new double [m_numTriangles];
  for ( k=0; k<m_numTriangles; k++)
    {
      const int * pV =&( ps_triangles[k].v[0]);
      triArea[k] = area2D( &m_xy[ 3*pV[0]],  &m_xy[ 3*pV[1]],  &m_xy[ 3*pV[2]]);
    }
  // if an edge is length 0, we must collapse the nodes, and 2 triangles
  // first construct a new 2d mesh
  // we will have to classify all triangles, edges
  // for each triangle we need to know its original triangle

  // first mark all triangles; then count positive, negative and zero area (some tolerance is needed)
  int numNeg = 0, numPos = 0, numZero = 0;
  // those that are negative will be reverted in the new array
  int * newTriId = new int [m_numTriangles];
  for (k=0; k<m_numTriangles ; k++)
    {
      if (triArea[k] > 0)
	{
	  //numPos++;
          newTriId[k] = numPos++;
	}
      else if (triArea[k] <0)
	{
	  //numNeg ++;
	  newTriId[k] = numNeg++;
	}
      else
	{
	  numZero ++;
	  newTriId[k] = -1;// receive no Id, as it will not be carried along
	}
    }
  // there are 2 groups: red (positive), blue (negative)
  // revert the negative ones, and decide the neighbors

  // red triangles are positive, blue are negative
  // the negative ones need to be reverted; revert the nodes, first, then edges
  // do we really need the edges? or just the neighboring triangles?
  // if a neighbor is the other sign or zero, will become boundary marker (large number, m_numTriangles+1) 
  // 
  // we need to find first 2 triangles (red and blue) that are intersecting
  // they will be the seeds
  m_redMesh = new PSTriangle2D [numPos] ;
  m_blueMesh = new  PSTriangle2D [numNeg];
  m_numPos = numPos;
  m_numNeg = numNeg;
  // do another loop, to really build the 2D mesh we start with
  // we may have potential starting triangles for the marching along, if the sign of one of the neighbors is 
  // different
  for (k=0; k<m_numTriangles ; k++)
    {
      PS3DTriangle & orgTria = ps_triangles[k];
      if (triArea[k] > 0)
	{
	  //numPos++;
          PSTriangle2D & redTria= m_redMesh[newTriId[k]];
	  redTria.oldId = k ; // index 
	  redTria.area = triArea[k];
          for (int j=0; j<3; j++)
	    {
	       
              // copy the edge information too
	      redTria.e[j] = orgTria.e[j]; 
              // qualify the neighbors
	      //
	      int t = orgTria.t[j];
	      redTria.v[j] = orgTria.v[j];
	      if (triArea[t]>0)
		{
		  redTria.t[j] = newTriId[t];// will be the index in red mesh 
		}
	      else
		{
		  redTria.t[j] = numPos; // marker for boundary
		}
	    }
	  
	}
      else if (triArea[k] <0)
	{
	  //numNeg ++;
	  PSTriangle2D & blueTria= m_blueMesh[newTriId[k]];
	  blueTria.oldId = k;
	  blueTria.area =triArea[k]; // this is for debugging I think
          for (int j=0; j<3; j++)
	    {
	      // copy the edge information too
	      blueTria.e[j] = orgTria.e[j];
	      // qualify the neighbors
	      // 
	      int t = orgTria.t[j];
	      blueTria.v[j] = orgTria.v[j];
	      if (triArea[t]<0)
		{
		  blueTria.t[j] = newTriId[t];// will be the index in red mesh 
		}
	      else
		{
		  blueTria.t[j] = numNeg; // marker for boundary
		}
	    }
     
	}
      else
	{
	  // numZero ++;
	  // newTriId[k] = -1;// receive no Id, as it will not be carried along
	  // nothing to do for null triangles
	}
    }
  // revert the blue triangles, so they become positive oriented too
  for (k=0; k<numNeg; k++)
    {
      // 
      PSTriangle2D & blueTri = m_blueMesh[k];
      int tmp = blueTri.v[1];
      blueTri.v[1] = blueTri.v[2];
      blueTri.v[2] = tmp;
      // this is really stupid: 
      // we should have switched triangle 1 and 3, not 2 and 3
      // hard to catch
      tmp = blueTri.t[0];
      blueTri.t[0] = blueTri.t[2];
      blueTri.t[2] = tmp;
      // reverse the edges too
      tmp = blueTri.e[0];
      blueTri.e[0] = -blueTri.e[2];
      blueTri.e[1] = -blueTri.e[1];
      blueTri.e[2] = -tmp; 
    }
  // at this point, we have red triangles and blue triangles, and 2d nodes array
  // we will start creating new triangles, step by step, pointing to the nodes in _finalNodes
  m_numCurrentNodes = m_numNodes;
  if (dbg)
    {
      std::ofstream fout("dbg.m"); 
      fout << "P=[ \n";
      for (int i=0; i<m_numNodes; i++)
	{
	  fout << m_xy[3*i] << " " << m_xy[3*i+1] << "\n";
	}
      fout << "];\n ";
      fout << "Ta=[ \n";
      for (int k=0; k<m_numPos; k++)
	{
	  PSTriangle2D & redTri = m_redMesh[k];
	  for (int j=0; j<3; j++)
	    fout << redTri.v[j]+1 << " " ;
	  for (int jj=0; jj<3; jj++)
	    {
	      fout << redTri.t[jj]+1 << " " ;
	    }
	  fout << "\n";
	}
      fout << "]; \n";
      fout << "Tb=[ \n";
      for (int kk=0; kk<m_numNeg; kk++)
	{
	  PSTriangle2D & blueTri = m_blueMesh[kk];
	  for (int j=0; j<3; j++)
	    fout << blueTri.v[j]+1 << " " ;
	  for (int jj=0; jj<3; jj++)
	    {
	      fout << blueTri.t[jj]+1 << " " ;
	    }
	  fout << "\n";
	}
      fout << "]; \n";
      fout.close();
      
    }
  delete [] edgeLen2D;
  delete [] triArea;
  delete [] newTriId;
  return 0; 
}

int ProjectShell::computeIntersections()
{
  // will start at 2 triangles, and advance an intersection front (2 queues, one for red triangles, one for blue ..)
  // will mark a red triangle that is processed, and add to the queue
  // for each red triangle will find all the blue triangles that are intersected
  // find first 2 triangles that intersect: these will be the seeds for intersection

  int startRed=0;
  int startBlue = 0;
  for (startBlue = 0; startBlue<m_numNeg; startBlue++)
    {
      double area = 0;
      // if area is > 0 , we have intersections
      double P[24]; // max 6 points, but it may grow bigger; why worry
      int nP = 0;
      int n[3];// sides
      computeIntersectionBetweenRedAndBlue(/* red */0, startBlue, P, nP, area,n);
      if (area>0)
	break; // found 2 triangles that intersect; these will be the seeds
    }
  if (startBlue==m_numNeg)
    {
      // can't find any triangle stop
      exit (1);
    }
  // on the red edges, we will keep a list of new points (in 2D)
  // they will be used to create or not new points for points P from intersection
  // (sometimes the points P are either on sides, or on vertices of blue or even red triangles)

  /*
    matlab code: 
    function M=InterfaceMatrix(Na,Ta,Nb,Tb);
    % INTERFACEMATRIX projection matrix for nonmatching triangular grids 
    %   M=InterfaceMatrix(Na,Ta,Nb,Tb); takes two triangular meshes Ta
    %   and Tb with associated nodal coordinates in Na and Nb and
    %   computes the interface projection matrix M

    bl=[1];                        % bl: list of triangles of Tb to treat
    bil=[1];                       % bil: list of triangles Ta to start with
    bd=zeros(size(Tb,1)+1,1);      % bd: flag for triangles in Tb treated 
    bd(end)=1;                     % guard, to treat boundaries
    bd(1)=1;                       % mark first triangle in b list.
    M=sparse(size(Nb,2),size(Na,2));
    while length(bl)>0
    bc=bl(1); bl=bl(2:end);      % bc: current triangle of Tb 
    al=bil(1); bil=bil(2:end);   % triangle of Ta to start with
    ad=zeros(size(Ta,1)+1,1);    % same as for bd
    ad(end)=1;
    ad(al)=1; 
    n=[0 0 0];                   % triangles intersecting with neighbors
    while length(al)>0
    ac=al(1); al=al(2:end);    % take next candidate
    [P,nc,Mc]=Intersect(Nb(:,Tb(bc,1:3)),Na(:,Ta(ac,1:3)));
    if ~isempty(P)             % intersection found
    M(Tb(bc,1:3),Ta(ac,1:3))=M(Tb(bc,1:3),Ta(ac,1:3))+Mc;
    t=Ta(ac,3+find(ad(Ta(ac,4:6))==0)); 
    al=[al t];               % add neighbors 
    ad(t)=1;
    n(find(nc>0))=ac;        % ac is starting candidate for neighbor  
    end
    end
    tmp=find(bd(Tb(bc,4:6))==0); % find non-treated neighbors
    idx=find(n(tmp)>0);          % take those which intersect
    t=Tb(bc,3+tmp(idx));
    bl=[bl t];                   % and add them
    bil=[bil n(tmp(idx))];       % with starting candidates Ta
    bd(t)=1;
    end
  */
  std::queue<int> blueQueue; // these are corresponding to Ta,
  blueQueue.push( startBlue);
  std::queue<int> redQueue;
  redQueue.push (startRed);
  // the flags are used for marking the triangles already considered
  int * blueFlag = new int [m_numNeg+1]; // number of blue triangles + 1, to account for the boundary
  int k=0;
  for (k=0; k<m_numNeg; k++)
    blueFlag[k] = 0;
  blueFlag[m_numNeg] = 1; // mark the "boundary"; stop at the boundary
  blueFlag[startBlue] = 1; // mark also the first one
  // also, red flag is declared outside the loop
  int * redFlag = new int [m_numPos+1];
  
  if (dbg)
    mout.open("patches.m");
  while( !blueQueue.empty() )
    {
      int n[3]; // flags for the side : indices in red mesh start from 0!!! (-1 means not found )
      for (k=0; k<3; k++)
	n[k] = -1; // a paired red not found yet for the neighbors of blue 
      int currentBlue = blueQueue.front();
      blueQueue.pop();
      for (k=0; k<m_numPos; k++)
	redFlag[k] = 0;
      redFlag[m_numPos] = 1; // to guard for the boundary
      int currentRed = redQueue.front(); // where do we check for redQueue???? 
      // red and blue queues are parallel
      redQueue.pop();// 
      redFlag[currentRed] = 1; // 
      std::queue<int> localRed;
      localRed.push(currentRed);
      while( !localRed.empty())
	{
	  // 
	  int redT = localRed.front();
	  localRed.pop();
	  double P[24], area;
          int nP = 0; // intersection points
	  int nc[3]= {0, 0, 0}; // means no intersection on the side (markers)
          computeIntersectionBetweenRedAndBlue(/* red */redT, currentBlue, P, nP, area,nc);
	  if (nP>0) 
	    {
	      // intersection found: output P and original triangles if nP > 2
	      if (dbg && area>0)
		{
		  //std::cout << "area: " << area<< " nP:"<<nP << " sources:" << redT+1 << ":" << m_redMesh[redT].oldId <<
		  //  " " << currentBlue+1<< ":" << m_blueMesh[currentBlue].oldId << std::endl;
		}
              if (dbg)
		{
		  mout << "pa=[\n";
		  
		  for (k=0; k<nP; k++)
		    {
		     
		      mout <<  P[2*k] << "\t " ;
		    }
	
		  mout << "\n";
		  for (k=0; k<nP; k++)
		    {
		    
		      mout << P[2*k+1] << "\t "; 
		    }
	
		  mout << " ]; \n";
		  mout << " patch(pa(1,:),pa(2,:),'m');       \n";
		}
	      // add neighbors to the localRed queue, if they are not marked
	      for (int nn= 0; nn<3; nn++)
		{
		  int neighbor = m_redMesh[redT].t[nn];
		  if (redFlag[neighbor] == 0)
		    {
		      localRed.push(neighbor);
		      redFlag[neighbor] =1; // flag it to not be added anymore
		    }
		  // n(find(nc>0))=ac;        % ac is starting candidate for neighbor
		  if (nc[nn]>0) // intersected  
		    n[nn] = redT;// start from 0!! 
		}
              if (nP>1) // this will also construct triangles, if needed
		findNodes(redT, currentBlue, P, nP);
	    }
	  
	}
      for (int j=0; j<3; j++)
	{
	  int blueNeigh = m_blueMesh[currentBlue].t[j];
	  if (blueFlag[blueNeigh]==0 && n[j] >=0 ) // not treated yet and marked as a neighbor
	    {
              // we identified triangle n[j] as intersecting with neighbor j of the blue triangle
	      blueQueue.push(blueNeigh);
	      redQueue.push(n[j]);
	      if (dbg)
		std::cout << "new triangles pushed: blue, red:" << blueNeigh+1 << " " << n[j]+1 << std::endl;
	      blueFlag[blueNeigh] = 1;
	    }
	}
    }
  delete [] redFlag;
  redFlag = NULL;

  delete [] blueFlag; // get rid of it
  blueFlag = NULL;
  if (dbg)
    mout.close();
  return 0;
}

// this is a local method
int EdgeIntersections(double * red, double * blue, int mark[3], double * points, int & nPoints)
{
  /* EDGEINTERSECTIONS computes edge intersections of two triangles
     [P,n]=EdgeIntersections(X,Y) computes for the two given triangles  * red 
     and blue ( stored column wise )
     (point coordinates are stored column-wise, in counter clock
     order) the points P where their edges intersect. In addition,
     in n the indices of which neighbors of red  are also intersecting
     with blue are given.
  */

  // points is an array with 12 slots   (12 * 2 doubles) 
  nPoints = 0; 
  mark[0] = mark[1] = mark[2]=0 ; // no neighbors of red involved yet
  /*for i=1:3                            % find all intersections of edges
    for j=1:3
    b=Y(:,j)-X(:,i);               
    A=[X(:,mod(i,3)+1)-X(:,i) -Y(:,mod(j,3)+1)+Y(:,j)];
    if rank(A)==2                   % edges not parallel
    r=A\b;
    if r(1)>=0 & r(1)<=1 & r(2)>=0 & r(2)<=1,  % intersection found
    k=k+1; P(:,k)=X(:,i)+r(1)*(X(:,mod(i,3)+1)-X(:,i)); n(i)=1;
    end;
    end; 
    end;
    end;*/
  for (int i=0; i<3; i++)
    {
      for (int j=0; j<3; j++)
	{
	  double b[2];
	  double a[2][2]; // 2*2
	  int iPlus1 = (i+1)%3;
          int jPlus1 = (j+1)%3;
	  for (int k=0; k<2; k++)
	    {
	      b[k] = blue[2*j+k] - red[2*i+k];
	      // row k of a: a(k, 0), a(k, 1)
	      a[ k ][0]    = red  [ 2*iPlus1 + k] - red [2*i + k];
	      a[ k ][1]    = blue [ 2 * j + k] - blue [2*jPlus1 +k];
               
	    }
	  double delta = a[0][0]*a[1][1] - a[0][1]*a[1][0];
	  if (delta != 0.)
	    {
	      // not parallel
	      double alfa = (b[0]*a[1][1]-a[0][1]*b[1])/delta;
	      double beta = (-b[0] * a[1][0] + b[1] * a[0][0])/delta;
	      if (0<= alfa && alfa <=1. && 0<= beta && beta <= 1.)
		{
		  // the intersection is good
		  for (int k=0; k<2; k++)
		    {
		      points[2*nPoints+k] = red[2*i+k] + alfa*(red[2*iPlus1+k]-red[2*i+k]);
		    }
		  mark[i] = 1; // so neighbor number i will be considered too.
		  nPoints++;
		}
	    }
	    
	}
    }
  return 0;
}

int borderPointsOfXinY(double * X, double * Y, double * P)
{
  // 2 triangles, 3 corners, is the corner of X in Y?
  // Y must have a positive area
  /*
    function P=PointsOfXInY(X,Y);
    % POINTSOFXINY finds corners of one triangle within another one
    %   P=PointsOfXInY(X,Y); computes for the two given triangles X
    %   and Y (point coordinates are stored column-wise, in counter clock
    %   order) the corners P of X which lie in the interior of Y.

    k=0;P=[];
    v0=Y(:,2)-Y(:,1); v1=Y(:,3)-Y(:,1);  % find interior points of X in Y
    d00=v0'*v0; d01=v0'*v1; d11=v1'*v1;  % using baricentric coordinates
    id=1/(d00*d11-d01*d01);
    for i=1:3 
    v2=X(:,i)-Y(:,1); d02=v0'*v2; d12=v1'*v2; 
    u=(d11*d02-d01*d12)*id; v=(d00*d12-d01*d02)*id;
    if u>=0 & v>=0 & u+v<=1            % also include nodes on the boundary
    k=k+1; P(:,k)=X(:,i);
    end;
    end;
  */
  int extraPoint = 0;
  for (int i=0; i<3; i++)
    {
      // compute twice the area of all 3 triangles formed by a side of Y and a corner of X; if one is negative, stop
      double A[2];
      for (int k=0; k<2; k++)
	A[k] = X[2*i+k];
      int inside=1;
      for (int j=0; j<3; j++)
	{
	  double B[2], C[2];
	  for (int k=0; k<2; k++)
	    {
	      B[k] = Y[2*j+k];
	      int j1 = (j+1)%3;
	      C[k] = Y[2*j1+k];
	    }

	  double area2 = (B[0]-A[0])*(C[1]-A[1])-(C[0]-A[0])*(B[1]-A[1]);
	  if (area2<0.)
	    {
	      inside = 0; break;
	    }
	}
      if (inside)
	{
	  P[extraPoint*2  ] = A[0];
	  P[extraPoint*2+1] = A[1];
	  extraPoint++;
	}
    }
  return extraPoint;
}

/* 
   function P=SortAndRemoveDoubles(P);
   % SORTANDREMOVEDOUBLES sort points and remove duplicates
   %   P=SortAndRemoveDoubles(P); orders polygon corners in P counter
   %   clock wise and removes duplicates

   ep=10*eps;                           % tolerance for identical nodes
   m=size(P,2); 
   if m>0                              
   c=sum(P')'/m;                      % order polygon corners counter 
   for i=1:m                          % clockwise
   d=P(:,i)-c; ao(i)=angle(d(1)+sqrt(-1)*d(2));
   end;
   [tmp,id]=sort(ao); 
   P=P(:,id);
   i=1;j=2;                           % remove duplicates
   while j<=m
   if norm(P(:,i)-P(:,j))>ep
   i=i+1;P(:,i)=P(:,j);j=j+1;
   else
   j=j+1;
   end;
   end;
   P=P(:,1:i);
   end;
*/

int swap (double * p , double * q)
{
  double tmp = *p;
  *p = *q;
  *q = tmp;
  return 0;
}
int SortAndRemoveDoubles(double * P, int & nP)
{
  if (nP<2)
    return 0; // nothing to do
  // center of gravity for the points
  double c[2] = {0., 0.};
  int k=0;
  for (k=0; k<nP; k++)
    {
      c[0]+=P[2*k];
      c[1]+=P[2*k+1];
    }
  c[0]/=nP;
  c[1]/=nP;
  double angle[12]; // could be at most 12 points; much less usually
  for (k=0; k<nP; k++)
    {
      double x = P[2*k] - c[0], y = P[2*k+1] - c[1] ;
      if ( x!= 0. || y!=0.)
        angle[k] = atan2(y, x);
      else
	{
	  angle[k] = 0;
	  // this means that the points are on a line, or all coincident // degenerate case
	}
    } 
  // sort according to angle; also eliminate close points
  int sorted = 1;
  do
    {
      sorted = 1;
      for(k=0; k<nP-1; k++)
	{
	  if (angle[k]>angle[k+1])
	    {
	      sorted = 0;
	      swap ( angle+k, angle+k+1);
	      swap ( P+(2*k), P+(2*k+2));
	      swap ( P+(2*k+1), P+(2*k+3));
	    }
	}
    }
  while (!sorted);
  // eliminate doubles

  int i=0, j=1; // the next one; j may advance faster than i
  // check the unit 
  // double epsilon = 1.e-5; // these are cm; 2 points are the same if the distance is less than 1.e-5 cm
  while (j<nP)
    {
      double d2 = dist2 ( &P[2*i], &P[2*j]);
      if (d2 > epsilon)
	{
	  i++;  
	  P[2*i] = P[2*j];
	  P[2*i+1] = P[2*j+1]; 
	}
      j++;
    }
  // test also the last point with the first one (index 0)
  
  double d2 = dist2(P, &P[2*i]); // check the first and last points (ordered from -pi to +pi)
  if (d2 > epsilon)
    {
      nP = i+1;
    }
  else
    nP = i; // effectively delete the last point (that would have been the same with first)
  if (nP==0)
    nP=1; // we should be left with at least one point we already tested if nP is 0 originally
  return 0;
}
// this method computed intersection between 2 triangles: will output n points, area, affected sides
int ProjectShell::computeIntersectionBetweenRedAndBlue(int red, int blue, double * P, 
						       int & nP, double & area,  int mark[3])
{
  // the points will be at most 9; they will describe a convex patch, after the points will be ordered and 
  // collapsed (eliminate doubles)
  // the area is not really required
  PSTriangle2D & redTri = m_redMesh[red];
  PSTriangle2D & blueTri = m_blueMesh[blue];
  double redTriangle[6];// column wise
  double blueTriangle[6];
  for (int i=0; i<3; i++)
    {
      // node i, 2 coordinates
      for (int k=0; k<2; k++)
	{
	  int node = redTri.v[i];
	  redTriangle[2*i+k] = m_xy[3*node+k];

	  node = blueTri.v[i];
	  blueTriangle[2*i+k] =  m_xy[3*node+k];
	}
    }
  /* Matlab source code:
     function [P,n,M]=Intersect(X,Y);
     % INTERSECT intersection of two triangles and mortar contribution
     %   [P,n,M]=Intersect(X,Y); computes for the two given triangles X and
     %   Y (point coordinates are stored column-wise, in counter clock
     %   order) the points P where they intersect, in n the indices of
     %   which neighbors of X are also intersecting with Y, and the local
     %   mortar matrix M of contributions of the P1 elements X on the P1
     %   element Y. The numerical challenges are handled by including
     %   points on the boundary and removing duplicates at the end.

     [P,n]=EdgeIntersections(X,Y);
     P1=PointsOfXInY(X,Y);
     if size(P1,2)>1                      % if two or more interior points
     n=[1 1 1];                         % the triangle is candidate for all 
     end                                  % neighbors
     P=[P P1];
     P=[P PointsOfXInY(Y,X)];
     P=SortAndRemoveDoubles(P);           % sort counter clock wise
     M=zeros(3,3);
     if size(P,2)>0
     for j=2:size(P,2)-1                % compute interface matrix
     M=M+MortarInt(P(:,[1 j j+1]),X,Y);
     end;
     patch(P(1,:),P(2,:),'m')           % draw intersection for illustration
     %  H=line([P(1,:) P(1,1)],[P(2,:),P(2,1)]);
     %  set(H,'LineWidth',3,'Color','m');
     pause(0)
     end;


  */
  //we do not really need the mortar matrix
 
  //int n[3]={0, 0, 0};// no intersection of side red with blue
  //double area= 0.;
  // X corresponds to blue, Y to red
  nP=0; // number of intersection points
  int ret = EdgeIntersections(blueTriangle, redTriangle, mark, P, nP);
  if (ret!=0)
    exit(1);// some unforeseen error
  
  int extraPoints = borderPointsOfXinY(blueTriangle, redTriangle, &(P[2*nP]));
  if (extraPoints>1)
    {
      mark[0] = mark[1] = mark[2]=1;
    }				      		     
  nP+=extraPoints;
  extraPoints =  borderPointsOfXinY(redTriangle, blueTriangle, &(P[2*nP]));
  nP+=extraPoints;

  // now sort and orient the points in P, such that they are forming a convex polygon
  // this will be the foundation of our new mesh
  //
  SortAndRemoveDoubles (P, nP); // nP should be at most 6 in the end 
  // if there are more than 3 points, some area will be positive
  area = 0.;
  if (nP>=3)
    {
      for (int k=1; k<nP-1; k++)
	area += area2D(P, &P[2*k], &P[2*k+2]);
    }
  return 0; // no error
}
int ProjectShell::findNodes(int red, int blue, double * iP, int nP)
{
  // first of all, check against red and blue vertices
  //
  if (dbg)
    {
      std::cout<< "red, blue, nP, P " << red << " " << blue << " " << nP <<"\n";
      for (int n=0; n<nP; n++)
	std::cout << " \t" << iP[2*n] << "\t" << iP[2*n+1] << "\n";

    }
  int * foundIds = new int [nP];
  for (int i=0; i<nP; i++)
    {
      double * pp = &iP[2*i];// iP+2*i
      int found = 0; 
      // first, are they on vertices from red or blue?
      PSTriangle2D & redTri = m_redMesh[red];
      int j=0;
      for (j=0; j<3&& !found; j++)
	{
          int node = redTri.v[j];
          double d2 = dist2( pp, m_xy+(3*node) );	
          if (dbg && i==0)
	    std::cout<< "  red node " << j << " " << node << " " << m_xy[3*node] << " " << m_xy[3*node+1] << " d2:" << d2 <<  " \n";
	  if (d2<epsilon)
	    {
              foundIds[i] = node; // no new node
	      found = 1;
	    }
	}
      PSTriangle2D & blueTri = m_blueMesh[blue];
      for (j=0; j<3 && !found; j++)
        {
          int node = blueTri.v[j];
          double d2 = dist2( pp, m_xy+(3*node) );	
          if (dbg && i==0) 
	    std::cout<< "  blu node " << j << " " << node << " " << m_xy[3*node] << " " << m_xy[3*node+1] << " d2:" << d2 <<  " \n";
	  if (d2<epsilon)
	    {
              foundIds[i] = node; // no new node
	      found = 1;
            }
        }
      if (!found)
        {
	  // find the edge it belongs, first
	  //
	  for (j=0; j<3; j++)
	    {
	      int edge = redTri.e[j];
	      int ae = abs(edge);
	      int v1 = ps_edges[ae-1].v[0];
	      int v2 = ps_edges[ae-1].v[1];
	      double area = area2D (&m_xy[3*v1], &m_xy[3*v2], pp);
	      if (dbg)
		std::cout << "   edge " << j << ": " << edge << " " << v1 << " " << v2 << "  area : " << area << "\n"; 
	      if ( fabs(area) < epsilon*epsilon  )
		{
		  // found the edge; now find if there is a point in the list here
		  std::list<int> & expts = ps_edges[ae-1].extraNodes;
		  // if the points pp is between extra points, then just give that id
		  // if not, create a new point, (check the id)
		  //
		  std::list<int>::iterator it;
		  for ( it = expts.begin(); it!=expts.end() && !found; it++)
		    {
		      int pnt = *it;
		      double d2 = dist2(pp, &m_xy[3*pnt]);
		      if (d2<epsilon)
			{
			  found = 1;
			  foundIds[i] = pnt;
			}
		    }
		  if (!found)
		    {
		      // create a new point in 2d (at the intersection)
		      foundIds [i] = m_num2dPoints;
		      expts.push_back(m_num2dPoints);
		      if (m_2dcapacity < m_num2dPoints+1)
			{
			  // exit, underestimate number of intersection points
			  std::cout << " underestimate capacity for 2d array\n";
			  // double the capacity of m_xy array
			     
			  double * new_xy = new double [6* m_2dcapacity];
			  int jj=0;
			  for (jj=0; jj<3*m_2dcapacity; jj++)
			    {
			      new_xy[jj] = m_xy[jj];
				
			    }
			  for (jj=3*m_2dcapacity-1; jj<6*m_2dcapacity; jj+=3)
			    new_xy[jj] = 0.;
			  // make 0 the z coordinate
			  m_2dcapacity *= 2;
			  delete [] m_xy;
			  m_xy = new_xy;
			}
		      m_xy[3*m_num2dPoints] = pp[0];
		      m_xy[3*m_num2dPoints+1] = pp[1];
		      m_num2dPoints++;
		      if (dbg)
			{
			  std::cout<< " new 2d " << m_num2dPoints - 1 << " : " << pp[0] << " " << pp[1] <<  "on edge " << ae << "\n";
			}
		      found = 1;
		    }
		} 
	    }
        }
      if (!found)
        {
	  std::cout << " a point pp is not on a red triangle " << *pp << " " << pp[1] << " red triangle " << red << " \n";
	  exit (1);
        }
    }
  // now we can build the triangles, from P array, with foundIds
  if (nP>=3)
    {
      // what are the triangles 
      FinalTriangle ftr;
      ftr.redTriangle = red;
      ftr.blueTriangle = blue;
      ftr.v[0] = foundIds[0];
      for (int i=1; i<=nP-2; i++)
	{
	  // triangle 0, i, i+1
	  ftr.v[1]=foundIds[i];
	  ftr.v[2] = foundIds[i+1];
	  m_finalMesh.push_back(ftr);
	  if (dbg)
	    {
	      std::cout << " triangle " << ftr.v[0] << " " << ftr.v[1] << " " << ftr.v[2] << "\n"; 
	    }
	}
    }
  delete [] foundIds;
  foundIds = NULL;
  return 0;
}


int ProjectShell::writeNewMesh(iMesh_Instance mesh)
{
  // here we will create new vertices, and use the coordinates in 2D
  //  m_num2dPoints is the number of nodes (some are not used, because they were probably   // collapsed
  //  we do not specifically merge red and blue nodes, but we are looking first for 
  //  nodes in red , then blue, then on red edges; some will not appear, according
  //  to the tolerance
  //
  iBase_EntityHandle * newVerts = NULL; // no vertices yet
  // iBase_INTERLEAVED
  int err= 0;
  int size1, size2;
  iMesh_createVtxArr(mesh,
		     /*in*/ m_num2dPoints,
		     /*in*/ iBase_INTERLEAVED,
		     /*in*/ m_xy,
		     /*in*/ m_num2dPoints*3,
		     /*inout*/ &newVerts,
		     /*inout*/ &size1,
		     /*inout*/ &size2,
		     /*out*/ &err); 
  if (err!=0)
    {
      std::cout<<"can't create vertices\n";
      exit (1);
    }
  // size of 
  //  then entity set, then elements (triangles)
  //
  int numTriangles = m_finalMesh.size();
  long int * adjacency = new long int [3*numTriangles];
  iBase_EntityHandle * conn = (iBase_EntityHandle *)adjacency;
  for (int L =0; L<numTriangles; L++)
    {
      FinalTriangle & tria = m_finalMesh[L];
 
      
      for (int k=0; k<3; k++)
	{
          int indexInV =  tria.v[k];
	  conn[3*L+k] = newVerts[indexInV];
	}
    }
  int numElements = numTriangles;
  iBase_EntitySetHandle orig_set;
  iMesh_createEntSet(mesh, 0, &orig_set, &err);
  int n = numTriangles;
  int junk1 = n, junk2 = n, junk3 = n, junk4 = n;
  int * stat = new int [numElements];
  int* ptr2 = stat;
  int ierr;
  iBase_EntityHandle *   	 new_entity_handles = NULL;
  iMesh_createEntArr( mesh,
		      iMesh_TRIANGLE,
		      conn, 3*numElements,
		      &new_entity_handles, &junk1, &junk2,
		      &ptr2, &junk3, &junk4,
		      &ierr );
  if (ierr!=0)
    {
      std::cout<<" can't create triangles\n";
      exit (1);
    }
  iMesh_addEntArrToSet  	(  mesh,
				   new_entity_handles,
				   numElements,
				   orig_set,
				   &ierr);
  if (ierr!=0)
    {
      std::cout<< " can't add to entity set \n";
      exit(1);
    }	 
  // now, look at the tags : red is positive, blue is negative oriented triangle
  iBase_TagHandle   	pos_tag_handle; 
  
  const char * tagName1 = "Positive";
  iMesh_createTag  	( mesh,
		          tagName1,
		          /*  size ? */ 1,
		          iBase_INTEGER ,
		          &pos_tag_handle,
		          &ierr,
		          std::strlen(tagName1) ) ;	 
	
  if (ierr!=0)
    {
      std::cout<< " can't create positive tag \n";
      exit(1);
    }
  int * tagValues = new int [numElements];
  for (int e=0; e<numElements; e++)
    {
      int redId = m_finalMesh[e].redTriangle;
      tagValues[e] = m_redMesh[ redId ].oldId;
       
    }
  // positive or negative triangles
  iMesh_setIntArrData  	(  mesh,
			   new_entity_handles,
			   numElements,
			   pos_tag_handle,
			   tagValues,
			   numElements,
			   &ierr); 
  if (ierr!=0)
    {
      std::cout<< " can't create positive tag field \n";
      exit(1);
    }
  // now blue tag
  // now, look at the tags : red is positive, blue is negative oriented triangle
  iBase_TagHandle   	neg_tag_handle; 
  
  const char * tagName2 = "Negative";
  iMesh_createTag  	( mesh,
		          tagName2,
		          /*  size ? */ 1,
		          iBase_INTEGER ,
		          &neg_tag_handle,
		          &ierr,
		          std::strlen(tagName1) ) ;	 
	
  if (ierr!=0)
    {
      std::cout<< " can't create negative tag \n";
      exit(1);
    }
  for (int e=0; e<numElements; e++)
    {
      int blueId = m_finalMesh[e].blueTriangle;
      tagValues[e] = m_blueMesh[ blueId ].oldId;
       
    }
  // positive or negative triangles
  iMesh_setIntArrData  	(  mesh,
			   new_entity_handles,
			   numElements,
			   neg_tag_handle,
			   tagValues,
			   numElements,
			   &ierr); 
   if (ierr!=0)
    {
      std::cout<< " can't create negative tag field \n";
      exit(1);
    }
   delete [] tagValues;
   delete [] adjacency;
  return 0;
}








