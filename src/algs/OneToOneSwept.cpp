#include "meshkit/OneToOneSwept.hpp"
#include "meshkit/MKCore.hpp"
#include "meshkit/iMesh.hh"
#include "meshkit/RegisterMeshOp.hpp"
#include <iostream>
#include <math.h>
#include <map>


#define Sign(u, v) ( (v)>=0.0 ? Abs(u) : -Abs(u) )
#define Max(u, v) ( (u)>(v)? (u) : (v) )
#define Abs(u) ((u)>0 ? (u) : (-u))
#define Min(u, v) ( (u)>(v)? (v) : (u) )


namespace MeshKit
{

//---------------------------------------------------------------------------//
//Entity Type initilization for OneToOneSwept meshing
moab::EntityType OneToOneSwept_tps[] = {moab::MBVERTEX, moab::MBQUAD, moab::MBHEX};
iBase_EntityType OneToOneSwept_mtp = iBase_REGION;

// static registration of this  mesh scheme
static RegisterMeshOp<OneToOneSwept,true> INIT("OneToOneSwept", &OneToOneSwept_mtp, 1, 
                                            OneToOneSwept_tps,
                                            sizeof(OneToOneSwept_tps)/sizeof(OneToOneSwept_tps[0]) );

//---------------------------------------------------------------------------//
// make an instance of the OneToOneSwept class
MeshOp *OneToOneSwept::factory(MKCore *mkcore, const MEntVector &me_vec)
{
	return new OneToOneSwept(mkcore, me_vec);
}

//---------------------------------------------------------------------------//
// construction function for OneToOneSwept class
OneToOneSwept::OneToOneSwept(MKCore *mk_core, const MEntVector &me_vec) : MeshScheme(mk_core, me_vec)
{

}

//---------------------------------------------------------------------------//
// deconstruction function for OneToOneSwept class
OneToOneSwept::~OneToOneSwept()
{

}


//---------------------------------------------------------------------------//
// return the type of entities this mesh creates
void OneToOneSwept::mesh_types(std::vector<moab::EntityType> &tps)
{
    tps.push_back(moab::MBVERTEX);
    tps.push_back(moab::MBHEX);
}

//---------------------------------------------------------------------------//
// setup function
void OneToOneSwept::setup_this()
{

}

//---------------------------------------------------------------------------//
// execute function: generate the all-hex mesh through sweeping from source 
// surface to target surface
void OneToOneSwept::execute_this()
{
	std::vector<double> coords;
	std::vector<moab::EntityHandle> nodes;

}

//---------------------------------------------------------------------------//
// set the source surface function
void OneToOneSwept::SetSourceSurface()
{

}

//---------------------------------------------------------------------------//
// set the target surface function
void OneToOneSwept::SetTargetSurface()
{

}

//---------------------------------------------------------------------------//
//function for obtaining the parametric coordinates from x,y,z coordinates
int OneToOneSwept::getUVCoords(iBase_EntityHandle gFaceHandle, Point3D pts3, Point2D &pts2)
{
	double xmin, ymin, zmin, xmax, ymax, zmax;

	iGeom::Error g_err = mk_core()->igeom_instance()->getEntBoundBox(gFaceHandle, xmin, ymin, zmin, xmax, ymax, zmax);
	IBERRCHK(g_err, "Trouble get the bounding box for the face entity.");
	if (pts3.px < xmin || pts3.px > xmax)
	{
        	cout << "Warning: Query point outside X Range [" << xmin << "," << xmax << "], x=" << pts3.px << endl;
	}
    	if (pts3.py < ymin || pts3.py > ymax)
        {
		cout << "Warning: Query point outside Y Range [" << ymin << "," << ymax << "], y=" << pts3.py << endl;
	}
    	if (pts3.pz < zmin || pts3.pz > zmax)
	{        
		cout << "Warning: Query point outside Z Range [" << zmin << "," << zmax << "], z=" << pts3.pz << endl;
	}
	g_err = mk_core()->igeom_instance()->getEntXYZtoUV(gFaceHandle, pts3.px, pts3.py, pts3.pz, pts2.pu, pts2.pv);
	IBERRCHK(g_err, "Trouble get the parametric coordinates from x,y,z coordinates.");

	return 1;
}

//---------------------------------------------------------------------------//
//function for obtaining the x,y,z coordinates from parametric coordinates
int OneToOneSwept::getXYZCoords(iBase_EntityHandle gFaceHandle, Point3D &pts3, double uv[2])
{
	double umin, umax, vmin, vmax;
	iGeom::Error g_err = mk_core()->igeom_instance()->getEntUVRange(gFaceHandle, umin, vmin, umax, vmax);
	IBERRCHK(g_err, "Trouble get the parametric coordinate range.");	
	if ((uv[0]<umin)||(uv[0]>umax))
		cout << "Warning: U exceeds the range" << endl;
	if ((uv[1]<vmin)||(uv[1]>vmax))
		cout << "Warning: V exceeds the range" << endl;
	g_err = mk_core()->igeom_instance()->getEntUVtoXYZ(gFaceHandle, uv[0], uv[1], pts3.px, pts3.py, pts3.pz);
	IBERRCHK(g_err, "Trouble get the x,y,z coordinates from parametric coordinates.");	

	return 1;
}

//---------------------------------------------------------------------------//
// map the mesh on the source surface to the target surface
int OneToOneSwept::TargetSurfProjection()
{
	iBase_EntitySetHandle targetSet;
	iRel::Error r_err = mk_core()->irel_pair()->getEntSetRelation(targetSurface, 0, targetSet);	
	if (!r_err)
	{//there exists an mesh on the target surfaces
		return 1;
	}
	
	iBase_TagHandle taghandle;
	iMesh::Error m_err = mk_core()->imesh_instance()->getTagHandle("source", taghandle);
	IBERRCHK(m_err, "Trouble get the tag handle 'source'.");
	
	int index=0, id, index_t;
	double sumin, sumax, tumin, tumax, su, tu;
	std::vector<iBase_EntityHandle> gsEdge, gtEdge, edgeNodes, geomEdgeEnds, meshEdgeEnds;
	std::vector<iBase_EntitySetHandle> edgeEndsMeshSets;
	vector<Point2D> sPtsUV(0), tPtsUV(0);   //store the boundary nodes on the sourface surface and target surface
	Point3D pts;
	Point2D Sc={0.0,0.0}, Tc={0.0,0.0}; //modified coordinates in order not to be singular matrix

	//define an array variable for storing the affine mapping matrix
	double A[2][2]= {{0, 0}, {0, 0}}; //affine map matrix
	double temp[2][2] = {{0, 0}, {0, 0}};//define an array variable for storing the sum of boundary nodes
	double b1[2] = {0,0}, b2[2] = {0,0};
		
	///////////////////////////////////////////////////////////////////////////////
    	// Step II: Collect all the geometric edges and discretize them.            //
    	///////////////////////////////////////////////////////////////////////////////

	
	//define the mesh set for various edges on the source surface
	std::vector<iBase_EntitySetHandle> edgeMeshSets(gsEdgeList.size());
	//loop over the various edges
	
	for (int i=0; i < gsEdgeList.size(); i++)
	{
		//get the mesh entityset for edge[i]
		r_err = mk_core()->irel_pair()->getEntSetRelation(gsEdgeList[i].gEdgeHandle, 0, edgeMeshSets[i]);
		IBERRCHK(r_err, "Trouble get the tag handle 'source'.");		
		
		//get the edge nodes for edge[i] mesh
		edgeNodes.clear();
		m_err = mk_core()->imesh_instance()->getEntities(edgeMeshSets[i], iBase_VERTEX, iMesh_POINT, edgeNodes);
		IBERRCHK(m_err, "Trouble get the mesh edge entities.");


		//find the corresponding relationship for both ends on the source suface and target surface
		//loop over the linking sides
		int index_a, index_b; //record the vertex numbering on the target surface
		index_a = cornerPairs[gsEdgeList[i].connect[0]->index];
		index_b = cornerPairs[gsEdgeList[i].connect[1]->index];
		
		//create the corresponding edge relationship between the source surface and target surface
		index_t = edgePairs[i];
		
		//detect whether there exists an mesh for the edge on the target surface
		iBase_EntitySetHandle entityset;
		r_err = mk_core()->irel_pair()->getEntSetRelation(gtEdgeList[index_t].gEdgeHandle, 0, entityset);		
		if (!r_err)
		{//there exists an mesh on the boundary of target surface
			std::vector<iBase_EntityHandle> tnodes;
			m_err = mk_core()->imesh_instance()->getEntities(entityset, iBase_VERTEX, iMesh_POINT, tnodes);
			IBERRCHK(m_err, "Trouble get the mesh entities.");

			//determine whether there are the same number of nodes between two corresponding edges
			if (edgeNodes.size()!=tnodes.size())
			{
				cout << "the node number between two corresponding edges is different" << endl;
				exit(1);
			}
			
			//determine the cooresponding relationship of mesh nodes between two corresponding edges
			bool isReverse1=false, isReverse2=false;
			double xyz1[3], xyz2[3], dist1, dist2;
			if (edgeNodes.size() > 1)
			{
				m_err = mk_core()->imesh_instance()->getVtxCoord(edgeNodes[0], xyz1[0], xyz1[1], xyz1[2]);
				IBERRCHK(m_err, "Trouble get the mesh node coordinates.");

				m_err = mk_core()->imesh_instance()->getVtxCoord(edgeNodes[1], xyz2[0], xyz2[1], xyz2[2]);
				IBERRCHK(m_err, "Trouble get the mesh node coordinates.");
				
				dist1 = sqrt( pow(gsEdgeList[i].connect[0]->xyzCoords[0]-xyz1[0],2) + pow(gsEdgeList[i].connect[0]->xyzCoords[1]-xyz1[1],2) + pow(gsEdgeList[i].connect[0]->xyzCoords[2]-xyz1[2],2));
				dist2 = sqrt( pow(gsEdgeList[i].connect[0]->xyzCoords[0]-xyz2[0],2) + pow(gsEdgeList[i].connect[0]->xyzCoords[1]-xyz2[1],2) + pow(gsEdgeList[i].connect[0]->xyzCoords[2]-xyz2[2],2));
				if (dist1 >= dist2)
				{
					isReverse1 = true;
				}
				
				m_err = mk_core()->imesh_instance()->getVtxCoord(tnodes[0], xyz1[0], xyz1[1], xyz1[2]);
				IBERRCHK(m_err, "Trouble get the mesh node coordinates.");
				m_err = mk_core()->imesh_instance()->getVtxCoord(tnodes[1], xyz2[0], xyz2[1], xyz2[2]);
				IBERRCHK(m_err, "Trouble get the mesh node coordinates.");
				
				dist1 = sqrt( pow(gVertexList[index_a].xyzCoords[0]-xyz1[0], 2) + pow(gVertexList[index_a].xyzCoords[1]-xyz1[1], 2) + pow(gVertexList[index_a].xyzCoords[2]-xyz1[2], 2));
				dist2 = sqrt( pow(gVertexList[index_a].xyzCoords[0]-xyz2[0], 2) + pow(gVertexList[index_a].xyzCoords[1]-xyz2[1], 2) + pow(gVertexList[index_a].xyzCoords[2]-xyz2[2], 2));
				if (dist1 >= dist2)
				{
					isReverse2 = true;
				}
			}
			for (int j=0; j < edgeNodes.size(); j++)
			{
				int tmpIndex;
				m_err = mk_core()->imesh_instance()->getIntData(edgeNodes[j], taghandle, tmpIndex);
				IBERRCHK(m_err, "Trouble get the int data of edge nodes.");
				
				index++;
				
				sPtsUV.resize(index);
				tPtsUV.resize(index);
				
				pts.px = NodeList[tmpIndex].xyzCoords[0];
				pts.py = NodeList[tmpIndex].xyzCoords[1];
				pts.pz = NodeList[tmpIndex].xyzCoords[2];
				
				getUVCoords(sourceSurface, pts, sPtsUV[index-1]);
				
				if (isReverse1!=isReverse2)
				{
					m_err = mk_core()->imesh_instance()->getVtxCoord(tnodes[edgeNodes.size()-j-1], TVertexList[tmpIndex].xyzCoords[0], TVertexList[tmpIndex].xyzCoords[1], TVertexList[tmpIndex].xyzCoords[2]);
					IBERRCHK(m_err, "Trouble get the edge node coordinates.");
					
					pts.px = TVertexList[tmpIndex].xyzCoords[0];
					pts.py = TVertexList[tmpIndex].xyzCoords[1];
					pts.pz = TVertexList[tmpIndex].xyzCoords[2];
					
					getUVCoords(targetSurface, pts, tPtsUV[index-1]);
					
					TVertexList[tmpIndex].uvCoords[0] = tPtsUV[index-1].pu;
					TVertexList[tmpIndex].uvCoords[1] = tPtsUV[index-1].pv;
					TVertexList[tmpIndex].gVertexHandle = tnodes[edgeNodes.size()-j-1];
					TVertexList[tmpIndex].onCorner = false;
					TVertexList[tmpIndex].onBoundary = true;
				}
				else
				{
					m_err = mk_core()->imesh_instance()->getVtxCoord(tnodes[j], TVertexList[tmpIndex].xyzCoords[0], TVertexList[tmpIndex].xyzCoords[1], TVertexList[tmpIndex].xyzCoords[2]);
					IBERRCHK(m_err, "Trouble get the edge node coordinates.");
					
					pts.px = TVertexList[tmpIndex].xyzCoords[0];
					pts.py = TVertexList[tmpIndex].xyzCoords[1];
					pts.pz = TVertexList[tmpIndex].xyzCoords[2];
					
					getUVCoords(targetSurface, pts, tPtsUV[index-1]);
					
					TVertexList[tmpIndex].uvCoords[0] = tPtsUV[index-1].pu;
					TVertexList[tmpIndex].uvCoords[1] = tPtsUV[index-1].pv;
					TVertexList[tmpIndex].gVertexHandle = tnodes[j];
					TVertexList[tmpIndex].onCorner = false;
					TVertexList[tmpIndex].onBoundary = true;
				}				
			}		
			
			continue;		
		}
		
		//get the parametric variable u for both ends on the boundary edge
		double sLeft, sRight, tLeft, tRight;
		iGeom::Error g_err = mk_core()->igeom_instance()->getEntXYZtoU(gsEdgeList[i].gEdgeHandle, gsEdgeList[i].connect[0]->xyzCoords[0], gsEdgeList[i].connect[0]->xyzCoords[1], gsEdgeList[i].connect[0]->xyzCoords[2], sLeft);
		IBERRCHK(g_err, "Trouble get the parametric coordinates from x,y,z coordinates.");

		g_err = mk_core()->igeom_instance()->getEntXYZtoU(gsEdgeList[i].gEdgeHandle, gsEdgeList[i].connect[1]->xyzCoords[0], gsEdgeList[i].connect[1]->xyzCoords[1], gsEdgeList[i].connect[1]->xyzCoords[2], sRight);
		IBERRCHK(g_err, "Trouble get the parametric coordinates from x,y,z coordinates.");

		g_err = mk_core()->igeom_instance()->getEntXYZtoU(gtEdgeList[index_t].gEdgeHandle, gVertexList[index_a].xyzCoords[0], gVertexList[index_a].xyzCoords[1], gVertexList[index_a].xyzCoords[2], tLeft);
		IBERRCHK(g_err, "Trouble get the parametric coordinates from x,y,z coordinates.");
		
		g_err = mk_core()->igeom_instance()->getEntXYZtoU(gtEdgeList[index_t].gEdgeHandle, gVertexList[index_b].xyzCoords[0], gVertexList[index_b].xyzCoords[1], gVertexList[index_b].xyzCoords[2], tRight);
		IBERRCHK(g_err, "Trouble get the parametric coordinates from x,y,z coordinates.");

		std::vector<iBase_EntityHandle> newNodehandle(edgeNodes.size());
		std::vector<iBase_EntityHandle> testNodeHandle;
		Point3D testPoint;
		//create the boundary node on the edges. This doesn't include the corners.
		for (int j=0; j < edgeNodes.size(); j++)
		{
			//get the cartesian coordinates for the edge nodes
			m_err = mk_core()->imesh_instance()->getVtxCoord(edgeNodes[j], pts.px, pts.py, pts.pz);
			IBERRCHK(g_err, "Trouble get the x,y,z coordinates.");

			//get the parametric coordinates for the vertex on the source surface
			index++;
			sPtsUV.resize(index);
			getUVCoords(sourceSurface, pts, sPtsUV[index-1]);

			//get vertex id for target surface vertex
			m_err = mk_core()->imesh_instance()->getIntData(edgeNodes[j], taghandle, id);
			IBERRCHK(m_err, "Trouble get the vertex id.");
			
			TVertexList[id].index = id;		
						
			//transform the cartesian coordinates into the parametric coordinates
			g_err = mk_core()->igeom_instance()->getEntXYZtoU(gsEdgeList[i].gEdgeHandle, pts.px, pts.py, pts.pz, su);
			IBERRCHK(g_err, "Trouble get the cartesian coordinates x, y, z.");

			//calculate the parametric coordinates on the target surface
			tu = tLeft + (tRight - tLeft) * (su - sLeft) / (sRight - sLeft);

			//transform the parametric coordinates into the cartesian coordinates on the target surface
			g_err = mk_core()->igeom_instance()->getEntUtoXYZ(gtEdgeList[index_t].gEdgeHandle, tu, pts.px, pts.py, pts.pz);
			IBERRCHK(g_err, "Trouble get the parametric coordinates from the cartesian coordinates x, y, z.");

			//create the vertex entity on the target boundary edge
			m_err = mk_core()->imesh_instance()->createVtx(pts.px, pts.py, pts.pz, newNodehandle[j]);
			IBERRCHK(m_err, "Trouble create a new node mesh entity.");

			//Get the parametric coordinates for new generated vertex on the target surface
			tPtsUV.resize(index);
			getUVCoords(targetSurface, pts, tPtsUV[index-1]);
			
			//add the new generated vertex on the target surface into the list
			TVertexList[id].xyzCoords[0] = pts.px; 
			TVertexList[id].xyzCoords[1] = pts.py;
			TVertexList[id].xyzCoords[2] = pts.pz;
			TVertexList[id].uvCoords[0] = tPtsUV[index-1].pu;
			TVertexList[id].uvCoords[1] = tPtsUV[index-1].pv;
			TVertexList[id].gVertexHandle = newNodehandle[j];
			TVertexList[id].onCorner = false;
			TVertexList[id].onBoundary = true;
			
		}
		
		//create entityset for storing the nodes on the target boundary edge and build association between the geometry and mesh
		r_err = mk_core()->irel_pair()->getEntSetRelation(gtEdgeList[index_t].gEdgeHandle, 0, entityset);
		if (r_err) //there is no entityset associated with gtEdgeList[index_t].gEdgeHandle
		{
			m_err = mk_core()->imesh_instance()->createEntSet(1, entityset);
			IBERRCHK(m_err, "Trouble create a new mesh entity set.");
		}
		
		//build the association
		m_err = mk_core()->imesh_instance()->addEntArrToSet(&newNodehandle[0], edgeNodes.size(), entityset);
		IBERRCHK(m_err, "Trouble add an array of nodes into the entity set.");
		
		//create the line segments on the boundary edge of target surface
		std::vector<iBase_EntityHandle> sedgeHandle;
		std::vector<iBase_EntityHandle> tedgeHandle(0);
		int edgeIndex=0;
		m_err = mk_core()->imesh_instance()->getEntities(edgeMeshSets[i], iBase_EDGE, iMesh_LINE_SEGMENT, sedgeHandle);
		IBERRCHK(m_err, "Trouble get the mesh line segments.");	
	
		for (int j=0; j< sedgeHandle.size(); j++)
		{
			int status;
			int nodeindex1, nodeindex2;
			std::vector<iBase_EntityHandle> connect(2);
			std::vector<iBase_EntityHandle> tmpNodes;
			m_err = mk_core()->imesh_instance()->getEntAdj(sedgeHandle[j], iBase_VERTEX, tmpNodes);
			IBERRCHK(m_err, "Trouble get the adjacent nodes of the edge.");

			assert(tmpNodes.size()==2);
			m_err = mk_core()->imesh_instance()->getIntData(tmpNodes[0], taghandle, nodeindex1);
			IBERRCHK(m_err, "Trouble get the int data for nodes.");
			m_err = mk_core()->imesh_instance()->getIntData(tmpNodes[1], taghandle, nodeindex2);
			IBERRCHK(m_err, "Trouble get the int data for nodes.");
			
			edgeIndex++;
			tedgeHandle.resize(edgeIndex);
			connect[0] = TVertexList[nodeindex1].gVertexHandle;
			connect[1] = TVertexList[nodeindex2].gVertexHandle;
			
			m_err = mk_core()->imesh_instance()->createEnt(iMesh_LINE_SEGMENT, &connect[0], 2, tedgeHandle[edgeIndex-1]);
			IBERRCHK(m_err, "Trouble create the line segment entities.");			
		}
		m_err = mk_core()->imesh_instance()->addEntArrToSet(&tedgeHandle[0], tedgeHandle.size(), entityset);
		IBERRCHK(m_err, "Trouble add an array of line segments to the entity set.");


		r_err = mk_core()->irel_pair()->getEntSetRelation(gtEdgeList[index_t].gEdgeHandle, 0, entityset);
		if (r_err) //there is no entityset associated with gtEdgeList[index_t].gEdgeHandle
		{
			r_err = mk_core()->irel_pair()->setEntSetRelation(gtEdgeList[index_t].gEdgeHandle, entityset);
			IBERRCHK(r_err, "Trouble set the association between the geometry and entityset.");
		}								
	}
	
	//Until now, all the nodes have been created on the boundary edge.
	//get the corner coordinates
	assert(NodeList.size()==TVertexList.size());
	for (int i=0; i < NodeList.size(); i++)
	{
		if (NodeList[i].onCorner||NodeList[i].onBoundary)
		{
			Point3D pts;
			index++;

			sPtsUV.resize(index);
			pts.px = NodeList[i].xyzCoords[0];
			pts.py = NodeList[i].xyzCoords[1];
			pts.pz = NodeList[i].xyzCoords[2];
			getUVCoords(sourceSurface, pts, sPtsUV[index-1]);
			
			pts.px = TVertexList[i].xyzCoords[0];
			pts.py = TVertexList[i].xyzCoords[1];
			pts.pz = TVertexList[i].xyzCoords[2];
			tPtsUV.resize(index);
			getUVCoords(targetSurface, pts, tPtsUV[index-1]);		
		}
	}
	
	
	//create the A matrix for mapping (affine) the inner nodes.
	//this affine mapping is done in the parametric domain
	//first get the modified coefficents
	for (int i=0; i < index; i++)
	{
		Sc.pu = Sc.pu + sPtsUV[i].pu;
		Sc.pv = Sc.pv + sPtsUV[i].pv;
		Tc.pu = Tc.pu + tPtsUV[i].pu;
		Tc.pv = Tc.pv + tPtsUV[i].pv;

	}

	Sc.pu = Sc.pu/double(index);
	Sc.pv = Sc.pv/double(index);
	Tc.pu = Tc.pu/double(index);
	Tc.pv = Tc.pv/double(index);

	//get the new coordinates for u&v in the parametric domain	
	for (int i=0; i < index; i++)
	{
		sPtsUV[i].pu = sPtsUV[i].pu - Sc.pu;
		sPtsUV[i].pv = sPtsUV[i].pv - Sc.pv;
		tPtsUV[i].pu = tPtsUV[i].pu - Tc.pu;
		tPtsUV[i].pv = tPtsUV[i].pv - Tc.pv;
	}
	
	//get sum of boundary nodes' coordinate
	for (int i=0; i < index; i++)
	{
		temp[0][0] = temp[0][0] + sPtsUV[i].pu*sPtsUV[i].pu;
		temp[0][1] = temp[0][1] + sPtsUV[i].pu*sPtsUV[i].pv;
		temp[1][0] = temp[1][0] + sPtsUV[i].pu*sPtsUV[i].pv;
		temp[1][1] = temp[1][1] + sPtsUV[i].pv*sPtsUV[i].pv;

		b1[0] = b1[0] + sPtsUV[i].pu*tPtsUV[i].pu;
		b1[1] = b1[1] + sPtsUV[i].pv*tPtsUV[i].pu;
		b2[0] = b2[0] + sPtsUV[i].pu*tPtsUV[i].pv;
		b2[1] = b2[1] + sPtsUV[i].pv*tPtsUV[i].pv;
	}

	//Solve the equation to get affine mapping matrix A.resize
	//check whether the equation has the solution
	assert((temp[0][0]*temp[1][1]-temp[0][1]*temp[1][0])!=0);
	A[0][0] = (temp[1][1]*b1[0] - temp[0][1]*b1[1])/(temp[0][0]*temp[1][1]-temp[0][1]*temp[1][0]);
	A[0][1] = (temp[0][0]*b1[1]-temp[1][0]*b1[0])/(temp[0][0]*temp[1][1]-temp[0][1]*temp[1][0]);
	
	A[1][0] = (temp[1][1]*b2[0] - temp[0][1]*b2[1])/(temp[0][0]*temp[1][1]-temp[0][1]*temp[1][0]);
	A[1][1] = (temp[0][0]*b2[1]-temp[1][0]*b2[0])/(temp[0][0]*temp[1][1]-temp[0][1]*temp[1][0]);
	
	//the affine mapping matrix A is obtained

	//mapping the inner nodes on the source surface onto the target surface
	iBase_EntitySetHandle entityset;  //this entityset is for storing the inner nodes on the target surface
	vector<iBase_EntityHandle>  newNodehandle(0), newEdgeHandle(0);
	int newIndex=0;

	r_err = mk_core()->irel_pair()->getEntSetRelation(targetSurface, 0, entityset);
	if (r_err) //there is no entityset associated with targetSurface
	{
		m_err = mk_core()->imesh_instance()->createEntSet(1, entityset);
		IBERRCHK(m_err, "Trouble create the entity set");
	}

	//create the inner nodes on the target surface
	for (int i=0; i < NodeList.size(); i++)
	{
		if ((!NodeList[i].onBoundary)&&(!NodeList[i].onCorner))//make sure that the node is the inner node
		{
			//first transform the vertex on the source surface: affine mapping
			double uv[2];
			Point3D pts; 
			uv[0] = A[0][0]*(NodeList[i].uvCoords[0] - Sc.pu) + A[0][1]*(NodeList[i].uvCoords[1]- Sc.pv) + Tc.pu;
			uv[1] = A[1][0]*(NodeList[i].uvCoords[0]- Sc.pu) + A[1][1]*(NodeList[i].uvCoords[1]- Sc.pv) + Tc.pv;
			
			getXYZCoords(targetSurface, pts, uv);//maybe there is a problem here,   notice
			//create the vertex on the target surface
			
			newIndex++;
			newNodehandle.resize(newIndex);
			m_err = mk_core()->imesh_instance()->createVtx(pts.px, pts.py, pts.pz, newNodehandle[newIndex-1]);
			IBERRCHK(m_err, "Trouble create the vertex entity.");

			//add the new generated vertex on the target surface into the list
			TVertexList[i].xyzCoords[0] = pts.px;
			TVertexList[i].xyzCoords[1] = pts.py;
			TVertexList[i].xyzCoords[2] = pts.pz;
			TVertexList[i].uvCoords[0] = uv[0];
			TVertexList[i].uvCoords[1] = uv[1];
			TVertexList[i].gVertexHandle = newNodehandle[newIndex-1];
			TVertexList[i].onCorner = false;
			TVertexList[i].onBoundary = false;
			TVertexList[i].index = NodeList[i].index;	
		}
	}
	//add the inner nodes to the entityset
	m_err = mk_core()->imesh_instance()->addEntArrToSet(&newNodehandle[0], newIndex, entityset);
	IBERRCHK(m_err, "Trouble add an array of nodes to the entityset.");

	//until now, all the nodes have been generated on the target surface

	//create the edge elements on the target surface
	int status;
	index =0;
	vector<iBase_EntityHandle> gEdgeHandle(0);
	for (int i=0; i < EdgeList.size(); i++)
	{
		if (!EdgeList[i].onBoundary)
		{
			//create the inner edge elements
			vector<iBase_EntityHandle> connect(2);
			connect[0] = TVertexList[EdgeList[i].connect[0]->index].gVertexHandle;
			connect[1] = TVertexList[EdgeList[i].connect[1]->index].gVertexHandle;
			
			index++;			
			gEdgeHandle.resize(index);

			m_err = mk_core()->imesh_instance()->createEnt(iMesh_LINE_SEGMENT, &connect[0], 2, gEdgeHandle[index-1]);
			IBERRCHK(m_err, "Trouble create the line segment entity.");
			
			//add the new generated inner edge elements on the target surface into the list
			TEdgeList[i].index = i;
			TEdgeList[i].connect[0] = &TVertexList[EdgeList[i].connect[0]->index];
			TEdgeList[i].connect[1] = &TVertexList[EdgeList[i].connect[1]->index];
			TEdgeList[i].gEdgeHandle = gEdgeHandle[index-1];
				
		}		
	}
	//add the inner edge elements to the entityset
	m_err = mk_core()->imesh_instance()->addEntArrToSet(&gEdgeHandle[0], index, entityset);
	IBERRCHK(m_err, "Trouble add an array of line segments to the entity set.");


	//determine the numbering order for quadrilateral nodes
	int sense_out, sense_out1, sense_out2, sense_out3, sense_out4;
	iGeom::Error g_err = mk_core()->igeom_instance()->getEgFcSense(gsEdgeList[0].gEdgeHandle, sourceSurface, sense_out1);
	IBERRCHK(g_err, "Trouble get the sense of edge with respect to the face.");
	g_err = mk_core()->igeom_instance()->getEgFcSense(gtEdgeList[edgePairs[gsEdgeList[0].index]].gEdgeHandle, targetSurface, sense_out2);
	IBERRCHK(g_err, "Trouble get the sense of edge with respect to the face.");
	g_err = mk_core()->igeom_instance()->getEgVtxSense(gsEdgeList[0].gEdgeHandle, gsEdgeList[0].connect[0]->gVertexHandle, gsEdgeList[0].connect[1]->gVertexHandle, sense_out3);
	IBERRCHK(g_err, "Trouble get the sense of vertex with respect to the edge.");	
	g_err = mk_core()->igeom_instance()->getEgVtxSense(gtEdgeList[edgePairs[gsEdgeList[0].index]].gEdgeHandle, gVertexList[cornerPairs[gsEdgeList[0].connect[0]->index]].gVertexHandle, gVertexList[cornerPairs[gsEdgeList[0].connect[1]->index]].gVertexHandle, sense_out4);
	IBERRCHK(g_err, "Trouble get the sense of vertex with respect to the edge.");	
	sense_out = sense_out1*sense_out2*sense_out3*sense_out4;


	//create the quadrilateral elements on the target surface
	vector<iBase_EntityHandle> mFaceHandle(FaceList.size());
	for (int i=0; i < FaceList.size(); i++)
	{
		vector<iBase_EntityHandle> connect(FaceList[i].getNumNodes());
		
		if (sense_out < 0)
		{
			connect[0] = TVertexList[(FaceList[i].getVertex(0))->index].gVertexHandle;
			connect[1] = TVertexList[(FaceList[i].getVertex(3))->index].gVertexHandle;
			connect[2] = TVertexList[(FaceList[i].getVertex(2))->index].gVertexHandle;
			connect[3] = TVertexList[(FaceList[i].getVertex(1))->index].gVertexHandle;
		}
		else
		{
			connect[0] = TVertexList[(FaceList[i].getVertex(0))->index].gVertexHandle;
			connect[1] = TVertexList[(FaceList[i].getVertex(1))->index].gVertexHandle;
			connect[2] = TVertexList[(FaceList[i].getVertex(2))->index].gVertexHandle;
			connect[3] = TVertexList[(FaceList[i].getVertex(3))->index].gVertexHandle;
		}
		m_err = mk_core()->imesh_instance()->createEnt(iMesh_QUADRILATERAL, &connect[0], FaceList[i].getNumNodes(), mFaceHandle[i]);
		IBERRCHK(g_err, "Trouble create the quadrilateral entity.");	
		
		//add the face elements on the target surface to the list
		TFaceList[i].index = FaceList[i].index;
		TFaceList[i].gFaceHandle = mFaceHandle[i];
		TFaceList[i].connect.resize(FaceList[i].getNumNodes());
		for (int j=0; j < FaceList[i].getNumNodes(); j++)
		{
			TFaceList[i].connect[j] = &TVertexList[FaceList[i].connect[j]->index];
		}	
	}
	//add the inner face elements to the entityset
	m_err = mk_core()->imesh_instance()->addEntArrToSet(&mFaceHandle[0], FaceList.size(), entityset);
	IBERRCHK(g_err, "Trouble add an array of quadrilateral entities to the entity set.");	
	
	//build the association
	r_err = mk_core()->irel_pair()->getEntSetRelation(targetSurface, 0, entityset);
	if (r_err) //there is no entityset associated with region[0]
	{
		r_err = mk_core()->irel_pair()->setEntSetRelation(targetSurface, entityset);
		IBERRCHK(g_err, "Trouble set the association between the target surface entity and mesh entity set.");	
	}
}


void OneToOneSwept::buildAssociation(iGeom_Instance &geom, iMesh_Instance &mesh, iRel_Instance &assoc, iRel_PairHandle &rel)
{

    int err, namelen;

    std::vector<iBase_EntitySetHandle> entitySets;
    iBase_EntitySetHandle geom_root_set, mesh_root_set;

    // Get the root sets of the geometry and mesh.
    mesh_root_set = mk_core()->igeom_instance()->getRootSet();
    geom_root_set = mk_core()->imesh_instance()->getRootSet();

    //get the global geometrical id tag, global mesh id tag, global dimension id tag
    iBase_TagHandle geom_id_tag, mesh_id_tag, geom_dim_tag;    
    iGeom::Error g_err = mk_core()->igeom_instance()->getTagHandle("GLOBAL_ID", geom_id_tag);
    IBERRCHK(g_err, "Trouble get global geometry dimension id tag.");
    iMesh::Error m_err = mk_core()->imesh_instance()->getTagHandle("GLOBAL_ID", mesh_id_tag);
    IBERRCHK(m_err, "Trouble get global mesh dimension id tag.");
    err = mk_core()->imesh_instance()->getTagHandle("GEOM_DIMENSION", geom_dim_tag);
    IBERRCHK(m_err, "Trouble get the geometric dimension id tag.");



    // Get all the entitySet in the mesh
    m_err = mk_core()->imesh_instance()->getEntSets(mesh_root_set, 0, entitySets);
    IBERRCHK(m_err, "Trouble get the geometric dimension id tag.");


    	//int ncount;
    	//iBase_EntityHandle gEntity;



    // Map all the geometric nodes
    //get size of nodes
    std::vector<iBase_EntityHandle> gNodes;
    g_err = mk_core()->igeom_instance()->getEntities(geom_root_set, iBase_VERTEX, gNodes);
    IBERRCHK(g_err, "Trouble get the geometric node entities.");
    
    int geom_id;
    std::map<int, iBase_EntityHandle> mapNodes;
    for (int i=0; i < gNodes.size(); i++)
    {
    	g_err = mk_core()->igeom_instance()->getIntData(gNodes[i], geom_id_tag, geom_id);
        IBERRCHK(g_err, "Trouble get the int data of node entities.");
        mapNodes[geom_id] = gNodes[i];
    }
    

    // Map all the geometric edges.resize
    std::vector<iBase_EntityHandle> gEdges;
    g_err = mk_core()->igeom_instance()->getEntities(geom_root_set, iBase_EDGE, gEdges);
    IBERRCHK(g_err, "Trouble get the geometric edge entities.");

    
    std::map<int, iBase_EntityHandle> mapEdges;
    for (int i = 0; i < gEdges.size(); i++)
    {
        g_err = mk_core()->igeom_instance()->getIntData(gEdges[i], geom_id_tag, geom_id);
        IBERRCHK(g_err, "Trouble get the int data of edge entities.");
        mapEdges[geom_id] = gEdges[i];
    }

    // Map all the geometric faces ...
    std::vector<iBase_EntityHandle> gFaces;
    g_err = mk_core()->igeom_instance()->getEntities(geom_root_set, iBase_FACE, gFaces);
    IBERRCHK(g_err, "Trouble get the geometric face entities.");

    std::map<int, iBase_EntityHandle> mapFaces;
    for (int i = 0; i < gFaces.size(); i++)
    {
        g_err = mk_core()->igeom_instance()->getIntData(gFaces[i], geom_id_tag, geom_id);
        IBERRCHK(g_err, "Trouble get the int data of face entities.");
        mapFaces[geom_id] = gFaces[i];
    }

    // Map all the geometric cells ...
    std::vector<iBase_EntityHandle> gCells;
    g_err = mk_core()->igeom_instance()->getEntities(geom_root_set, iBase_REGION, gCells);
    IBERRCHK(g_err, "Trouble get the geometric cell entities.");

    std::map<int, iBase_EntityHandle> mapCells;
    for (int i = 0; i < gCells.size(); i++)
    {
        g_err = mk_core()->igeom_instance()->getIntData(gCells[i], geom_id_tag, geom_id);
        IBERRCHK(g_err, "Trouble get the int data of cell entities.");
        mapCells[geom_id] = gCells[i];
    }
    
    ///////////////////////////////////////////////////////////////////////////////
    // Create Vertex Assocations:
    ///////////////////////////////////////////////////////////////////////////////
    cout << " Building Vertex Associations " << endl;
    int numNodes, ncount = 0, geom_dim;
    iBase_EntityHandle gEntity;
    m_err = mk_core()->imesh_instance()->getNumOfType(mesh_root_set, iBase_VERTEX, numNodes);
    IBERRCHK(m_err, "Trouble get the number of vertices.");
    std::vector<iBase_EntityHandle> mNodes;
    
    int numAssociations = 0;
    
    for (int i = 0; i < entitySets.size(); i++)
    {
        mNodes.clear();
        m_err = mk_core()->imesh_instance()->getEntities(entitySets[i], iBase_VERTEX, iMesh_ALL_TOPOLOGIES, mNodes);
        IBERRCHK(m_err, "Trouble get the number of vertices.");	
        if (mNodes.size() && (mNodes.size() != numNodes))
        {
            ncount += mNodes.size();

            m_err = mk_core()->imesh_instance()->getEntSetIntData(entitySets[i], mesh_id_tag, geom_id);
            IBERRCHK(m_err, "Trouble set the mesh id tag with int value for entity set.");
	
            m_err = mk_core()->imesh_instance()->getEntSetIntData(entitySets[i], geom_dim_tag, geom_dim);
            IBERRCHK(m_err, "Trouble set the geom dim tag with int value for mesh entity set.");

            gEntity = 0;
            switch (geom_dim)
            {
            case 0:
            	if (mapNodes.find(geom_id) != mapNodes.end())
                {
                    gEntity = mapNodes[geom_id];
                    numAssociations++;
                }
                else
                {
                    cout << "Fatal Error: Geometric Edge not found : " << geom_id << endl;
                    exit(0);
                }
                break;
            	
            case 1:

                if (mapEdges.find(geom_id) != mapEdges.end())
                {
                    gEntity = mapEdges[geom_id];
                    numAssociations++;
                }
                else
                {
                    cout << "Fatal Error: Geometric Edge not found : " << geom_id << endl;
                    exit(0);
                }
                break;
            case 2:
                if (mapFaces.find(geom_id) != mapFaces.end())
                    gEntity = mapFaces[geom_id];
                else
                {
                    cout << "Fatal Error: Geometric Face not found : " << geom_id << endl;
                    exit(0);
                }
                break;
            case 3:
                if (mapCells.find(geom_id) != mapCells.end())
                    gEntity = mapCells[geom_id];
                else
                {
                    cout << "Fatal Error: Geometric Cell not found : " << geom_id << endl;
                    exit(0);
                }
                break;
            default:
                cout << "Error: Invalid geometric dimension " << geom_dim << endl;
                exit(0);
            }

            if (gEntity)
            {
                iRel::Error r_err = mk_core()->irel_pair()->setEntSetRelation(gEntity, entitySets[i]);
                IBERRCHK(r_err, "Trouble set the association between the entity handle and mesh entity set.");	
            }   
        }
    }

    
    
    ///////////////////////////////////////////////////////////////////////////////
    // Create Edge Assocations:
    ///////////////////////////////////////////////////////////////////////////////
    cout << " Building Edge Associations " << endl;

    int numEdges;
    m_err = mk_core()->imesh_instance()->getNumOfType(mesh_root_set, iBase_EDGE, numEdges);
    IBERRCHK(m_err, "Trouble get the number of edges.");	

    std::vector<iBase_EntityHandle> mEdges;

    numAssociations = 0;
    ncount = 0;
    for (int i = 0; i < entitySets.size(); i++)
    {
        mEdges.clear();
        m_err = mk_core()->imesh_instance()->getEntities(entitySets[i], iBase_EDGE, iMesh_ALL_TOPOLOGIES, mEdges);
        IBERRCHK(m_err, "Trouble get the edge entities from the entity set.");

        if (mEdges.size() && (mEdges.size() != numEdges))
        {
            ncount += mEdges.size();

            m_err = mk_core()->imesh_instance()->getEntSetIntData(entitySets[i], mesh_id_tag, geom_id);
            IBERRCHK(m_err, "Trouble set the geom_id for entity set in the edge association.");

            m_err = mk_core()->imesh_instance()->getEntSetIntData(entitySets[i], geom_dim_tag, geom_dim);
            IBERRCHK(m_err, "Trouble set the geom_dim for entity set in the edge association.");
	    
            gEntity = 0;
            switch (geom_dim)
            {
            
            case 1:

                if (mapEdges.find(geom_id) != mapEdges.end())
                {
                    gEntity = mapEdges[geom_id];
                    numAssociations++;
                }
                else
                {
                    cout << "Fatal Error: Geometric Edge not found : " << geom_id << endl;
                    exit(0);
                }
                break;
            case 2:
                if (mapFaces.find(geom_id) != mapFaces.end())
                    gEntity = mapFaces[geom_id];
                else
                {
                    cout << "Fatal Error: Geometric Face not found : " << geom_id << endl;
                    exit(0);
                }
                break;
            case 3:
                if (mapCells.find(geom_id) != mapCells.end())
                    gEntity = mapCells[geom_id];
                else
                {
                    cout << "Fatal Error: Geometric Cell not found : " << geom_id << endl;
                    exit(0);
                }
                break;
            default:
                cout << "Error: Invalid geometric dimension " << geom_dim << endl;
                exit(0);
            }

            if (gEntity)
            {
                iRel::Error r_err = mk_core()->irel_pair()->setEntSetRelation(gEntity, entitySets[i]);
                IBERRCHK(r_err, "Trouble set the association between the entity handle and mesh entity set.");	
            }
        }
    }

    if (numAssociations != mapEdges.size())
        cout << "Warning: There are more edge entitySet than geometric edges " << endl;
    
    //////////////////////////////////////////////////////////////////////////////
    // Face Association
    //////////////////////////////////////////////////////////////////////////////
    cout << " Building Face Associations " << endl;

    std::vector<iBase_EntityHandle> mFaces;

    int numFaces;
    m_err = mk_core()->imesh_instance()->getNumOfType(mesh_root_set, iBase_FACE, numFaces);
    IBERRCHK(m_err, "Trouble get the number of faces.");

    int mesh_faceid;

    ncount = 0;
    numAssociations = 0;
    for (int i = 0; i < entitySets.size(); i++)
    {
        mFaces.clear();
        m_err = mk_core()->imesh_instance()->getEntities(entitySets[i], iBase_FACE, iMesh_ALL_TOPOLOGIES, mFaces);
        IBERRCHK(m_err, "Trouble get the mesh face entities.");

        if (mFaces.size() && (mFaces.size() != numFaces))
        {
            ncount += mFaces.size();
            m_err = mk_core()->imesh_instance()->getEntSetIntData(entitySets[i], mesh_id_tag, geom_id);
            IBERRCHK(m_err, "Trouble set the geom_id for entity set in the face association.");

            m_err = mk_core()->imesh_instance()->getEntSetIntData(entitySets[i], geom_dim_tag, geom_dim);
            IBERRCHK(m_err, "Trouble set the geom_dim for entity set in the face association.");

            gEntity = 0;
            switch (geom_dim)
            {
            case 2:
                if (mapFaces.find(geom_id) != mapFaces.end())
                {
                    gEntity = mapFaces[geom_id];
                    numAssociations++;
                }
                else
                {
                    cout << "Fatal Error: Geometric face not found " << geom_id << endl;
                    exit(0);
                }
                break;
            case 3:
                if (mapCells.find(geom_id) != mapCells.end())
                    gEntity = mapCells[geom_id];
                else
                { IBERRCHK(m_err, "Trouble set the geom_id for entity set in the face association.");
                    cout << "Fatal Error: Geometric face not found " << geom_id << endl;
                    exit(0);
                }
                break;
            }
            if (gEntity)
	    {
                iRel::Error r_err = mk_core()->irel_pair()->setEntSetRelation(gEntity, entitySets[i]);
                IBERRCHK(r_err, "Trouble set the association between the entity handle and mesh entity set.");	
	    }        
	}
    }
    
    if (numAssociations != mapFaces.size())
    {
	cout << "Warning: There are more face entitySet than geometric faces " << endl;
    }
	
    //////////////////////////////////////////////////////////////////////////////
    // Cell Association
    //////////////////////////////////////////////////////////////////////////////

    std::vector<iBase_EntityHandle> mCells;

    int mesh_cellid;

    int numCells;
    m_err = mk_core()->imesh_instance()->getNumOfType(mesh_root_set, iBase_REGION, numCells);
    IBERRCHK(m_err, "Trouble get the number of cells.");

    ncount = 0; IBERRCHK(m_err, "Trouble set the geom_id for entity set in the face association.");
    for (int i = 0; i < entitySets.size(); i++)
    {
        mCells.clear();
        m_err = mk_core()->imesh_instance()->getEntities(entitySets[i], iBase_REGION, iMesh_ALL_TOPOLOGIES, mCells);
        IBERRCHK(m_err, "Trouble get the mesh cell entities.");

        if (mCells.size() && (mCells.size() != numCells))
        {
            ncount += mCells.size();
            m_err = mk_core()->imesh_instance()->getEntSetIntData(entitySets[i], mesh_id_tag, geom_id);
            IBERRCHK(m_err, "Trouble set the geom_id for entity set in the cell association.");

            if (mapCells.find(geom_id) != mapCells.end())
            {
                if (mapCells.find(geom_id) != mapCells.end())
                {
                    gEntity = mapCells[geom_id];
                    iRel::Error r_err = mk_core()->irel_pair()->setEntSetRelation(gEntity, entitySets[i]);
                    IBERRCHK(r_err, "Trouble set the association between the entity handle and mesh entity set.");	
                }
            }
        }
    }
}


}

