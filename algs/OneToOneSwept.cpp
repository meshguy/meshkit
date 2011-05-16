#include "OneToOneSwept.hpp"
#include <iostream>
#include <math.h>
#include <map>


#define Sign(u, v) ( (v)>=0.0 ? Abs(u) : -Abs(u) )
#define Max(u, v) ( (u)>(v)? (u) : (v) )
#define Abs(u) ((u)>0 ? (u) : (-u))
#define Min(u, v) ( (u)>(v)? (v) : (u) )


OneToOneSwept::OneToOneSwept(iGeom_Instance &geometry, iMesh_Instance &Mesh, iRel_Instance &association, iRel_PairHandle &irel)
{
	geom = geometry;
	mesh = Mesh;
	assoc = association;
	rel   = irel;
	
	MeshSetting();
	
	buildAssociation(geom, mesh, assoc, rel);
	
	int err;
	iGeom_getRootSet(geom, &geom_root_set, &err);
    	assert(!err);
	iMesh_getRootSet(mesh, &mesh_root_set, &err);
    	assert(!err);
    	iMesh_createEntSet(mesh, 1, &volumeSet, &err);
    	assert(!err);
	const char *tag = "GLOBAL_ID";
    	int namelen = strlen(tag);
    	iGeom_getTagHandle(geom, tag, &geom_id_tag, &err, namelen);
	assert(!err);
	iMesh_getTagHandle(mesh, tag, &mesh_id_tag, &err, namelen);
	assert(!err);
}

int OneToOneSwept::MeshSetting()
{
    int err;

    SimpleArray<int> adjTable;
    iMesh_getAdjTable(mesh, ARRAY_INOUT(adjTable), &err);
    assert(!err);
    if (adjTable[5] == 0) adjTable[5] = 1;
    if (adjTable[4] == 0) adjTable[4] = 1;
    if (adjTable[10] == 0) adjTable[10] = 1;
    if (adjTable[8] == 0) adjTable[8] = 1;
    if (adjTable[6] == 0) adjTable[6] = 1;
    
    iMesh_setAdjTable(mesh, ARRAY_IN(adjTable), &err);
    if (err)
    {
        char descr[1000];
        int len = 1000;
        iMesh_getDescription(mesh, descr, len);
        cout << descr << endl;
        exit(0);
    }

    return 1;
}

int OneToOneSwept::getList()
{
	//find the corner node list, inner node list and edge node list for the mesh on the source surface
	//confused points: how to find the corner list, inner node list
	int err, entity_index=0;
	SimpleArray<iBase_EntityHandle> Nodes, Edges, Faces;
	iBase_EntitySetHandle SourceSets;
	
	//get the vertex list on the source surface
	iRel_getEntSetRelation(assoc, rel, sourceSurface, 0, &SourceSets, &err);
	Nodes.clear();
	//get inner nodes not boundary nodes
	iMesh_getEntities(mesh, SourceSets, iBase_VERTEX, iMesh_POINT, ARRAY_INOUT(Nodes), &err);
	assert(!err);
	
	iBase_TagHandle taghandle=0;
	iMesh_createTag(mesh, "source", 1, iBase_INTEGER, &taghandle, &err, strlen("source"));
 	assert(!err);
	
	int testnum;
	iMesh_getNumOfType(mesh, SourceSets, iBase_FACE, &testnum, &err);
	assert(!err);	
	
	NodeList.resize(Nodes.size());
	for (int i=0; i < Nodes.size(); i++)
	{
		entity_index++;
		NodeList[entity_index-1].gVertexHandle = Nodes[i];
		NodeList[entity_index-1].index = entity_index-1;
		NodeList[entity_index-1].id = NodeList[entity_index-1].index;
				
		iMesh_getVtxCoord(mesh, Nodes[i], &NodeList[entity_index-1].xyzCoords[0], &NodeList[entity_index-1].xyzCoords[1], &NodeList[entity_index-1].xyzCoords[2], &err);
		assert(!err);
		
		Point3D pts3={NodeList[entity_index-1].xyzCoords[0], NodeList[entity_index-1].xyzCoords[1], NodeList[entity_index-1].xyzCoords[2]};
		Point2D pts2;	
		getUVCoords(sourceSurface, pts3, pts2);
		NodeList[entity_index-1].uvCoords[0] = pts2.pu;
		NodeList[entity_index-1].uvCoords[1] = pts2.pv;

		NodeList[entity_index-1].onBoundary = false;
		NodeList[entity_index-1].onCorner = false;
		
		//set the int data to the entity
		iMesh_setIntData(mesh, Nodes[i], taghandle, NodeList[entity_index-1].id, &err);
		assert(!err);
	}	
	
	//loop over the edges and find the boundary nodes	
	for (unsigned int i=0; i < gsEdgeList.size(); i++)
	{
		iBase_EntitySetHandle tmpSet;
		iRel_getEntSetRelation(assoc, rel, gsEdgeList[i].gEdgeHandle, 0, &tmpSet, &err);
		assert(!err);
		Nodes.clear();
		iMesh_getEntities(mesh, tmpSet, iBase_VERTEX, iMesh_POINT, ARRAY_INOUT(Nodes), &err);
		assert(!err);
		
		NodeList.resize(NodeList.size() + Nodes.size());
		
		for (int j=0; j < Nodes.size(); j++)
		{
			entity_index++;
			NodeList[entity_index-1].gVertexHandle = Nodes[j];
			NodeList[entity_index-1].index = entity_index-1;
			NodeList[entity_index-1].id = entity_index-1;	
			iMesh_getVtxCoord(mesh, Nodes[j], &NodeList[entity_index-1].xyzCoords[0], &NodeList[entity_index-1].xyzCoords[1], &NodeList[entity_index-1].xyzCoords[2], &err);
			assert(!err);
		
			Point3D pts3={NodeList[entity_index-1].xyzCoords[0], NodeList[entity_index-1].xyzCoords[1], NodeList[entity_index-1].xyzCoords[2]};
			Point2D pts2;	
			getUVCoords(sourceSurface, pts3, pts2);
			NodeList[entity_index-1].uvCoords[0] = pts2.pu;
			NodeList[entity_index-1].uvCoords[1] = pts2.pv;

			NodeList[entity_index-1].onBoundary = true;
			NodeList[entity_index-1].onCorner = false;
			
			iMesh_setIntData(mesh, Nodes[j], taghandle, NodeList[entity_index-1].id, &err);
			assert(!err);			
		}	
	}	
	
	//loop over the corners and find the corner nodes
	for (unsigned int i=0; i < gVertexList.size(); i=i+2)
	{
		iBase_EntitySetHandle tmpSet;
		iRel_getEntSetRelation(assoc, rel, gVertexList[i].gVertexHandle, 0, &tmpSet, &err);
		assert(!err);
		Nodes.clear();
		iMesh_getEntities(mesh, tmpSet, iBase_VERTEX, iMesh_POINT, ARRAY_INOUT(Nodes), &err);
		assert(!err);
		
		NodeList.resize(NodeList.size()+Nodes.size());
		
		for (int j=0; j < Nodes.size(); j++)
		{
			entity_index++;
			NodeList[entity_index-1].gVertexHandle = Nodes[j];
			NodeList[entity_index-1].index = entity_index-1;
			NodeList[entity_index-1].id = NodeList[entity_index-1].index;
			iMesh_getVtxCoord(mesh, Nodes[j], &NodeList[entity_index-1].xyzCoords[0], &NodeList[entity_index-1].xyzCoords[1], &NodeList[entity_index-1].xyzCoords[2], &err);
			assert(!err);
		
			Point3D pts3={NodeList[entity_index-1].xyzCoords[0], NodeList[entity_index-1].xyzCoords[1], NodeList[entity_index-1].xyzCoords[2]};
			Point2D pts2;	
			getUVCoords(sourceSurface, pts3, pts2);
			NodeList[entity_index-1].uvCoords[0] = pts2.pu;
			NodeList[entity_index-1].uvCoords[1] = pts2.pv;

			NodeList[entity_index-1].onBoundary = false;
			NodeList[entity_index-1].onCorner = true;
			
			iMesh_setIntData(mesh, Nodes[j], taghandle, NodeList[entity_index-1].id, &err);
			assert(!err);
		}		
	}	
	
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	
	//get the edge list for the source surf mesh
	Edges.clear();
	iMesh_getEntities(mesh, SourceSets, iBase_EDGE, iMesh_LINE_SEGMENT, ARRAY_INOUT(Edges), &err);
	assert(!err);
	
	EdgeList.resize(Edges.size());
	for (int i=0; i < Edges.size(); i++)
	{
		EdgeList[i].gEdgeHandle = Edges[i];
		EdgeList[i].index = i;
		iMesh_getIntData(mesh, Edges[i], mesh_id_tag, &EdgeList[i].EdgeID, &err);
		assert(!err);

		//get the nodes for the edge[i], use the function iMesh_isEntContained
		Nodes.clear();
		iMesh_getEntAdj(mesh, Edges[i], iMesh_POINT, ARRAY_INOUT(Nodes), &err);
		assert(!err);
		
		//loop over the nodes on the edge elements		
		for (int j=0; j < Nodes.size(); j++)
		{
			int tmpIndex=-1;
			iMesh_getIntData(mesh, Nodes[j], taghandle, &tmpIndex, &err);
			assert(!err);
			//find the corresponding nodes on the vertex list
			EdgeList[i].connect[j] = &NodeList[tmpIndex];
		}

		//determine whether the edge is a boundary edge element or an inner edge element
		if (isEdgeBoundary(Edges[i]))
			EdgeList[i].onBoundary = true;
		else
			EdgeList[i].onBoundary = false;
	}
	
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	Faces.clear();
	iMesh_getEntities(mesh, SourceSets, iBase_FACE, iMesh_ALL_TOPOLOGIES, ARRAY_INOUT(Faces), &err);
	assert(!err);
	
	FaceList.resize(Faces.size());
	for (int i=0; i < Faces.size(); i++)
	{
		FaceList[i].gFaceHandle = Faces[i];
		FaceList[i].index = i;
		iMesh_getIntData(mesh, Faces[i], mesh_id_tag, &FaceList[i].FaceID, &err);
		assert(!err);
		//get the nodes on the face elements
		Nodes.clear();
		iMesh_getEntAdj(mesh, Faces[i], iBase_VERTEX, ARRAY_INOUT(Nodes), &err);
		assert(!err);

		FaceList[i].connect.resize(Nodes.size());
		for (int j=0; j < Nodes.size(); j++)
		{			
			int tmpIndex=-1;
			iMesh_getIntData(mesh, Nodes[j], taghandle, &tmpIndex, &err);
			assert(!err);
			
			FaceList[i].connect[j] = &NodeList[tmpIndex];
		}
	}
	
	//initialize the mesh size on the target surface
	TVertexList.resize(NodeList.size());
	TEdgeList.resize(EdgeList.size());
	TFaceList.resize(FaceList.size());
		
	return 1;
}

int OneToOneSwept::isVertexCorner(iBase_EntityHandle gNodeHandle, iBase_EntityHandle gFaceHandle)
{
	return 1;
}

int OneToOneSwept::isEdgeBoundary(iBase_EntityHandle gEdgeHandle)
{
	int err;
	SimpleArray<iBase_EntityHandle> Faces;	
	Faces.clear();
	iMesh_getEntAdj(mesh, gEdgeHandle, iMesh_POLYGON, ARRAY_INOUT(Faces), &err);
	assert(!err);
	if (Faces.size()==1)
		return 1;
	else
		return 0;
}

void OneToOneSwept::Execute()
{	
	cout << "Meshing...." << endl;

	//TestGeom();
	
	getList();
	
	FindCorners();
	TargetSurfProjection();	
	InnerLayerMeshing();
}

int OneToOneSwept::FindCorners()
{
	int err;
	Point3D pts;
	Point2D uv;
	iBase_TagHandle taghandle;
	iMesh_getTagHandle(mesh, "source", &taghandle, &err, strlen("source"));
	assert(!err);
	
	//Collect all the feature nodes on the target surface
	for (unsigned int i=0; i < gVertexList.size(); i=i+2)
	{
		//find the corresponding the node on the source surface
		int ID;
		int tindex = cornerPairs[i];//get the corresponding corner on the target surface
		
		iBase_EntitySetHandle tmpSet;
		iRel_getEntSetRelation(assoc, rel, gVertexList[i].gVertexHandle, 0, &tmpSet, &err);
		assert(!err);
		SimpleArray<iBase_EntityHandle> Nodes;
		iMesh_getEntities(mesh, tmpSet, iBase_VERTEX, iMesh_POINT, ARRAY_INOUT(Nodes), &err);
		assert(!err);
		assert(Nodes.size()==1);
		
		iMesh_getIntData(mesh, Nodes[0], taghandle, &ID, &err);
		assert(!err);
		
		TVertexList[ID].index = ID;

		TVertexList[ID].xyzCoords[0] = gVertexList[tindex].xyzCoords[0];
		TVertexList[ID].xyzCoords[1] = gVertexList[tindex].xyzCoords[1];
		TVertexList[ID].xyzCoords[2] = gVertexList[tindex].xyzCoords[2];
		
		TVertexList[ID].onCorner = true;
		TVertexList[ID].onBoundary = false;
		
		//determine whether there is a mesh vertex at this position
		iRel_getEntSetRelation(assoc, rel, gVertexList[tindex].gVertexHandle, 0, &tmpSet, &err);
		if (err)
		{//there is no existing mesh vertex at this location
			iMesh_createVtx(mesh, gVertexList[tindex].xyzCoords[0], gVertexList[tindex].xyzCoords[1], gVertexList[tindex].xyzCoords[2], &TVertexList[ID].gVertexHandle, &err);
			assert(!err);
			iBase_EntitySetHandle entityset;
			iMesh_createEntSet(mesh, 1, &entityset, &err);
			assert(!err);
			iMesh_addEntToSet(mesh, TVertexList[ID].gVertexHandle, entityset, &err);
			assert(!err);
		
			iRel_setEntSetRelation(assoc, rel, gVertexList[tindex].gVertexHandle, entityset, &err);
			assert(!err);
			pts.px = gVertexList[tindex].xyzCoords[0];
			pts.py = gVertexList[tindex].xyzCoords[1];
			pts.pz = gVertexList[tindex].xyzCoords[2];
			getUVCoords(targetSurface, pts, uv);	
		}
		else
		{//there is an existing vertex at this location
			iRel_getEntSetRelation(assoc, rel, gVertexList[tindex].gVertexHandle, 0, &tmpSet, &err);
			assert(!err);
			
			SimpleArray<iBase_EntityHandle> tmpNodes;
			iMesh_getEntities(mesh, tmpSet, iBase_VERTEX, iMesh_POINT, ARRAY_INOUT(tmpNodes), &err);
			assert(!err);
			assert(Nodes.size()==1);
			
			pts.px = gVertexList[tindex].xyzCoords[0];
			pts.py = gVertexList[tindex].xyzCoords[1];
			pts.pz = gVertexList[tindex].xyzCoords[2];
			getUVCoords(targetSurface, pts, uv);
			TVertexList[ID].gVertexHandle = tmpNodes[0];	
		}
		TVertexList[ID].uvCoords[0] = uv.pu;
		TVertexList[ID].uvCoords[1] = uv.pv;
		
		
	}

	return 1;
}

OneToOneSwept::~OneToOneSwept()
{
	cout << "*********************************************************************************" << endl;
	int num_vertex, num_edge, num_face, num_cell, err;
	iMesh_getNumOfType(mesh, mesh_root_set, iBase_VERTEX, &num_vertex, &err);
	assert(!err);
	iMesh_getNumOfType(mesh, mesh_root_set, iBase_EDGE, &num_edge, &err);
	assert(!err);
	iMesh_getNumOfType(mesh, mesh_root_set, iBase_FACE, &num_face, &err);
	assert(!err);
	iMesh_getNumOfType(mesh, mesh_root_set, iBase_REGION, &num_cell, &err);
	assert(!err);
	cout << "Number of Vertices: " << num_vertex << endl;
	cout << "Number of Edges   : " << num_edge << endl;
	cout << "Number of faces   : " << num_face << endl;
	cout << "Number of cells   : " << num_cell << endl;
	cout << "*********************************************************************************" << endl;
}

int OneToOneSwept::TargetSurfProjection()
{
	cout << "Target surface meshing..." << endl;
	int err;
	iBase_EntitySetHandle targetSet;
	iRel_getEntSetRelation(assoc, rel, targetSurface, 0, &targetSet, &err);
	if (!err)
	{//there exists an mesh on the target surfaces
		return 1;
	}
	iBase_TagHandle taghandle;
	iMesh_getTagHandle(mesh, "source", &taghandle, &err, strlen("source"));
	assert(!err);
	
	int index=0, id, index_t;
	double su, tu;
	SimpleArray<iBase_EntityHandle> gsEdge, gtEdge, edgeNodes, geomEdgeEnds, meshEdgeEnds;
	vector<iBase_EntitySetHandle> edgeEndsMeshSets;
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
	vector<iBase_EntitySetHandle> edgeMeshSets(gsEdgeList.size());
	//loop over the various edges
	
	for (unsigned int i=0; i < gsEdgeList.size(); i++)
	{
		//get the mesh entityset for edge[i]
		
		iRel_getEntSetRelation(assoc, rel, gsEdgeList[i].gEdgeHandle, 0, &edgeMeshSets[i], &err);
		assert(!err);
		
		//get the edge nodes for edge[i] mesh
		edgeNodes.clear();
		iMesh_getEntities(mesh, edgeMeshSets[i], iBase_VERTEX, iMesh_POINT, ARRAY_INOUT(edgeNodes), &err);
		assert(!err);
		//find the corresponding relationship for both ends on the source suface and target surface
		//loop over the linking sides
		int index_a, index_b; //record the vertex numbering on the target surface
		index_a = cornerPairs[gsEdgeList[i].connect[0]->index];
		index_b = cornerPairs[gsEdgeList[i].connect[1]->index];
		
		//create the corresponding edge relationship between the source surface and target surface
		index_t = edgePairs[i];
		
		//detect whether there exists an mesh for the edge on the target surface
		iBase_EntitySetHandle entityset;
		iRel_getEntSetRelation(assoc, rel, gtEdgeList[index_t].gEdgeHandle, 0, &entityset, &err);
		if (!err)
		{//there exists an mesh on the boundary of target surface
			SimpleArray<iBase_EntityHandle> tnodes;
			iMesh_getEntities(mesh, entityset, iBase_VERTEX, iMesh_POINT, ARRAY_INOUT(tnodes), &err);
			assert(!err);
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
				iMesh_getVtxCoord(mesh, edgeNodes[0], &xyz1[0], &xyz1[1], &xyz1[2], &err);
				assert(!err);
				iMesh_getVtxCoord(mesh, edgeNodes[1], &xyz2[0], &xyz2[1], &xyz2[2], &err);
				assert(!err);
				
				dist1 = sqrt( pow(gsEdgeList[i].connect[0]->xyzCoords[0]-xyz1[0],2) + pow(gsEdgeList[i].connect[0]->xyzCoords[1]-xyz1[1],2) + pow(gsEdgeList[i].connect[0]->xyzCoords[2]-xyz1[2],2));
				dist2 = sqrt( pow(gsEdgeList[i].connect[0]->xyzCoords[0]-xyz2[0],2) + pow(gsEdgeList[i].connect[0]->xyzCoords[1]-xyz2[1],2) + pow(gsEdgeList[i].connect[0]->xyzCoords[2]-xyz2[2],2));
				if (dist1 >= dist2)
				{
					isReverse1 = true;
				}
				iMesh_getVtxCoord(mesh, tnodes[0], &xyz1[0], &xyz1[1], &xyz1[2], &err);
				assert(!err);
				iMesh_getVtxCoord(mesh, tnodes[1], &xyz2[0], &xyz2[1], &xyz2[2], &err);
				assert(!err);
				
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
				iMesh_getIntData(mesh, edgeNodes[j], taghandle, &tmpIndex, &err);
				assert(!err);
				
				index++;
				
				sPtsUV.resize(index);
				tPtsUV.resize(index);
				
				pts.px = NodeList[tmpIndex].xyzCoords[0];
				pts.py = NodeList[tmpIndex].xyzCoords[1];
				pts.pz = NodeList[tmpIndex].xyzCoords[2];
				
				getUVCoords(sourceSurface, pts, sPtsUV[index-1]);
				
				if (isReverse1!=isReverse2)
				{
					iMesh_getVtxCoord(mesh, tnodes[edgeNodes.size()-j-1], &TVertexList[tmpIndex].xyzCoords[0], &TVertexList[tmpIndex].xyzCoords[1], &TVertexList[tmpIndex].xyzCoords[2], &err);
					assert(!err);
					
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
					iMesh_getVtxCoord(mesh, tnodes[j], &TVertexList[tmpIndex].xyzCoords[0], &TVertexList[tmpIndex].xyzCoords[1], &TVertexList[tmpIndex].xyzCoords[2], &err);
					assert(!err);
					
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
		iGeom_getEntXYZtoU(geom, gsEdgeList[i].gEdgeHandle, gsEdgeList[i].connect[0]->xyzCoords[0], gsEdgeList[i].connect[0]->xyzCoords[1], gsEdgeList[i].connect[0]->xyzCoords[2], &sLeft, &err);
		assert(!err);
		iGeom_getEntXYZtoU(geom, gsEdgeList[i].gEdgeHandle, gsEdgeList[i].connect[1]->xyzCoords[0], gsEdgeList[i].connect[1]->xyzCoords[1], gsEdgeList[i].connect[1]->xyzCoords[2], &sRight, &err);
		assert(!err);
		iGeom_getEntXYZtoU(geom, gtEdgeList[index_t].gEdgeHandle, gVertexList[index_a].xyzCoords[0], gVertexList[index_a].xyzCoords[1], gVertexList[index_a].xyzCoords[2], &tLeft, &err);
		assert(!err);
		iGeom_getEntXYZtoU(geom, gtEdgeList[index_t].gEdgeHandle, gVertexList[index_b].xyzCoords[0], gVertexList[index_b].xyzCoords[1], gVertexList[index_b].xyzCoords[2], &tRight, &err);
		assert(!err);

		vector<iBase_EntityHandle> newNodehandle(edgeNodes.size());
		SimpleArray<iBase_EntityHandle> testNodeHandle;
		//create the boundary node on the edges. This doesn't include the corners.
		for (int j=0; j < edgeNodes.size(); j++)
		{
			//get the cartesian coordinates for the edge nodes
			iMesh_getVtxCoord(mesh, edgeNodes[j], &pts.px, &pts.py, &pts.pz, &err);
			assert(!err);

			//get the parametric coordinates for the vertex on the source surface
			index++;
			sPtsUV.resize(index);
			getUVCoords(sourceSurface, pts, sPtsUV[index-1]);

			//get vertex id for target surface vertex
			iMesh_getIntData(mesh, edgeNodes[j], taghandle, &id, &err);
			assert(!err);
			
			TVertexList[id].index = id;		
						
			//transform the cartesian coordinates into the parametric coordinates
			iGeom_getEntXYZtoU(geom, gsEdgeList[i].gEdgeHandle, pts.px, pts.py, pts.pz, &su, &err);
			assert(!err);

			//calculate the parametric coordinates on the target surface
			tu = tLeft + (tRight - tLeft) * (su - sLeft) / (sRight - sLeft);

			//transform the parametric coordinates into the cartesian coordinates on the target surface
			iGeom_getEntUtoXYZ(geom, gtEdgeList[index_t].gEdgeHandle, tu, &pts.px, &pts.py, &pts.pz, &err);
			assert(!err);

			//create the vertex entity on the target boundary edge
			iMesh_createVtx(mesh, pts.px, pts.py, pts.pz, &newNodehandle[j], &err);
			assert(!err);

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
		iRel_getEntSetRelation(assoc, rel, gtEdgeList[index_t].gEdgeHandle, 0, &entityset, &err);
		if (err) //there is no entityset associated with gtEdgeList[index_t].gEdgeHandle
		{
			iMesh_createEntSet(mesh, 1, &entityset, &err);
			assert(!err);
		}
		
		//build the association
		iMesh_addEntArrToSet(mesh, &newNodehandle[0], edgeNodes.size(), entityset, &err);
		assert(!err);
		
		//create the line segments on the boundary edge of target surface
		SimpleArray<iBase_EntityHandle> sedgeHandle;
		vector<iBase_EntityHandle> tedgeHandle(0);
		int edgeIndex=0;
		iMesh_getEntities(mesh, edgeMeshSets[i], iBase_EDGE, iMesh_LINE_SEGMENT, ARRAY_INOUT(sedgeHandle), &err);
		assert(!err);
		for (int j=0; j< sedgeHandle.size(); j++)
		{
			int status;
			int nodeindex1, nodeindex2;
			vector<iBase_EntityHandle> connect(2);
			SimpleArray<iBase_EntityHandle> tmpNodes;
			iMesh_getEntAdj(mesh, sedgeHandle[j], iBase_VERTEX, ARRAY_INOUT(tmpNodes), &err);
			assert(!err);
			assert(tmpNodes.size()==2);
			iMesh_getIntData(mesh, tmpNodes[0], taghandle, &nodeindex1, &err);
			assert(!err);
			iMesh_getIntData(mesh, tmpNodes[1], taghandle, &nodeindex2, &err);
			assert(!err);
			
			edgeIndex++;
			tedgeHandle.resize(edgeIndex);
			connect[0] = TVertexList[nodeindex1].gVertexHandle;
			connect[1] = TVertexList[nodeindex2].gVertexHandle;
			
			iMesh_createEnt(mesh, iMesh_LINE_SEGMENT, &connect[0], 2, &tedgeHandle[edgeIndex-1], &status, &err);
			assert(!err);
			
		}
		iMesh_addEntArrToSet(mesh, &tedgeHandle[0], tedgeHandle.size(), entityset, &err);
		assert(!err);
		
		iRel_getEntSetRelation(assoc, rel, gtEdgeList[index_t].gEdgeHandle, 0, &entityset, &err);
		if (err) //there is no entityset associated with gtEdgeList[index_t].gEdgeHandle
		{
			iRel_setEntSetRelation(assoc, rel, gtEdgeList[index_t].gEdgeHandle, entityset, &err);
			assert(!err);
		}								
	}
	
	//Until now, all the nodes have been created on the boundary edge.
	//get the corner coordinates
	assert(NodeList.size()==TVertexList.size());
	for (unsigned int i=0; i < NodeList.size(); i++)
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
	iRel_getEntSetRelation(assoc, rel, targetSurface, 0, &entityset, &err);
	if (err) //there is no entityset associated with targetSurface
	{
		iMesh_createEntSet(mesh, 1, &entityset, &err);
		assert(!err);
	}

	//create the inner nodes on the target surface
	for (unsigned int i=0; i < NodeList.size(); i++)
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
			iMesh_createVtx(mesh, pts.px, pts.py, pts.pz, &newNodehandle[newIndex-1], &err);
			assert(!err);

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
	iMesh_addEntArrToSet(mesh, &newNodehandle[0], newIndex, entityset, &err);
	assert(!err);

	//until now, all the nodes have been generated on the target surface

	//create the edge elements on the target surface
	int status;
	index =0;
	vector<iBase_EntityHandle> gEdgeHandle(0);
	for (unsigned int i=0; i < EdgeList.size(); i++)
	{
		if (!EdgeList[i].onBoundary)
		{
			//create the inner edge elements
			vector<iBase_EntityHandle> connect(2);
			connect[0] = TVertexList[EdgeList[i].connect[0]->index].gVertexHandle;
			connect[1] = TVertexList[EdgeList[i].connect[1]->index].gVertexHandle;
			
			index++;			
			gEdgeHandle.resize(index);

			iMesh_createEnt(mesh, iMesh_LINE_SEGMENT, &connect[0], 2, &gEdgeHandle[index-1], &status, &err);
			assert(!err);
			
			//add the new generated inner edge elements on the target surface into the list
			TEdgeList[i].index = i;
			TEdgeList[i].connect[0] = &TVertexList[EdgeList[i].connect[0]->index];
			TEdgeList[i].connect[1] = &TVertexList[EdgeList[i].connect[1]->index];
			TEdgeList[i].gEdgeHandle = gEdgeHandle[index-1];
				
		}		
	}
	//add the inner edge elements to the entityset
	iMesh_addEntArrToSet(mesh, &gEdgeHandle[0], index, entityset, &err);
	assert(!err);



	//determine the numbering order for quadrilateral nodes
	int sense_out, sense_out1, sense_out2, sense_out3, sense_out4;
	iGeom_getEgFcSense(geom, gsEdgeList[0].gEdgeHandle, sourceSurface, &sense_out1, &err);
	assert(!err);
	
	iGeom_getEgFcSense(geom, gtEdgeList[edgePairs[gsEdgeList[0].index]].gEdgeHandle, targetSurface, &sense_out2, &err);
	assert(!err);
	
	iGeom_getEgVtxSense( geom, gsEdgeList[0].gEdgeHandle, gsEdgeList[0].connect[0]->gVertexHandle, gsEdgeList[0].connect[1]->gVertexHandle, &sense_out3, &err );
	assert(!err);
	
	iGeom_getEgVtxSense( geom, gtEdgeList[edgePairs[gsEdgeList[0].index]].gEdgeHandle, gVertexList[cornerPairs[gsEdgeList[0].connect[0]->index]].gVertexHandle, gVertexList[cornerPairs[gsEdgeList[0].connect[1]->index]].gVertexHandle, &sense_out4, &err );
	assert(!err);
	sense_out = sense_out1*sense_out2*sense_out3*sense_out4;


	//create the quadrilateral elements on the target surface
	vector<iBase_EntityHandle> mFaceHandle(FaceList.size());
	for (unsigned int i=0; i < FaceList.size(); i++)
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
		iMesh_createEnt(mesh, iMesh_QUADRILATERAL, &connect[0], FaceList[i].getNumNodes(), &mFaceHandle[i], &status, &err);
		assert(!err);
		
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
	iMesh_addEntArrToSet(mesh, &mFaceHandle[0], FaceList.size(), entityset, &err);
	assert(!err);
	
	//build the association
	iRel_getEntSetRelation(assoc, rel, targetSurface, 0, &entityset, &err);
	if (err) //there is no entityset associated with region[0]
	{
		iRel_setEntSetRelation(assoc, rel, targetSurface, entityset, &err);
		assert(!err);	
	}
	return 1;
}

int OneToOneSwept::getXYZCoords(iBase_EntityHandle gFaceHandle, Point3D &pts3, double uv[2])
{
	int err;
	double umin, umax, vmin, vmax;
	iGeom_getEntUVRange(geom, gFaceHandle, &umin, &vmin, &umax, &vmax, &err);
	assert(!err);
	if ((uv[0]<umin)||(uv[0]>umax))
		cout << "Warning: U exceeds the range" << endl;
	if ((uv[1]<vmin)||(uv[1]>vmax))
		cout << "Warning: V exceeds the range" << endl;
	
	iGeom_getEntUVtoXYZ(geom, gFaceHandle, uv[0], uv[1], &pts3.px, &pts3.py, &pts3.pz, &err);
    	assert(!err);

	return 1;
}

int OneToOneSwept::getUVCoords(iBase_EntityHandle gFaceHandle, Point3D pts3, Point2D &pts2)
{
	int err;
	double xmin, ymin, zmin, xmax, ymax, zmax;

	iGeom_getEntBoundBox(geom, gFaceHandle, &xmin, &ymin, &zmin, &xmax, &ymax, &zmax, &err);
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
	iGeom_getEntXYZtoUV(geom, gFaceHandle, pts3.px, pts3.py, pts3.pz, &pts2.pu, &pts2.pv, &err);
    	assert(!err);

	return 1;	
}

void OneToOneSwept::SurfaceSpecifying()
{
	int err, choice, index=0;
	
	SimpleArray<iBase_EntityHandle> gFaces, gRegions, gEdges, gsEdges, gtEdges, gsNodes, gtNodes, gNodes, gVolume;
	iBase_EntitySetHandle notLinkSides;
	
	iGeom_getEntities(geom, geom_root_set, iBase_REGION, ARRAY_INOUT(gVolume), &err);
	assert(!err);
	//make sure this is one-to-one sweeeping
	assert(gVolume.size()==1);
	volEntity = gVolume[0];
	iRel_setEntSetRelation(assoc, rel, gVolume[0], volumeSet, &err);
	assert(!err);
	
	gFaces.clear();
	iGeom_getEntities(geom, geom_root_set, iBase_FACE, ARRAY_INOUT(gFaces), &err);
	assert(!err);
	gRegions.clear();
	iGeom_getEntities(geom, geom_root_set, iBase_REGION, ARRAY_INOUT(gRegions), &err);
	assert(!err);
	
	for (int i=0; i < gFaces.size(); i++)
	{	
		gNodes.clear();
		iGeom_getEntAdj(geom, gFaces[i], iBase_VERTEX, ARRAY_INOUT(gNodes), &err);
		assert(!err);
		SimpleArray<double> coords;
		cout << "-------------- Face " <<  i << " --------------" << endl;
		coords.clear();
		iGeom_getVtxArrCoords(geom, &gNodes[0], gNodes.size(), 1, ARRAY_INOUT(coords), &err);
		assert(!err);
		for (int j=0; j < gNodes.size(); j++)
		{
			cout << "Node [" << j << "]'s coordinates are x:" << coords[3*j] << ", y:" << coords[3*j+1] << ", z:" << coords[3*j+2] << endl;
		}
	}

	cout << "Please specify the source surface" << endl;
	cin >> choice;
	sourceSurface = gFaces[choice];
	sourceRegion = gRegions[choice];
	cout << "The source surface has been set up" << endl;

	cout << "Please specify the target surface" << endl;
	cin >> choice;
	targetSurface = gFaces[choice];
	targetRegion = gRegions[choice];
	cout << "The target surface has been set up" << endl;

	//specify the linking surface
	cout << "Please specify the linking surfaces" << endl;
	index=0;
	for (int i=0; i < (gFaces.size()-2); i++)
	{
		cin >> choice;
		index++;
		gLinkFaceList.resize(index);
		gLinkFaceList[index-1].gFaceHandle = gFaces[choice];
		gLinkFaceList[index-1].index = index-1;
		iGeom_getIntData(geom, gFaces[choice], geom_id_tag, &gLinkFaceList[index-1].FaceID, &err);
		assert(!err);
	}
	cout << "Linking surfaces have been set up" << endl;

	//create geometrical vertex list
	gsNodes.clear();
	iGeom_getEntAdj(geom, sourceSurface, iBase_VERTEX, ARRAY_INOUT(gsNodes), &err);
	assert(!err);
	gtNodes.clear();
	iGeom_getEntAdj(geom, targetSurface, iBase_VERTEX, ARRAY_INOUT(gtNodes), &err);
	assert(!err);
	assert(gsNodes.size()==gtNodes.size());
	index=0;
	for (int i=0; i < gsNodes.size(); i++)
	{
		Point3D pts;
		Point2D uv;
		iGeom_getVtxCoord(geom, gsNodes[i], &pts.px, &pts.py, &pts.pz, &err);
		assert(!err);
		iGeom_getEntXYZtoUV(geom, sourceSurface, pts.px, pts.py, pts.pz, &uv.pu, &uv.pv, &err);
		assert(!err);
		
		index++;
		gVertexList.resize(index);
		gVertexList[index-1].gVertexHandle = gsNodes[i];
		gVertexList[index-1].index = 2*i;
		iGeom_getIntData(geom, gsNodes[i], geom_id_tag, &gVertexList[index-1].id, &err);
		assert(!err);
		gVertexList[index-1].xyzCoords[0] = pts.px;
		gVertexList[index-1].xyzCoords[1] = pts.py;
		gVertexList[index-1].xyzCoords[2] = pts.pz;
		gVertexList[index-1].uvCoords[0] = uv.pu;
		gVertexList[index-1].uvCoords[1] = uv.pv;
		
		iGeom_getVtxCoord(geom, gtNodes[i], &pts.px, &pts.py, &pts.pz, &err);
		assert(!err);
		iGeom_getEntXYZtoUV(geom, targetSurface, pts.px, pts.py, pts.pz, &uv.pu, &uv.pv, &err);
		assert(!err);

		index++;
		gVertexList.resize(index);
		gVertexList[index-1].gVertexHandle = gtNodes[i];
		gVertexList[index-1].index = 2*i+1;
		iGeom_getIntData(geom, gtNodes[i], geom_id_tag, &gVertexList[index-1].id, &err);
		assert(!err);
		gVertexList[index-1].xyzCoords[0] = pts.px;
		gVertexList[index-1].xyzCoords[1] = pts.py;
		gVertexList[index-1].xyzCoords[2] = pts.pz;
		gVertexList[index-1].uvCoords[0] = uv.pu;
		gVertexList[index-1].uvCoords[1] = uv.pv;
	}

	//create geometrical edge list
	gsEdges.clear();
	iGeom_getEntAdj(geom, sourceSurface, iBase_EDGE, ARRAY_INOUT(gsEdges), &err);
	assert(!err);
	gtEdges.clear();
	iGeom_getEntAdj(geom, targetSurface, iBase_EDGE, ARRAY_INOUT(gtEdges), &err);
	assert(!err);

	assert(gsEdges.size()==gtEdges.size());
	index=0;

	for (int i=0; i < gsEdges.size(); i++)
	{
		Point3D pts[2];	
		index++;
		gsEdgeList.resize(index);
		gsEdgeList[index-1].gEdgeHandle = gsEdges[i];
		gsEdgeList[index-1].index = i;
		iGeom_getIntData(geom, gsEdges[i], geom_id_tag, &gsEdgeList[index-1].EdgeID, &err);
		assert(!err);
		//create the linking between nodes on the edge and vertex list
		gsNodes.clear();		
		iGeom_getEntAdj(geom, gsEdges[i], iBase_VERTEX, ARRAY_INOUT(gsNodes), &err);
		assert(!err);
		assert(gsNodes.size()==2);

		iGeom_getVtxCoord(geom, gsNodes[0], &pts[0].px, &pts[0].py, &pts[0].pz, &err);
		assert(!err);
		iGeom_getVtxCoord(geom, gsNodes[1], &pts[1].px, &pts[1].py, &pts[1].pz, &err);
		assert(!err);

		for (unsigned int j=0; j < gVertexList.size(); j=j+2)
		{
			if (fabs(sqrt(pow(pts[0].px-gVertexList[j].xyzCoords[0],2)+pow(pts[0].py-gVertexList[j].xyzCoords[1],2)+pow(pts[0].pz-gVertexList[j].xyzCoords[2],2)))<1e-10)
			{
				gsEdgeList[index-1].connect[0] = &gVertexList[j];
				continue;
			}
			if (fabs(sqrt(pow(pts[1].px-gVertexList[j].xyzCoords[0],2)+pow(pts[1].py-gVertexList[j].xyzCoords[1],2)+pow(pts[1].pz-gVertexList[j].xyzCoords[2],2)))<1e-10)
			{
				gsEdgeList[index-1].connect[1] = &gVertexList[j];
				continue;
			}
		}		
		
		gtEdgeList.resize(index);
		gtEdgeList[index-1].gEdgeHandle = gtEdges[i];
		gtEdgeList[index-1].index = i;
		iGeom_getIntData(geom,  gtEdges[i], geom_id_tag, &gtEdgeList[index-1].EdgeID, &err);
		assert(!err);
		//create the linking between nodes on the edge and vertex list
		gtNodes.clear();		
		iGeom_getEntAdj(geom, gtEdges[i], iBase_VERTEX, ARRAY_INOUT(gtNodes), &err);
		assert(!err);
		assert(gtNodes.size()==2);

		iGeom_getVtxCoord(geom, gtNodes[0], &pts[0].px, &pts[0].py, &pts[0].pz, &err);
		assert(!err);
		iGeom_getVtxCoord(geom, gtNodes[1], &pts[1].px, &pts[1].py, &pts[1].pz, &err);
		assert(!err);

		for (unsigned int j=1; j < gVertexList.size(); j=j+2)
		{
			if (fabs(sqrt(pow(pts[0].px-gVertexList[j].xyzCoords[0],2)+pow(pts[0].py-gVertexList[j].xyzCoords[1],2)+pow(pts[0].pz-gVertexList[j].xyzCoords[2],2)))<1e-10)
			{
				gtEdgeList[index-1].connect[0] = &gVertexList[j];
				continue;
			}
			if (fabs(sqrt(pow(pts[1].px-gVertexList[j].xyzCoords[0],2)+pow(pts[1].py-gVertexList[j].xyzCoords[1],2)+pow(pts[1].pz-gVertexList[j].xyzCoords[2],2)))<1e-10)
			{
				gtEdgeList[index-1].connect[1] = &gVertexList[j];
				continue;
			}
		}
	}

	//this part needs to be improved
	cout << "Please specify the link sides" << endl;

	iGeom_createEntSet(geom, 1, &notLinkSides, &err);
	assert(!err);

	//add the edges on the source surface and target surface to the entitysets
	iGeom_addEntArrToSet(geom, &gsEdges[0], gsEdges.size(), notLinkSides, &err);
	assert(!err);
	iGeom_addEntArrToSet(geom, &gtEdges[0], gtEdges.size(), notLinkSides, &err);
	assert(!err);
	
	gEdges.clear();
	iGeom_getEntities(geom, geom_root_set, iBase_EDGE, ARRAY_INOUT(gEdges), &err);	
	assert(!err);
	int is_contained;
	index=0;
	gLinkSides.resize(index);
	for (int i=0; i < gEdges.size(); i++)
	{
		iGeom_isEntContained(geom, notLinkSides, gEdges[i], &is_contained, &err);
		assert(!err);
		if (!is_contained)
		{
			Point3D pts[2];
			index++;
			gLinkSides.resize(index);
			gLinkSides[index-1].gEdgeHandle = gEdges[i];
			gLinkSides[index-1].index = index-1;
			iGeom_getIntData(geom, gEdges[i], geom_id_tag, &gLinkSides[index-1].EdgeID, &err);
			assert(!err);
			//get ends for linking sides
			gNodes.clear();		
			iGeom_getEntAdj(geom, gEdges[i], iBase_VERTEX, ARRAY_INOUT(gNodes), &err);
			assert(!err);
			assert(gNodes.size()==2);

			iGeom_getVtxCoord(geom, gNodes[0], &pts[0].px, &pts[0].py, &pts[0].pz, &err);
			assert(!err);
			iGeom_getVtxCoord(geom, gNodes[1], &pts[1].px, &pts[1].py, &pts[1].pz, &err);
			assert(!err);		
	
			for (unsigned int j=0; j < gVertexList.size(); j++)
			{
				if (fabs(sqrt(pow(pts[0].px-gVertexList[j].xyzCoords[0],2)+pow(pts[0].py-gVertexList[j].xyzCoords[1],2)+pow(pts[0].pz-gVertexList[j].xyzCoords[2],2)))<1e-10)
				{
					gLinkSides[index-1].connect[0] = &gVertexList[j];
					continue;
				}
				if (fabs(sqrt(pow(pts[1].px-gVertexList[j].xyzCoords[0],2)+pow(pts[1].py-gVertexList[j].xyzCoords[1],2)+pow(pts[1].pz-gVertexList[j].xyzCoords[2],2)))<1e-10)
				{
					gLinkSides[index-1].connect[1] = &gVertexList[j];
					continue;
				}
			}	
		}

	}

	iGeom_destroyEntSet(geom, notLinkSides, &err);
	assert(!err);

	cout << "The linking sides have been set up" << endl;


	//specify the element size
	cout << "Please specify the element size" << endl;
	cin >> numLayers;

	cout << "The element size has been specified" << endl;

	//create the corner relationship between source surface and target surface
	for (unsigned int i=0; i < gVertexList.size(); i=i+2)
	{
		int index_a;
		for (unsigned int j=0; j < gLinkSides.size(); j++)
		{
			if (gVertexList[i].id == gLinkSides[j].connect[0]->id)
			{
				index_a = gLinkSides[j].connect[1]->index;
				continue;
			}
			if (gVertexList[i].id == gLinkSides[j].connect[1]->id)
			{
				index_a = gLinkSides[j].connect[0]->index;
				continue;
			}						
		}
		cornerPairs[i]=index_a;
	}	
	
	//create the edge relationship between the source surface and target surface
	for (unsigned int i=0; i < gsEdgeList.size(); i++)
	{
		int index_a = cornerPairs[gsEdgeList[i].connect[0]->index];
		int index_b = cornerPairs[gsEdgeList[i].connect[1]->index];	 	
		for (unsigned int j=0; j < gtEdgeList.size(); j++)
		{
			if (((gVertexList[index_a].id == gtEdgeList[j].connect[0]->id)&&(gVertexList[index_b].id == gtEdgeList[j].connect[1]->id))||((gVertexList[index_a].id == gtEdgeList[j].connect[1]->id)&&(gVertexList[index_b].id == gtEdgeList[j].connect[0]->id)))
			{	
				edgePairs[i]=j;
				break;
			}
		}
	}

	//add an array of vertex to the linking surf list
	for (unsigned int i=0; i < gLinkFaceList.size(); i++)
	{
		gNodes.clear();
		iGeom_getEntAdj(geom, gFaces[i], iBase_VERTEX, ARRAY_INOUT(gNodes), &err);
		assert(!err);
		gLinkFaceList[i].connect.resize(gNodes.size());
		for (int j=0; j < gNodes.size(); j++)
		{
			Point3D pts;
			iGeom_getVtxCoord(geom, gNodes[j], &pts.px, &pts.py, &pts.pz, &err);
			assert(!err);		
			for (unsigned int k=0; k < gVertexList.size(); k++)
			{
				if (is_SamePoint(pts.px, pts.py, pts.pz, gVertexList[k].xyzCoords[0], gVertexList[k].xyzCoords[1], gVertexList[k].xyzCoords[2]))
				{
					gLinkFaceList[i].connect[j] = &gVertexList[k];
				}
				else
					continue;
			}

		}
	}
}

int OneToOneSwept::is_SamePoint(double ptx1, double pty1, double ptz1, double ptx2, double pty2, double ptz2)
{
	double distance;
	distance = sqrt(pow((ptx1-ptx2),2) + pow((pty1-pty2),2) + pow((ptz1-ptz2),2));
	if (distance < 1e-10)
		return 1;
	else 
		return 0;
}

void OneToOneSwept::buildAssociation(iGeom_Instance &geom, iMesh_Instance &mesh, iRel_Instance &assoc, iRel_PairHandle &rel)
{

    int err, namelen;

    SimpleArray<iBase_EntitySetHandle> entitySets;
    iBase_EntitySetHandle geom_root_set, mesh_root_set;

    // Get the root sets of the geometry and mesh.
    iMesh_getRootSet(mesh, &mesh_root_set, &err);
    iGeom_getRootSet(geom, &geom_root_set, &err);

    iBase_TagHandle geom_id_tag, mesh_id_tag, geom_dim_tag;
    const char *tag1 = "GLOBAL_ID";
    namelen = strlen(tag1);
    iGeom_getTagHandle(geom, tag1, &geom_id_tag, &err, namelen);
    iMesh_getTagHandle(mesh, tag1, &mesh_id_tag, &err, namelen);

    const char *tag2 = "GEOM_DIMENSION";
    namelen = strlen(tag2);
    iMesh_getTagHandle(mesh, tag2, &geom_dim_tag, &err, namelen);
    assert(!err);

    iRel_create(0, &assoc, &err, 0);

    iRel_createPair(assoc,
                    geom, iRel_ENTITY, iRel_IGEOM_IFACE, iRel_ACTIVE,
                    mesh, iRel_BOTH, iRel_IMESH_IFACE, iRel_ACTIVE, &rel, &err);
    assert(!err);

    // Get all the entitySet in the mesh
    iMesh_getEntSets(mesh, mesh_root_set, 0, ARRAY_INOUT(entitySets), &err);

    std::cout << "The size of all the EntitySets is " << entitySets.size() << std::endl;

    int ncount;
    iBase_EntityHandle gEntity;

    // Map all the geometric nodes
    SimpleArray<iBase_EntityHandle> gNodes;
    iGeom_getEntities(geom, geom_root_set, iBase_VERTEX, ARRAY_INOUT(gNodes), &err);
    assert(!err);
    
    int geom_id, geom_dim;
    std::map<int, iBase_EntityHandle> mapNodes;
    for (int i=0; i < gNodes.size(); i++)
    {
    	iGeom_getIntData(geom, gNodes[i], geom_id_tag, &geom_id, &err);
        mapNodes[geom_id] = gNodes[i];
    }
    

    // Map all the geometric edges.resize
    SimpleArray<iBase_EntityHandle> gEdges;
    iGeom_getEntities(geom, geom_root_set, iBase_EDGE, ARRAY_INOUT(gEdges), &err);
    assert(!err);

    
    std::map<int, iBase_EntityHandle> mapEdges;
    for (int i = 0; i < gEdges.size(); i++)
    {
        iGeom_getIntData(geom, gEdges[i], geom_id_tag, &geom_id, &err);
        mapEdges[geom_id] = gEdges[i];
    }

    // Map all the geometric faces ...
    SimpleArray<iBase_EntityHandle> gFaces;
    iGeom_getEntities(geom, geom_root_set, iBase_FACE, ARRAY_INOUT(gFaces), &err);
    assert(!err);

    std::map<int, iBase_EntityHandle> mapFaces;
    for (int i = 0; i < gFaces.size(); i++)
    {
        iGeom_getIntData(geom, gFaces[i], geom_id_tag, &geom_id, &err);
        mapFaces[geom_id] = gFaces[i];
    }

    // Map all the geometric cells ...
    SimpleArray<iBase_EntityHandle> gCells;
    iGeom_getEntities(geom, geom_root_set, iBase_REGION, ARRAY_INOUT(gCells), &err);
    assert(!err);

    std::map<int, iBase_EntityHandle> mapCells;
    for (int i = 0; i < gCells.size(); i++)
    {
        iGeom_getIntData(geom, gCells[i], geom_id_tag, &geom_id, &err);
        mapCells[geom_id] = gCells[i];
    }
    
    ///////////////////////////////////////////////////////////////////////////////
    // Create Vertex Assocations:
    ///////////////////////////////////////////////////////////////////////////////
    cout << " Building Vertex Associations " << endl;
    int numNodes;
    iMesh_getNumOfType(mesh, mesh_root_set, iBase_VERTEX, &numNodes, &err);
    assert(!err);
    SimpleArray<iBase_EntityHandle> mNodes;
    
    int numAssociations = 0;
    
    for (int i = 0; i < entitySets.size(); i++)
    {
        mNodes.clear();
        iMesh_getEntities(mesh, entitySets[i], iBase_VERTEX, iMesh_ALL_TOPOLOGIES,
                          ARRAY_INOUT(mNodes), &err);

        if (mNodes.size() && (mNodes.size() != numNodes))
        {
            ncount += mNodes.size();

            iMesh_getEntSetIntData(mesh, entitySets[i], mesh_id_tag, &geom_id, &err);
            assert(!err);

            iMesh_getEntSetIntData(mesh, entitySets[i], geom_dim_tag, &geom_dim, &err);
            assert(!err);

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
                iRel_setEntSetRelation(assoc, rel, gEntity, entitySets[i], &err);
        }
    }

    
    
    ///////////////////////////////////////////////////////////////////////////////
    // Create Edge Assocations:
    ///////////////////////////////////////////////////////////////////////////////
    cout << " Building Edge Associations " << endl;

    int numEdges;
    iMesh_getNumOfType(mesh, mesh_root_set, iBase_EDGE, &numEdges, &err);

    SimpleArray<iBase_EntityHandle> mEdges;

    numAssociations = 0;
    for (int i = 0; i < entitySets.size(); i++)
    {
        mEdges.clear();
        iMesh_getEntities(mesh, entitySets[i], iBase_EDGE, iMesh_ALL_TOPOLOGIES,
                          ARRAY_INOUT(mEdges), &err);

        if (mEdges.size() && (mEdges.size() != numEdges))
        {
            ncount += mEdges.size();

            iMesh_getEntSetIntData(mesh, entitySets[i], mesh_id_tag, &geom_id, &err);
            assert(!err);

            iMesh_getEntSetIntData(mesh, entitySets[i], geom_dim_tag, &geom_dim, &err);
            assert(!err);
	    
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
                iRel_setEntSetRelation(assoc, rel, gEntity, entitySets[i], &err);
        }
    }

    if (numAssociations != int(mapEdges.size()))
        cout << "Warning: There are more edge entitySet than geometric edges " << endl;
    
    //////////////////////////////////////////////////////////////////////////////
    // Face Association
    //////////////////////////////////////////////////////////////////////////////
    cout << " Building Face Associations " << endl;

    SimpleArray<iBase_EntityHandle> mFaces;

    int numFaces;
    iMesh_getNumOfType(mesh, mesh_root_set, iBase_FACE, &numFaces, &err);



    numAssociations = 0;
    for (int i = 0; i < entitySets.size(); i++)
    {
        mFaces.clear();
        iMesh_getEntities(mesh, entitySets[i], iBase_FACE, iMesh_ALL_TOPOLOGIES,
                          ARRAY_INOUT(mFaces), &err);

        if (mFaces.size() && (mFaces.size() != numFaces))
        {
            ncount += mFaces.size();
            iMesh_getEntSetIntData(mesh, entitySets[i], mesh_id_tag, &geom_id, &err);
            assert(!err);
            iMesh_getEntSetIntData(mesh, entitySets[i], geom_dim_tag, &geom_dim, &err);
            assert(!err);

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
                {
                    cout << "Fatal Error: Geometric face not found " << geom_id << endl;
                    exit(0);
                }
                break;
            }
            if (gEntity)
	    {
                iRel_setEntSetRelation(assoc, rel, gEntity, entitySets[i], &err);
		assert(!err);
	    }        
	}
    }
    
    if (numAssociations != int(mapFaces.size()))
    {
	cout << "Warning: There are more face entitySet than geometric faces " << endl;
    }
	
    //////////////////////////////////////////////////////////////////////////////
    // Cell Association
    //////////////////////////////////////////////////////////////////////////////

    SimpleArray<iBase_EntityHandle> mCells;

   

    int numCells;
    iMesh_getNumOfType(mesh, mesh_root_set, iBase_REGION, &numCells, &err);

    ncount = 0;
    for (int i = 0; i < entitySets.size(); i++)
    {
        mCells.clear();
        iMesh_getEntities(mesh, entitySets[i], iBase_REGION, iMesh_ALL_TOPOLOGIES,
                          ARRAY_INOUT(mCells), &err);

        if (mCells.size() && (mCells.size() != numCells))
        {
            ncount += mCells.size();
            iMesh_getEntSetIntData(mesh, entitySets[i], mesh_id_tag, &geom_id, &err);
            assert(!err);

            if (mapCells.find(geom_id) != mapCells.end())
            {
                if (mapCells.find(geom_id) != mapCells.end())
                {
                    gEntity = mapCells[geom_id];
                    iRel_setEntSetRelation(assoc, rel, gEntity, entitySets[i], &err);
                }
            }
        }
    }
}

int OneToOneSwept::SaveMesh(const char *FileName)
{
	int namelen = strlen(FileName), err;
    	iMesh_save(mesh, mesh_root_set, FileName, NULL, &err, namelen, 0);
	assert(!err);
	if(err)
	{
		return 1;//fail
	}
	else
	{
		return 0;//succeed
	}
	
}

int OneToOneSwept::InnerLayerMeshing()
{
	//first discretize the linking sides
	//temporarily discretize the linking sides based on the equally spacing
	int err, linking_num_pts=-1;
	
	
	//create the vertex on the linking surface
	//determine whether the different linking sides have the same number of layers
	for (unsigned int i=0; i < gLinkSides.size(); i++)
	{
		iBase_EntitySetHandle testset;
		iRel_getEntSetRelation(assoc, rel, gLinkSides[i].gEdgeHandle, 0, &testset, &err);
		if (!err)
		{//there is the mesh for the linking sides
			SimpleArray<iBase_EntityHandle> lvertices;
			iMesh_getEntities(mesh, testset, iBase_VERTEX, iMesh_ALL_TOPOLOGIES, ARRAY_INOUT(lvertices), &err);
			assert(!err);
			if (linking_num_pts==-1)
				linking_num_pts = lvertices.size();
			else
			{
				if (lvertices.size()!=linking_num_pts)
				{//there are the inconsistent mesh for the different linking surfs
					cout << "There are inconsistent meshes on the different linking surfs\nThe sweeeping algorithm will terminate" << endl;
					exit(1);
				}
			}
			numLayers = lvertices.size()+1;
		}
	}
	
	assert((gsEdgeList.size()==gtEdgeList.size())&&(gsEdgeList.size()==gLinkFaceList.size()));
	vector<vector <Vertex> > linkVertexList(numLayers-1, vector<Vertex>(NodeList.size()));
	
	LinkSurfMeshing(linkVertexList);
	
	//create the inner vertex on the different layers
	InnerNodesProjection(linkVertexList);
	
	
	//create the quadrilateral face elements on the different layers and cell elements
	CreateElements(linkVertexList); 	
	
	return 1;				
}

int OneToOneSwept::CreateElements(vector<vector <Vertex> > &linkVertexList)
{
	//create the quadrilateral elements on the different layers
	//it is not necessary to create the quadrilateral elements on the different layers. Hex element can be created based on the eight nodes
	int err, status;
	
	vector<iBase_EntityHandle> mVolumeHandle(FaceList.size());
		
	for (int m=0; m < numLayers-1; m++)
	{
		if (m==0)
		{
			for (unsigned int i=0; i < FaceList.size(); i++)
			{
				vector<iBase_EntityHandle> connect(8);
				
				connect[0] = NodeList[(FaceList[i].getVertex(0))->index].gVertexHandle;
				connect[1] = NodeList[(FaceList[i].getVertex(1))->index].gVertexHandle;
				connect[2] = linkVertexList[m][(FaceList[i].getVertex(1))->index].gVertexHandle;
				connect[3] = linkVertexList[m][(FaceList[i].getVertex(0))->index].gVertexHandle;
				
				connect[4] = NodeList[(FaceList[i].getVertex(3))->index].gVertexHandle;
				connect[5] = NodeList[(FaceList[i].getVertex(2))->index].gVertexHandle;
				connect[6] = linkVertexList[m][(FaceList[i].getVertex(2))->index].gVertexHandle;
				connect[7] = linkVertexList[m][(FaceList[i].getVertex(3))->index].gVertexHandle;
				iMesh_createEnt(mesh, iMesh_HEXAHEDRON, &connect[0], 8, &mVolumeHandle[i], &status, &err);
				assert(!err);
			}
			iMesh_addEntArrToSet(mesh, &mVolumeHandle[0], FaceList.size(), volumeSet, &err);
			assert(!err);
			if (m == (numLayers-2))
			{
				for (unsigned int i=0; i < FaceList.size(); i++)
				{
					vector<iBase_EntityHandle> connect(8);
					
					connect[0] = linkVertexList[m][(FaceList[i].getVertex(0))->index].gVertexHandle;
					connect[1] = linkVertexList[m][(FaceList[i].getVertex(1))->index].gVertexHandle;
					connect[2] = TVertexList[(FaceList[i].getVertex(1))->index].gVertexHandle;
					connect[3] = TVertexList[(FaceList[i].getVertex(0))->index].gVertexHandle;
				
					connect[4] = linkVertexList[m][(FaceList[i].getVertex(3))->index].gVertexHandle;
					connect[5] = linkVertexList[m][(FaceList[i].getVertex(2))->index].gVertexHandle;
					connect[6] = TVertexList[(FaceList[i].getVertex(2))->index].gVertexHandle;
					connect[7] = TVertexList[(FaceList[i].getVertex(3))->index].gVertexHandle;
					
					iMesh_createEnt(mesh, iMesh_HEXAHEDRON, &connect[0], 8, &mVolumeHandle[i], &status, &err);
					assert(!err);				
				}
				iMesh_addEntArrToSet(mesh, &mVolumeHandle[0], FaceList.size(), volumeSet, &err);
				assert(!err);
			}
		}
		else
		{
			
			for (unsigned int i=0; i < FaceList.size(); i++)
			{
				vector<iBase_EntityHandle> connect(8);
				
				connect[0] = linkVertexList[m-1][(FaceList[i].getVertex(0))->index].gVertexHandle;
				connect[1] = linkVertexList[m-1][(FaceList[i].getVertex(1))->index].gVertexHandle;
				connect[2] = linkVertexList[m][(FaceList[i].getVertex(1))->index].gVertexHandle;
				connect[3] = linkVertexList[m][(FaceList[i].getVertex(0))->index].gVertexHandle;
				
				connect[4] = linkVertexList[m-1][(FaceList[i].getVertex(3))->index].gVertexHandle;
				connect[5] = linkVertexList[m-1][(FaceList[i].getVertex(2))->index].gVertexHandle;
				connect[6] = linkVertexList[m][(FaceList[i].getVertex(2))->index].gVertexHandle;
				connect[7] = linkVertexList[m][(FaceList[i].getVertex(3))->index].gVertexHandle;
				iMesh_createEnt(mesh, iMesh_HEXAHEDRON, &connect[0], 8, &mVolumeHandle[i], &status, &err);
				assert(!err);				
			}
			iMesh_addEntArrToSet(mesh, &mVolumeHandle[0], FaceList.size(), volumeSet, &err);
			assert(!err);
			if (m == (numLayers-2))
			{
				for (unsigned int i=0; i < FaceList.size(); i++)
				{
					vector<iBase_EntityHandle> connect(8);
					
					connect[0] = linkVertexList[m][(FaceList[i].getVertex(0))->index].gVertexHandle;
					connect[1] = linkVertexList[m][(FaceList[i].getVertex(1))->index].gVertexHandle;
					connect[2] = TVertexList[(FaceList[i].getVertex(1))->index].gVertexHandle;
					connect[3] = TVertexList[(FaceList[i].getVertex(0))->index].gVertexHandle;
				
					connect[4] = linkVertexList[m][(FaceList[i].getVertex(3))->index].gVertexHandle;
					connect[5] = linkVertexList[m][(FaceList[i].getVertex(2))->index].gVertexHandle;
					connect[6] = TVertexList[(FaceList[i].getVertex(2))->index].gVertexHandle;
					connect[7] = TVertexList[(FaceList[i].getVertex(3))->index].gVertexHandle;
					iMesh_createEnt(mesh, iMesh_HEXAHEDRON, &connect[0], 8, &mVolumeHandle[i], &status, &err);
					assert(!err);				
				}
				iMesh_addEntArrToSet(mesh, &mVolumeHandle[0], FaceList.size(), volumeSet, &err);
				assert(!err);
			}
		}
	}
	
	return 1;
}


int OneToOneSwept::InnerNodesProjection(vector<vector <Vertex> > &linkVertexList)
{
	int err, numPts=0;
	Point3D sPtsCenter = {0, 0, 0}, tPtsCenter={0, 0, 0};
	vector<Point3D> PtsCenter(numLayers-1);
	vector<Point3D> sBoundaryNodes(0), tBoundaryNodes(0);
	vector<vector <Point3D> > iBoundaryNodes(numLayers-1, vector<Point3D>(0));
	
	//calculate the center coordinates
	for (int i=0; i< numLayers-1; i++)
	{
		PtsCenter[i].px = 0;
		PtsCenter[i].py = 0;
		PtsCenter[i].pz = 0;
	}
	for (unsigned int i=0; i < NodeList.size(); i++)
	{
		if (NodeList[i].onBoundary || NodeList[i].onCorner)
		{
			sPtsCenter.px = sPtsCenter.px + NodeList[i].xyzCoords[0];
			sPtsCenter.py = sPtsCenter.py + NodeList[i].xyzCoords[1];
			sPtsCenter.pz = sPtsCenter.pz + NodeList[i].xyzCoords[2];
			
			tPtsCenter.px = tPtsCenter.px + TVertexList[i].xyzCoords[0];
			tPtsCenter.py = tPtsCenter.py + TVertexList[i].xyzCoords[1];
			tPtsCenter.pz = tPtsCenter.pz + TVertexList[i].xyzCoords[2];
			
			numPts++;
			
			sBoundaryNodes.resize(numPts);
			tBoundaryNodes.resize(numPts);
			
			sBoundaryNodes[numPts-1].px = NodeList[i].xyzCoords[0];
			sBoundaryNodes[numPts-1].py = NodeList[i].xyzCoords[1];
			sBoundaryNodes[numPts-1].pz = NodeList[i].xyzCoords[2];	
			
			tBoundaryNodes[numPts-1].px = TVertexList[i].xyzCoords[0];
			tBoundaryNodes[numPts-1].py = TVertexList[i].xyzCoords[1];
			tBoundaryNodes[numPts-1].pz = TVertexList[i].xyzCoords[2];
			
			for (int j=0; j< numLayers-1; j++)
			{
				iBoundaryNodes[j].resize(numPts);
				PtsCenter[j].px = PtsCenter[j].px + linkVertexList[j][i].xyzCoords[0];
				PtsCenter[j].py = PtsCenter[j].py + linkVertexList[j][i].xyzCoords[1];
				PtsCenter[j].pz = PtsCenter[j].pz + linkVertexList[j][i].xyzCoords[2];
				
				iBoundaryNodes[j][numPts-1].px = linkVertexList[j][i].xyzCoords[0];
				iBoundaryNodes[j][numPts-1].py = linkVertexList[j][i].xyzCoords[1];
				iBoundaryNodes[j][numPts-1].pz = linkVertexList[j][i].xyzCoords[2];
			}
		}
	}
	sPtsCenter.px = sPtsCenter.px/double(numPts);
	sPtsCenter.py = sPtsCenter.py/double(numPts);
	sPtsCenter.pz = sPtsCenter.pz/double(numPts);
	
	tPtsCenter.px = tPtsCenter.px/double(numPts);
	tPtsCenter.py = tPtsCenter.py/double(numPts);
	tPtsCenter.pz = tPtsCenter.pz/double(numPts);
	
	//calculate the center coordinates for the ith layer 
	for (int i=0; i< numLayers-1; i++)
	{
		PtsCenter[i].px = PtsCenter[i].px/double(numPts);
		PtsCenter[i].py = PtsCenter[i].py/double(numPts);
		PtsCenter[i].pz = PtsCenter[i].pz/double(numPts);
	}
	
	//loop over different layers
	for (int i=0; i < numLayers-1; i++)
	{
		double sA[3][3], tA[3][3];
		double stransMatrix[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}}, ttransMatrix[3][3]={{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
		double sInvMatrix[3][3], tInvMatrix[3][3];
		double sb[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}}, tb[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
		vector<Point3D> sBNodes(numPts);	//boundary nodes on the source surface
		vector<Point3D> tBNodes(numPts);	//boundary nodes on the target surface
		vector<vector <Point3D> > isBNodes(numLayers-1, vector<Point3D>(numPts)), itBNodes(numLayers-1, vector<Point3D>(numPts));  //boundary nodes on the different layer
		
		
		//transform the coordinates
		for (int j=0; j < numPts; j++)
		{
			//transform the coordinates on the source suface
			sBNodes[j].px = sBoundaryNodes[j].px;
			sBNodes[j].py = sBoundaryNodes[j].py;
			sBNodes[j].pz = sBoundaryNodes[j].pz;
			
			tBNodes[j].px = tBoundaryNodes[j].px;
			tBNodes[j].py = tBoundaryNodes[j].py;
			tBNodes[j].pz = tBoundaryNodes[j].pz;
			
			isBNodes[i][j].px = iBoundaryNodes[i][j].px;
			isBNodes[i][j].py = iBoundaryNodes[i][j].py;
			isBNodes[i][j].pz = iBoundaryNodes[i][j].pz;
			
			itBNodes[i][j].px = iBoundaryNodes[i][j].px;
			itBNodes[i][j].py = iBoundaryNodes[i][j].py;
			itBNodes[i][j].pz = iBoundaryNodes[i][j].pz;
			
			sBNodes[j].px = sBNodes[j].px - 2*sPtsCenter.px + PtsCenter[i].px;
			sBNodes[j].py = sBNodes[j].py - 2*sPtsCenter.py + PtsCenter[i].py;
			sBNodes[j].pz = sBNodes[j].pz - 2*sPtsCenter.pz + PtsCenter[i].pz;
			
			tBNodes[j].px = tBNodes[j].px - 2*tPtsCenter.px + PtsCenter[i].px;
			tBNodes[j].py = tBNodes[j].py - 2*tPtsCenter.py + PtsCenter[i].py;
			tBNodes[j].pz = tBNodes[j].pz - 2*tPtsCenter.pz + PtsCenter[i].pz;
					
			//transform the coordinates on the different layers
			isBNodes[i][j].px = isBNodes[i][j].px - sPtsCenter.px;
			isBNodes[i][j].py = isBNodes[i][j].py - sPtsCenter.py;
			isBNodes[i][j].pz = isBNodes[i][j].pz - sPtsCenter.pz;
			
			itBNodes[i][j].px = itBNodes[i][j].px - tPtsCenter.px;
			itBNodes[i][j].py = itBNodes[i][j].py - tPtsCenter.py;
			itBNodes[i][j].pz = itBNodes[i][j].pz - tPtsCenter.pz;
		}
		
		//calculate the temporary matrix
		for (int j=0; j < numPts; j++)
		{
			//first row entries in the temporary matrix
			stransMatrix[0][0] = stransMatrix[0][0] + sBNodes[j].px*sBNodes[j].px;
			stransMatrix[0][1] = stransMatrix[0][1] + sBNodes[j].px*sBNodes[j].py;
			stransMatrix[0][2] = stransMatrix[0][2] + sBNodes[j].px*sBNodes[j].pz;
			//second row entries in the temporary matrix
			stransMatrix[1][0] = stransMatrix[1][0] + sBNodes[j].py*sBNodes[j].px;
			stransMatrix[1][1] = stransMatrix[1][1] + sBNodes[j].py*sBNodes[j].py;
			stransMatrix[1][2] = stransMatrix[1][2] + sBNodes[j].py*sBNodes[j].pz;
			//third row entries in the temporary matrix
			stransMatrix[2][0] = stransMatrix[2][0] + sBNodes[j].pz*sBNodes[j].px;
			stransMatrix[2][1] = stransMatrix[2][1] + sBNodes[j].pz*sBNodes[j].py;
			stransMatrix[2][2] = stransMatrix[2][2] + sBNodes[j].pz*sBNodes[j].pz;
			
			//first row entries in the temporary matrix
			ttransMatrix[0][0] = ttransMatrix[0][0] + tBNodes[j].px*tBNodes[j].px;
			ttransMatrix[0][1] = ttransMatrix[0][1] + tBNodes[j].px*tBNodes[j].py;
			ttransMatrix[0][2] = ttransMatrix[0][2] + tBNodes[j].px*tBNodes[j].pz;
			//second row entries in the temporary matrix
			ttransMatrix[1][0] = ttransMatrix[1][0] + tBNodes[j].py*tBNodes[j].px;
			ttransMatrix[1][1] = ttransMatrix[1][1] + tBNodes[j].py*tBNodes[j].py;
			ttransMatrix[1][2] = ttransMatrix[1][2] + tBNodes[j].py*tBNodes[j].pz;
			//third row entries in the temporary matrix
			ttransMatrix[2][0] = ttransMatrix[2][0] + tBNodes[j].pz*tBNodes[j].px;
			ttransMatrix[2][1] = ttransMatrix[2][1] + tBNodes[j].pz*tBNodes[j].py;
			ttransMatrix[2][2] = ttransMatrix[2][2] + tBNodes[j].pz*tBNodes[j].pz;
			
			
			//first row entries in the temporary matrix
			sb[0][0] = sb[0][0] + isBNodes[i][j].px*sBNodes[j].px;
			sb[0][1] = sb[0][1] + isBNodes[i][j].px*sBNodes[j].py;
			sb[0][2] = sb[0][2] + isBNodes[i][j].px*sBNodes[j].pz;
			//second row entries in the temporary matrix
			sb[1][0] = sb[1][0] + isBNodes[i][j].py*sBNodes[j].px;
			sb[1][1] = sb[1][1] + isBNodes[i][j].py*sBNodes[j].py;
			sb[1][2] = sb[1][2] + isBNodes[i][j].py*sBNodes[j].pz;
			//third row entries in the temporary matrix
			sb[2][0] = sb[2][0] + isBNodes[i][j].pz*sBNodes[j].px;
			sb[2][1] = sb[2][1] + isBNodes[i][j].pz*sBNodes[j].py;
			sb[2][2] = sb[2][2] + isBNodes[i][j].pz*sBNodes[j].pz;
			
			//first row entries in the temporary matrix
			tb[0][0] = tb[0][0] + itBNodes[i][j].px*tBNodes[j].px;
			tb[0][1] = tb[0][1] + itBNodes[i][j].px*tBNodes[j].py;
			tb[0][2] = tb[0][2] + itBNodes[i][j].px*tBNodes[j].pz;
			//second row entries in the temporary matrix
			tb[1][0] = tb[1][0] + itBNodes[i][j].py*tBNodes[j].px;
			tb[1][1] = tb[1][1] + itBNodes[i][j].py*tBNodes[j].py;
			tb[1][2] = tb[1][2] + itBNodes[i][j].py*tBNodes[j].pz;
			//third row entries in the temporary matrix
			tb[2][0] = tb[2][0] + itBNodes[i][j].pz*tBNodes[j].px;
			tb[2][1] = tb[2][1] + itBNodes[i][j].pz*tBNodes[j].py;
			tb[2][2] = tb[2][2] + itBNodes[i][j].pz*tBNodes[j].pz;
			
		}		
		
		//first determine the affine mapping matrix is singular or not
		double sdetMatrix = stransMatrix[0][2]*stransMatrix[1][1]*stransMatrix[2][0] - stransMatrix[0][1]*stransMatrix[1][2]*stransMatrix[2][0] - stransMatrix[0][2]*stransMatrix[1][0]*stransMatrix[2][1] + stransMatrix[0][0]*stransMatrix[1][2]*stransMatrix[2][1] + stransMatrix[0][1]*stransMatrix[1][0]*stransMatrix[2][2] - stransMatrix[0][0]*stransMatrix[1][1]*stransMatrix[2][2];
		double tdetMatrix = ttransMatrix[0][2]*ttransMatrix[1][1]*ttransMatrix[2][0] - ttransMatrix[0][1]*ttransMatrix[1][2]*ttransMatrix[2][0] - ttransMatrix[0][2]*ttransMatrix[1][0]*ttransMatrix[2][1] + ttransMatrix[0][0]*ttransMatrix[1][2]*ttransMatrix[2][1] + ttransMatrix[0][1]*ttransMatrix[1][0]*ttransMatrix[2][2] - ttransMatrix[0][0]*ttransMatrix[1][1]*ttransMatrix[2][2];
		//transMatrix[0][0]*(transMatrix[1][1]*transMatrix[2][2]-transMatrix[2][1]*transMatrix[1][2]) - transMatrix[0][1]*(transMatrix[1][0]*transMatrix[2][2] - transMatrix[2][0]*transMatrix[1][2]) + transMatrix[0][2]*(transMatrix[1][0]*transMatrix[2][1]-transMatrix[1][1]*transMatrix[2][0]);
		assert(pow(sdetMatrix, 2)>1.0e-20);
		assert(pow(tdetMatrix, 2)>1.0e-20);
		
		////solve the affine mapping matrix, make use of inverse matrix to get affine mapping matrix
		sInvMatrix[0][0] = (stransMatrix[2][1]*stransMatrix[1][2] - stransMatrix[1][1]*stransMatrix[2][2])/sdetMatrix;
		
		sInvMatrix[0][1] = (stransMatrix[0][1]*stransMatrix[2][2] - stransMatrix[0][2]*stransMatrix[2][1])/sdetMatrix;
		
		sInvMatrix[0][2] = (stransMatrix[0][2]*stransMatrix[1][1] - stransMatrix[0][1]*stransMatrix[1][2])/sdetMatrix;
		
		sInvMatrix[1][0] = (stransMatrix[1][0]*stransMatrix[2][2] - stransMatrix[1][2]*stransMatrix[2][0])/sdetMatrix;
		
		sInvMatrix[1][1] = (stransMatrix[0][2]*stransMatrix[2][0] - stransMatrix[0][0]*stransMatrix[2][2])/sdetMatrix;
		
		sInvMatrix[1][2] = (stransMatrix[0][0]*stransMatrix[1][2] - stransMatrix[0][2]*stransMatrix[1][0])/sdetMatrix;
		
		sInvMatrix[2][0] = (stransMatrix[1][1]*stransMatrix[2][0] - stransMatrix[1][0]*stransMatrix[2][1])/sdetMatrix;
		
		sInvMatrix[2][1] = (stransMatrix[0][0]*stransMatrix[2][1] - stransMatrix[0][1]*stransMatrix[2][0])/sdetMatrix;
		
		sInvMatrix[2][2] = (stransMatrix[0][1]*stransMatrix[1][0] - stransMatrix[0][0]*stransMatrix[1][1])/sdetMatrix;
		
		//projection from the target surface
		tInvMatrix[0][0] = (ttransMatrix[2][1]*ttransMatrix[1][2] - ttransMatrix[1][1]*ttransMatrix[2][2])/tdetMatrix;
		
		tInvMatrix[0][1] = (ttransMatrix[0][1]*ttransMatrix[2][2] - ttransMatrix[0][2]*ttransMatrix[2][1])/tdetMatrix;
		
		tInvMatrix[0][2] = (ttransMatrix[0][2]*ttransMatrix[1][1] - ttransMatrix[0][1]*ttransMatrix[1][2])/tdetMatrix;
		
		tInvMatrix[1][0] = (ttransMatrix[1][0]*ttransMatrix[2][2] - ttransMatrix[1][2]*ttransMatrix[2][0])/tdetMatrix;
		
		tInvMatrix[1][1] = (ttransMatrix[0][2]*ttransMatrix[2][0] - ttransMatrix[0][0]*ttransMatrix[2][2])/tdetMatrix;
		
		tInvMatrix[1][2] = (ttransMatrix[0][0]*ttransMatrix[1][2] - ttransMatrix[0][2]*ttransMatrix[1][0])/tdetMatrix;
		
		tInvMatrix[2][0] = (ttransMatrix[1][1]*ttransMatrix[2][0] - ttransMatrix[1][0]*ttransMatrix[2][1])/tdetMatrix;
		
		tInvMatrix[2][1] = (ttransMatrix[0][0]*ttransMatrix[2][1] - ttransMatrix[0][1]*ttransMatrix[2][0])/tdetMatrix;
		
		tInvMatrix[2][2] = (ttransMatrix[0][1]*ttransMatrix[1][0] - ttransMatrix[0][0]*ttransMatrix[1][1])/tdetMatrix;
		
		
		sA[0][0] = sInvMatrix[0][0]*sb[0][0] + sInvMatrix[0][1]*sb[0][1] + sInvMatrix[0][2]*sb[0][2];
		sA[0][1] = sInvMatrix[1][0]*sb[0][0] + sInvMatrix[1][1]*sb[0][1] + sInvMatrix[1][2]*sb[0][2];
		sA[0][2] = sInvMatrix[2][0]*sb[0][0] + sInvMatrix[2][1]*sb[0][1] + sInvMatrix[2][2]*sb[0][2];
		
		sA[1][0] = sInvMatrix[0][0]*sb[1][0] + sInvMatrix[0][1]*sb[1][1] + sInvMatrix[0][2]*sb[1][2];
		sA[1][1] = sInvMatrix[1][0]*sb[1][0] + sInvMatrix[1][1]*sb[1][1] + sInvMatrix[1][2]*sb[1][2];
		sA[1][2] = sInvMatrix[2][0]*sb[1][0] + sInvMatrix[2][1]*sb[1][1] + sInvMatrix[2][2]*sb[1][2];
		
		sA[2][0] = sInvMatrix[0][0]*sb[2][0] + sInvMatrix[0][1]*sb[2][1] + sInvMatrix[0][2]*sb[2][2];
		sA[2][1] = sInvMatrix[1][0]*sb[2][0] + sInvMatrix[1][1]*sb[2][1] + sInvMatrix[1][2]*sb[2][2];
		sA[2][2] = sInvMatrix[2][0]*sb[2][0] + sInvMatrix[2][1]*sb[2][1] + sInvMatrix[2][2]*sb[2][2];
		
		tA[0][0] = tInvMatrix[0][0]*tb[0][0] + tInvMatrix[0][1]*tb[0][1] + tInvMatrix[0][2]*tb[0][2];
		tA[0][1] = tInvMatrix[1][0]*tb[0][0] + tInvMatrix[1][1]*tb[0][1] + tInvMatrix[1][2]*tb[0][2];
		tA[0][2] = tInvMatrix[2][0]*tb[0][0] + tInvMatrix[2][1]*tb[0][1] + tInvMatrix[2][2]*tb[0][2];
		
		tA[1][0] = tInvMatrix[0][0]*tb[1][0] + tInvMatrix[0][1]*tb[1][1] + tInvMatrix[0][2]*tb[1][2];
		tA[1][1] = tInvMatrix[1][0]*tb[1][0] + tInvMatrix[1][1]*tb[1][1] + tInvMatrix[1][2]*tb[1][2];
		tA[1][2] = tInvMatrix[2][0]*tb[1][0] + tInvMatrix[2][1]*tb[1][1] + tInvMatrix[2][2]*tb[1][2];
		
		tA[2][0] = tInvMatrix[0][0]*tb[2][0] + tInvMatrix[0][1]*tb[2][1] + tInvMatrix[0][2]*tb[2][2];
		tA[2][1] = tInvMatrix[1][0]*tb[2][0] + tInvMatrix[1][1]*tb[2][1] + tInvMatrix[1][2]*tb[2][2];
		tA[2][2] = tInvMatrix[2][0]*tb[2][0] + tInvMatrix[2][1]*tb[2][1] + tInvMatrix[2][2]*tb[2][2];
		
		//calculate the inner nodes for different layers
		for (unsigned int j=0; j < NodeList.size(); j++)
		{
			if (!(NodeList[j].onBoundary || NodeList[j].onCorner))
			{
				Point3D spts, tpts, pts;
				double s;
				iBase_EntityHandle nodeHandle;
				spts.px = sA[0][0]*(NodeList[j].xyzCoords[0] - 2*sPtsCenter.px + PtsCenter[i].px) + sA[0][1]*(NodeList[j].xyzCoords[1] - 2*sPtsCenter.py + PtsCenter[i].py) + sA[0][2]*(NodeList[j].xyzCoords[2] - 2*sPtsCenter.pz + PtsCenter[i].pz) + sPtsCenter.px;
				spts.py = sA[1][0]*(NodeList[j].xyzCoords[0] - 2*sPtsCenter.px + PtsCenter[i].px) + sA[1][1]*(NodeList[j].xyzCoords[1] - 2*sPtsCenter.py + PtsCenter[i].py) + sA[1][2]*(NodeList[j].xyzCoords[2] - 2*sPtsCenter.pz + PtsCenter[i].pz) + sPtsCenter.py;
				spts.pz = sA[2][0]*(NodeList[j].xyzCoords[0] - 2*sPtsCenter.px + PtsCenter[i].px) + sA[2][1]*(NodeList[j].xyzCoords[1] - 2*sPtsCenter.py + PtsCenter[i].py) + sA[2][2]*(NodeList[j].xyzCoords[2] - 2*sPtsCenter.pz + PtsCenter[i].pz) + sPtsCenter.pz;
				
				tpts.px = tA[0][0]*(TVertexList[j].xyzCoords[0] - 2*tPtsCenter.px + PtsCenter[i].px) + tA[0][1]*(TVertexList[j].xyzCoords[1] - 2*tPtsCenter.py + PtsCenter[i].py) + tA[0][2]*(TVertexList[j].xyzCoords[2] - 2*tPtsCenter.pz + PtsCenter[i].pz) + tPtsCenter.px;
				tpts.py = tA[1][0]*(TVertexList[j].xyzCoords[0] - 2*tPtsCenter.px + PtsCenter[i].px) + tA[1][1]*(TVertexList[j].xyzCoords[1] - 2*tPtsCenter.py + PtsCenter[i].py) + tA[1][2]*(TVertexList[j].xyzCoords[2] - 2*tPtsCenter.pz + PtsCenter[i].pz) + tPtsCenter.py;
				tpts.pz = tA[2][0]*(TVertexList[j].xyzCoords[0] - 2*tPtsCenter.px + PtsCenter[i].px) + tA[2][1]*(TVertexList[j].xyzCoords[1] - 2*tPtsCenter.py + PtsCenter[i].py) + tA[2][2]*(TVertexList[j].xyzCoords[2] - 2*tPtsCenter.pz + PtsCenter[i].pz) + tPtsCenter.pz;
				
				s = (i+1)*1/double(numLayers);
				pts.px = linear_interpolation(s, spts.px, tpts.px);
				pts.py = linear_interpolation(s, spts.py, tpts.py);
				pts.pz = linear_interpolation(s, spts.pz, tpts.pz);
				
				linkVertexList[i][j].xyzCoords[0] = pts.px;
				linkVertexList[i][j].xyzCoords[1] = pts.py;
				linkVertexList[i][j].xyzCoords[2] = pts.pz;
				
				
				iMesh_createVtx(mesh, pts.px, pts.py, pts.pz, &nodeHandle, &err);
				assert(!err);	
				linkVertexList[i][j].gVertexHandle = nodeHandle;
				iMesh_addEntToSet(mesh, nodeHandle, volumeSet, &err);
				assert(!err);			
			}	
		}								
	}
	
	return 1;
}

int OneToOneSwept::LinkSurfMeshing(vector<vector <Vertex> > &linkVertexList)
{
	//discretize the linking sides
	cout << "Linking surface meshing..." << endl;
	int err, LineIndex=0, faceIndex=0;
	vector<iBase_EntityHandle> lineHandle(0), faceHandle(0);
	iBase_TagHandle taghandle;
	iMesh_getTagHandle(mesh, "source", &taghandle, &err, strlen("source"));
	assert(!err);
	for (unsigned int i=0; i < gLinkSides.size(); i++)
	{
		int sIndex, stIndex;
		double u0, u1;
		vector<iBase_EntityHandle> mNodeHandle(numLayers-1);
		//find the node
		for (unsigned int j=0; j < gVertexList.size(); j=j+2)
		{
			if ((gLinkSides[i].connect[0]->id == gVertexList[j].id)||(gLinkSides[i].connect[1]->id == gVertexList[j].id))
			{
				sIndex = j;
				break;
			}
		}
		
		//find the corresponding mesh vertex on the target surface in order to get index for mesh vertex
		iBase_EntitySetHandle CornerSet;
		SimpleArray<iBase_EntityHandle> tmpNodes;
		iRel_getEntSetRelation(assoc, rel, gVertexList[sIndex].gVertexHandle, 0, &CornerSet, &err);
		assert(!err);
		iMesh_getEntities(mesh, CornerSet, iBase_VERTEX, iMesh_POINT, ARRAY_INOUT(tmpNodes), &err);
		assert(!err);
		assert(tmpNodes.size()==1);
		iMesh_getIntData(mesh, tmpNodes[0], taghandle, &stIndex, &err);
		assert(!err);
		
	
		iGeom_getEntXYZtoU(geom, gLinkSides[i].gEdgeHandle, gVertexList[sIndex].xyzCoords[0], gVertexList[sIndex].xyzCoords[1], gVertexList[sIndex].xyzCoords[2], &u0, &err);
		assert(!err);
		iGeom_getEntXYZtoU(geom, gLinkSides[i].gEdgeHandle, gVertexList[cornerPairs[sIndex]].xyzCoords[0], gVertexList[cornerPairs[sIndex]].xyzCoords[1], gVertexList[cornerPairs[sIndex]].xyzCoords[2], &u1, &err);
		assert(!err);
		
		
		/////////////////////////////////////////////////////////////////////////////////////////////////////////
		//detect whether linking side[i] is discretized or not
		iBase_EntitySetHandle TestSet;
		iRel_getEntSetRelation(assoc, rel, gLinkSides[i].gEdgeHandle, 0, &TestSet, &err);
		if (!err)
		{//linking side gLinkSides[i] has been discretized
			SimpleArray<iBase_EntityHandle> vertices;
			iMesh_getEntities(mesh, TestSet, iBase_VERTEX, iMesh_ALL_TOPOLOGIES, ARRAY_INOUT(vertices), &err);
			assert(!err);
			double xyz1[3], xyz2[3], dist1, dist2;
			bool isReverse = false;
			
			if (vertices.size()>1)
			{
				iMesh_getVtxCoord(mesh, vertices[0], &xyz1[0], &xyz1[1], &xyz1[2], &err);
				assert(!err);
				iMesh_getVtxCoord(mesh, vertices[1], &xyz2[0], &xyz2[1], &xyz2[2], &err);
				assert(!err);
				dist1 = sqrt( pow(NodeList[stIndex].xyzCoords[0]-xyz1[0], 2) + pow( NodeList[stIndex].xyzCoords[1]-xyz1[1], 2) + pow(NodeList[stIndex].xyzCoords[2]-xyz1[2], 2));
				dist2 = sqrt( pow(NodeList[stIndex].xyzCoords[0]-xyz2[0], 2) + pow( NodeList[stIndex].xyzCoords[1]-xyz2[1], 2) + pow(NodeList[stIndex].xyzCoords[2]-xyz2[2], 2));
				if (dist1 > dist2)
					isReverse = true;
			}
			
			for (int m=0; m < vertices.size(); m++)
			{
				iMesh_getVtxCoord(mesh, vertices[m], &xyz1[0], &xyz1[1], &xyz1[2], &err);
				assert(!err);
				if (!isReverse)
				{
					linkVertexList[m][stIndex].xyzCoords[0] = xyz1[0];
					linkVertexList[m][stIndex].xyzCoords[1] = xyz1[1];
					linkVertexList[m][stIndex].xyzCoords[2] = xyz1[2];
					linkVertexList[m][stIndex].index = stIndex;
					linkVertexList[m][stIndex].gVertexHandle = vertices[m];
				}
				else
				{
					linkVertexList[numLayers-m-2][stIndex].xyzCoords[0] = xyz1[0];
					linkVertexList[numLayers-m-2][stIndex].xyzCoords[1] = xyz1[1];
					linkVertexList[numLayers-m-2][stIndex].xyzCoords[2] = xyz1[2];
					linkVertexList[numLayers-m-2][stIndex].index = stIndex;
					linkVertexList[numLayers-m-2][stIndex].gVertexHandle = vertices[m];
				}				
			}
			continue;
		}
		/////////////////////////////////////////////////////////////////////////////////////////////////////////
		
		
		for (int j=0; j < (numLayers-1); j++)
		{
			//discretize the linking sides
			double t, u;
			Point3D pts;
			
			t=(j+1)*1/double(numLayers);
			u = linear_interpolation(t, u0, u1);
			iGeom_getEntUtoXYZ(geom, gLinkSides[i].gEdgeHandle, u, &pts.px, &pts.py, &pts.pz, &err);
			assert(!err);
			iMesh_createVtx(mesh, pts.px, pts.py, pts.pz, &mNodeHandle[j], &err);
			assert(!err);
			
			//add the new generated vertex to the list
			linkVertexList[j][stIndex].xyzCoords[0] = pts.px;
			linkVertexList[j][stIndex].xyzCoords[1] = pts.py;
			linkVertexList[j][stIndex].xyzCoords[2] = pts.pz;
			linkVertexList[j][stIndex].index = stIndex;
			linkVertexList[j][stIndex].gVertexHandle = mNodeHandle[j];	
		}
		
		//create the line segments on the linking sides
		vector<iBase_EntityHandle> connect(2*numLayers);
		SimpleArray<iBase_EntityHandle> edgeHandle;
		SimpleArray<int> status;
		connect[0] = NodeList[stIndex].gVertexHandle;
		for (int j=0; j < (numLayers-1); j++)
		{
			connect[2*j+1] = mNodeHandle[j];
			connect[2*j+2] = mNodeHandle[j];		
		}
		connect[2*numLayers-1] = TVertexList[stIndex].gVertexHandle;
		iMesh_createEntArr(mesh, iMesh_LINE_SEGMENT, &connect[0], 2*numLayers, ARRAY_INOUT(edgeHandle), ARRAY_INOUT(status), &err);
		assert(!err);
		
		//create the entityset	for the gLinkSides[i]
		iBase_EntitySetHandle entityset;
		iRel_getEntSetRelation(assoc, rel, gLinkSides[i].gEdgeHandle, 0, &entityset, &err);
		if (err)  //there is no entityset associated with gLinkSides[i].gEdgeHandle
		{
			iMesh_createEntSet(mesh, 1, &entityset, &err);
			assert(!err);
		}
		iMesh_addEntArrToSet(mesh, &mNodeHandle[0], numLayers-1, entityset, &err);
		assert(!err);
		iMesh_addEntArrToSet(mesh, &edgeHandle[0], edgeHandle.size(), entityset, &err);
		
		//create the irel between the linking sides and entityset
		iRel_getEntSetRelation(assoc, rel, gLinkSides[i].gEdgeHandle, 0, &entityset, &err);
		if (err)  //there is no entityset associated with gLinkSides[i].gEdgeHandle
		{
			iRel_setEntSetRelation(assoc, rel, gLinkSides[i].gEdgeHandle, entityset, &err);
			assert(!err);
		}
	}
	//finish the discretization of linking sides
	
	
	
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//check whether this linking surf is meshed or not
	vector<bool> isMeshed(gLinkFaceList.size());
	for (unsigned int i=0; i < gLinkFaceList.size(); i++)
	{
		int t1=1, t2=1, t3=1;
		int sEdgeIndex, gLeftIndex, gRightIndex;
		iBase_EntitySetHandle tmpSet;
		SimpleArray<iBase_EntityHandle> aEdges, aNodes, aFaces, gEdges;
		SimpleArray<iBase_EntityHandle> bNodes1, bNodes2, bNodes3;
		iBase_EntitySetHandle aEdgeSet;
		iBase_EntityHandle preNode, nextNode;
		vector<iBase_EntityHandle> Connect;
		double xyzCoords[3];
		
		gEdges.clear();
		iGeom_getEntAdj(geom, gLinkFaceList[i].gFaceHandle, iBase_EDGE, ARRAY_INOUT(gEdges), &err);
		assert(!err);
		//find the corresponding edge on the linking surface, which is shared by source surface
		for (int j=0; j< gEdges.size(); j++)
		{
			int edgeId;
			
			iGeom_getIntData(geom, gEdges[j], geom_id_tag, &edgeId, &err);
			assert(!err);
			//actually, the num of linking sides is equal to the num of edges on the source surface
			for (unsigned int k=0; ((k < gsEdgeList.size())&&(t1)); k++)
			{
				if (edgeId == gsEdgeList[k].EdgeID)
				{
					sEdgeIndex=k;
					t1=0;
					break;
				}
			}
			for (unsigned int k=0; (k < gLinkSides.size())&&(t2||t3); k++)
			{
				if ((edgeId == gLinkSides[k].EdgeID)&&(t2))
				{
					gLeftIndex=k;
					t2=0;
					break;
				}
				if ((edgeId == gLinkSides[k].EdgeID)&&(t3))
				{
					gRightIndex=k;
					t3=0;
					break;
				}
			}			
		}
		if (!((gsEdgeList[sEdgeIndex].connect[0]->id==gLinkSides[gLeftIndex].connect[0]->id)||(gsEdgeList[sEdgeIndex].connect[0]->id==gLinkSides[gLeftIndex].connect[1]->id)))
		{
			int temp=gRightIndex;
			gRightIndex = gLeftIndex;
			gLeftIndex = temp;
		}
		
		
		
		iRel_getEntSetRelation(assoc, rel, gLinkFaceList[i].gFaceHandle, 0, &tmpSet, &err);
		if(!err)
		{
			isMeshed[i] = true;
			
			iRel_getEntSetRelation(assoc, rel, gsEdgeList[sEdgeIndex].gEdgeHandle, 0, &aEdgeSet, &err);
			assert(!err);
			
			iMesh_getEntities(mesh, aEdgeSet, iBase_VERTEX, iMesh_ALL_TOPOLOGIES, ARRAY_INOUT(aNodes), &err);
			assert(!err);
			
			
			for (int mm=0; mm < aNodes.size(); mm++)
			{
				int VertexId;
				iMesh_getIntData(mesh, aNodes[mm], taghandle, &VertexId, &err);
				assert(!err);
				
				nextNode = aNodes[mm];
				preNode = aNodes[mm];				
				
				iMesh_getVtxCoord(mesh, nextNode, &xyzCoords[0], &xyzCoords[1], &xyzCoords[2], &err);
				assert(!err);
				
				for (int m2=0; m2 < numLayers-1; m2++)
				{					
					int tmpId;
					aFaces.clear();
					iMesh_getEntAdj(mesh, nextNode, iBase_FACE, ARRAY_INOUT(aFaces), &err);
					assert(!err);
					
					iMesh_getIntData(mesh, nextNode, mesh_id_tag, &tmpId, &err);
					assert(!err);				
					
					int tmpIndex=0;
					
					Connect.clear();
					if (m2==0)
					{
						tmpIndex = 0;
						for (int m1=0; m1 < aFaces.size(); m1++)
						{
							int is_contained;
							iMesh_isEntContained(mesh, tmpSet, aFaces[m1], &is_contained, &err);
							assert(!err);
							if (is_contained)
							{
								tmpIndex++;
								Connect.resize(tmpIndex);
								Connect[tmpIndex-1] = aFaces[m1];
							}						
						}
					}
					else
					{
						tmpIndex = 0;
						Connect.clear();
						int preIndex;
						
						iMesh_getIntData(mesh, preNode, mesh_id_tag, &preIndex, &err);
						assert(!err);
						
						for (int m1=0; m1 < aFaces.size(); m1++)
						{
							bool found = false;
							bNodes3.clear();
							iMesh_getEntAdj(mesh, aFaces[m1], iBase_VERTEX, ARRAY_INOUT(bNodes3), &err);
							assert(!err);
							
							
							for (int m3=0; m3 < bNodes3.size(); m3++)
							{
								int tmp;
								iMesh_getIntData(mesh, bNodes3[m3], mesh_id_tag, &tmp, &err);
								assert(!err);
								
								iMesh_getVtxCoord(mesh, bNodes3[m3], &xyzCoords[0], &xyzCoords[1], &xyzCoords[2], &err);
								assert(!err);
								if (preIndex == tmp)
								{
									found = true;
									break;
								}
								
							}
							if (!found)
							{
								tmpIndex++;
								Connect.resize(tmpIndex);
								Connect[tmpIndex-1] = aFaces[m1];
							}
						}
					}
					
					bNodes1.clear();
					iMesh_getEntAdj(mesh, Connect[0], iBase_VERTEX, ARRAY_INOUT(bNodes1), &err);
					assert(!err);
					bNodes2.clear();
					iMesh_getEntAdj(mesh, Connect[1], iBase_VERTEX, ARRAY_INOUT(bNodes2), &err);
					assert(!err);
					assert((bNodes1.size()==bNodes1.size())&&(bNodes1.size()==4));
					
					for (int m1=0; m1 < bNodes1.size(); m1++)
					{
						int tmp;
						iMesh_getIntData(mesh, bNodes1[m1], mesh_id_tag, &tmp, &err);
						assert(!err);
						if (tmpId==tmp)
						{
							tmpIndex = m1;
							Connect[0] = bNodes1[m1];
							break;
						}
					}
					
					int option1= (tmpIndex+1)%4, option2= (3+tmpIndex)%4, index1, index2;
					iMesh_getIntData(mesh, bNodes1[option1], mesh_id_tag, &index1, &err);
					assert(!err);
					iMesh_getIntData(mesh, bNodes1[option2], mesh_id_tag, &index2, &err);
					assert(!err);
					
					for (int m1 =0; m1 < bNodes2.size(); m1++)
					{
						int tmp;
						iMesh_getIntData(mesh, bNodes2[m1], mesh_id_tag, &tmp, &err);
						assert(!err);
						
						if (index1==tmp)
						{
							Connect[1] = bNodes2[m1];
							linkVertexList[m2][VertexId].gVertexHandle = bNodes2[m1];
							linkVertexList[m2][VertexId].index = tmp;
							
							iMesh_getVtxCoord(mesh, linkVertexList[m2][VertexId].gVertexHandle, &linkVertexList[m2][VertexId].xyzCoords[0], &linkVertexList[m2][VertexId].xyzCoords[1], &linkVertexList[m2][VertexId].xyzCoords[2], &err);							
							assert(!err);
							break;
						}
						if (index2==tmp)
						{
							Connect[1] = bNodes2[m1];
							linkVertexList[m2][VertexId].gVertexHandle = bNodes2[m1];
							linkVertexList[m2][VertexId].index = tmp;
							
							iMesh_getVtxCoord(mesh, linkVertexList[m2][VertexId].gVertexHandle, &linkVertexList[m2][VertexId].xyzCoords[0], &linkVertexList[m2][VertexId].xyzCoords[1], &linkVertexList[m2][VertexId].xyzCoords[2], &err);
							assert(!err);
							break;
						}
					}					
					//preNode = Connect[1];
					nextNode = Connect[1];
					if (m2 > 0)
						preNode = linkVertexList[m2-1][VertexId].gVertexHandle;
				}
				
			}
		}
		else
		{
			isMeshed[i] = false;	
		}
		
	
	
	}
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	//discretize the linking surface		
	for (unsigned int i=0; i < gLinkFaceList.size(); i++)
	{
		//find the linking sides
		if (!isMeshed[i])
		{
			int t1=1, t2=1, t3=1, leftIndex, RightIndex, nodehandleIndex=0;//temporary variable for controlling the loop
			//leftIndex and RightIndex are used to store the index for mesh corners
			double suLeft, suRight;
			//suLeft------parametrical coordinatiMesh_addEntArrToSet(mesh, &mNodeHandle[0], numLayers-1, entityset, &err);es on the source boundary edge
			Point2D pt00, pt01, pt10, pt11;  //parametric coordinates for four corners on the linking surface
			SimpleArray<iBase_EntityHandle> gEdges, mNodes;
			int sEdgeIndex, gLeftIndex, gRightIndex;
			iBase_EntitySetHandle edgeNodeSet;
			vector<iBase_EntityHandle>  nodeHandle(0);
		
			gEdges.clear();
			iGeom_getEntAdj(geom, gLinkFaceList[i].gFaceHandle, iBase_EDGE, ARRAY_INOUT(gEdges), &err);
			assert(!err);
			//find the corresponding edge on the linking surface, which is shared by source surface
			for (int j=0; j< gEdges.size(); j++)
			{
				int edgeId;
			
				iGeom_getIntData(geom, gEdges[j], geom_id_tag, &edgeId, &err);
				assert(!err);
				//actually, the num of linking sides is equal to the num of edges on the source surface
				for (unsigned int k=0; ((k < gsEdgeList.size())&&(t1)); k++)
				{
					if (edgeId == gsEdgeList[k].EdgeID)
					{
						sEdgeIndex=k;
						t1=0;
						break;
					}
				}
				for (unsigned int k=0; (k < gLinkSides.size())&&(t2||t3); k++)
				{
					if ((edgeId == gLinkSides[k].EdgeID)&&(t2))
					{
						gLeftIndex=k;
						t2=0;
						break;
					}
					if ((edgeId == gLinkSides[k].EdgeID)&&(t3))
					{
						gRightIndex=k;
						t3=0;
						break;
					}
				}			
			}
			if (!((gsEdgeList[sEdgeIndex].connect[0]->id==gLinkSides[gLeftIndex].connect[0]->id)||(gsEdgeList[sEdgeIndex].connect[0]->id==gLinkSides[gLeftIndex].connect[1]->id)))
			{
				int temp=gRightIndex;
				gRightIndex = gLeftIndex;
				gLeftIndex = temp;
			}
		
			//parametrize the linking surface
			//loop over the layers
			iRel_getEntSetRelation(assoc, rel, gsEdgeList[sEdgeIndex].gEdgeHandle, 0, &edgeNodeSet, &err);
			assert(!err);
			mNodes.clear();
			iMesh_getEntities(mesh, edgeNodeSet, iBase_VERTEX, iMesh_POINT, ARRAY_INOUT(mNodes), &err);
			assert(!err);
		
			//get the parametric coordinates for four corners on the linking surface
			iGeom_getEntXYZtoUV(geom, gLinkFaceList[i].gFaceHandle, gsEdgeList[sEdgeIndex].connect[0]->xyzCoords[0], gsEdgeList[sEdgeIndex].connect[0]->xyzCoords[1], gsEdgeList[sEdgeIndex].connect[0]->xyzCoords[2], &pt00.pu, &pt00.pv, &err);
			assert(!err);
			iGeom_getEntXYZtoUV(geom, gLinkFaceList[i].gFaceHandle, gsEdgeList[sEdgeIndex].connect[1]->xyzCoords[0], gsEdgeList[sEdgeIndex].connect[1]->xyzCoords[1], gsEdgeList[sEdgeIndex].connect[1]->xyzCoords[2], &pt10.pu, &pt10.pv, &err);
			assert(!err);

			iGeom_getEntXYZtoUV(geom, gLinkFaceList[i].gFaceHandle, gVertexList[cornerPairs[gsEdgeList[sEdgeIndex].connect[0]->index]].xyzCoords[0], gVertexList[cornerPairs[gsEdgeList[sEdgeIndex].connect[0]->index]].xyzCoords[1], gVertexList[cornerPairs[gsEdgeList[sEdgeIndex].connect[0]->index]].xyzCoords[2], &pt01.pu, &pt01.pv, &err);
			assert(!err);
			iGeom_getEntXYZtoUV(geom, gLinkFaceList[i].gFaceHandle, gVertexList[cornerPairs[gsEdgeList[sEdgeIndex].connect[1]->index]].xyzCoords[0], gVertexList[cornerPairs[gsEdgeList[sEdgeIndex].connect[1]->index]].xyzCoords[1], gVertexList[cornerPairs[gsEdgeList[sEdgeIndex].connect[1]->index]].xyzCoords[2], &pt11.pu, &pt11.pv, &err);
			assert(!err);
		
			iGeom_getEntXYZtoU(geom, gsEdgeList[sEdgeIndex].gEdgeHandle, gsEdgeList[sEdgeIndex].connect[0]->xyzCoords[0], gsEdgeList[sEdgeIndex].connect[0]->xyzCoords[1], gsEdgeList[sEdgeIndex].connect[0]->xyzCoords[2], &suLeft, &err);
			assert(!err);
			iGeom_getEntXYZtoU(geom, gsEdgeList[sEdgeIndex].gEdgeHandle, gsEdgeList[sEdgeIndex].connect[1]->xyzCoords[0], gsEdgeList[sEdgeIndex].connect[1]->xyzCoords[1], gsEdgeList[sEdgeIndex].connect[1]->xyzCoords[2], &suRight, &err);
			assert(!err);
		
			//find the corresponding the bottom mesh corner and top mesh corner
			iBase_EntitySetHandle tmpSet;
			SimpleArray<iBase_EntityHandle> tmpNodes;
			iRel_getEntSetRelation(assoc, rel, gsEdgeList[sEdgeIndex].connect[0]->gVertexHandle, 0, &tmpSet, &err);
			assert(!err);
			iMesh_getEntities(mesh, tmpSet, iBase_VERTEX, iMesh_POINT, ARRAY_INOUT(tmpNodes), &err);
			assert(!err);
			assert(tmpNodes.size()==1);
			iMesh_getIntData(mesh, tmpNodes[0], taghandle, &leftIndex, &err);
			assert(!err);
		
			//get the index for mesh vertices on the source edge
			iRel_getEntSetRelation(assoc, rel, gsEdgeList[sEdgeIndex].connect[1]->gVertexHandle, 0, &tmpSet, &err);
			assert(!err);
			iMesh_getEntities(mesh, tmpSet, iBase_VERTEX, iMesh_POINT, ARRAY_INOUT(tmpNodes), &err);
			assert(!err);
			assert(tmpNodes.size()==1);
			iMesh_getIntData(mesh, tmpNodes[0], taghandle, &RightIndex, &err);
			assert(!err);
		
			vector<int> sIndex(mNodes.size());  //sIndex is used to store the index for mesh nodes on the source surface boundary
			for (int j=0; j < mNodes.size(); j++)  //inner node on the edge
			{
				int VertexId;
				double u;//parametric coordinates on the source boundary edge
			
				Point2D ptuv[2];
				iMesh_getIntData(mesh, mNodes[j], taghandle, &VertexId, &err);
				assert(!err);
			
				sIndex[j] = VertexId;
			
				
			
				//now the geometrical bottom mesh corner is NodeList[leftIndex] and top mesh corner is TVertexList[RightIndex] 
				iGeom_getEntXYZtoUV(geom, gLinkFaceList[i].gFaceHandle, NodeList[sIndex[j]].xyzCoords[0], NodeList[sIndex[j]].xyzCoords[1], NodeList[sIndex[j]].xyzCoords[2], &ptuv[0].pu, &ptuv[0].pv, &err);
				assert(!err);
		
				iGeom_getEntXYZtoUV(geom, gLinkFaceList[i].gFaceHandle, TVertexList[sIndex[j]].xyzCoords[0], TVertexList[sIndex[j]].xyzCoords[1], TVertexList[sIndex[j]].xyzCoords[2], &ptuv[1].pu, &ptuv[1].pv, &err);
				assert(!err);
		
				//get the parametric coordinates mNodes[j] on the source boundary edge gsEdgeList[sEdgeIndex]
				iGeom_getEntXYZtoU(geom, gsEdgeList[sEdgeIndex].gEdgeHandle, NodeList[sIndex[j]].xyzCoords[0], NodeList[sIndex[j]].xyzCoords[1], NodeList[sIndex[j]].xyzCoords[2], &u, &err);
				assert(!err);
		
				//calculate the inner node in the linking surface
				for (int k=0; k < (numLayers-1); k++)
				{
					Point2D pt_0s, pt_1s, pt_r0, pt_r1;
					double r, s;
					s=(k+1)*1/double(numLayers);
			
					iGeom_getEntXYZtoUV(geom, gLinkFaceList[i].gFaceHandle, linkVertexList[k][leftIndex].xyzCoords[0], linkVertexList[k][leftIndex].xyzCoords[1], linkVertexList[k][leftIndex].xyzCoords[2], &pt_0s.pu, &pt_0s.pv, &err);
					assert(!err);
					iGeom_getEntXYZtoUV(geom, gLinkFaceList[i].gFaceHandle, linkVertexList[k][RightIndex].xyzCoords[0], linkVertexList[k][RightIndex].xyzCoords[1], linkVertexList[k][RightIndex].xyzCoords[2], &pt_1s.pu, &pt_1s.pv, &err);
					assert(!err);
			
					r = (u - suLeft)/(suRight - suLeft);
			
					pt_r0.pu = ptuv[0].pu;
					pt_r0.pv = ptuv[0].pv;
					pt_r1.pu = ptuv[1].pu;
					pt_r1.pv = ptuv[1].pv;
			
					linkVertexList[k][sIndex[j]].uvCoords[0] = linear_interpolation(s, pt_r0.pu, pt_r1.pu);
					linkVertexList[k][sIndex[j]].uvCoords[1] = linear_interpolation(r, pt_0s.pv, pt_1s.pv);
			
					//get the Physical coordinates for inner nodes on the linking surface
					iGeom_getEntUVtoXYZ(geom, gLinkFaceList[i].gFaceHandle, linkVertexList[k][sIndex[j]].uvCoords[0], linkVertexList[k][sIndex[j]].uvCoords[1], &linkVertexList[k][sIndex[j]].xyzCoords[0], &linkVertexList[k][sIndex[j]].xyzCoords[1], &linkVertexList[k][sIndex[j]].xyzCoords[2], &err);
					assert(!err);
			
					//create the vertex in the mesh
					nodehandleIndex++;
					nodeHandle.resize(nodehandleIndex);
					iMesh_createVtx(mesh, linkVertexList[k][sIndex[j]].xyzCoords[0], linkVertexList[k][sIndex[j]].xyzCoords[1], linkVertexList[k][sIndex[j]].xyzCoords[2], &nodeHandle[nodehandleIndex-1], &err);
					assert(!err);
					//add new generated vertex to the list
					linkVertexList[k][sIndex[j]].gVertexHandle = nodeHandle[nodehandleIndex-1];
			
					//create the line segments on the linking surface in the vertical direction
					int status;
					vector<iBase_EntityHandle> connect(2);
					if (k==0)
					{
						LineIndex++;
						lineHandle.resize(LineIndex);
						connect[0] = NodeList[sIndex[j]].gVertexHandle;
						connect[1] = nodeHandle[nodehandleIndex-1];
						iMesh_createEnt(mesh, iMesh_LINE_SEGMENT, &connect[0], 2, &lineHandle[LineIndex-1], &status, &err);
						assert(!err);
				
						if (k==(numLayers-2))
						{
							LineIndex++;
							lineHandle.resize(LineIndex);
							connect[0] = nodeHandle[nodehandleIndex-1];
							connect[1] = TVertexList[sIndex[j]].gVertexHandle;
							iMesh_createEnt(mesh, iMesh_LINE_SEGMENT, &connect[0], 2, &lineHandle[LineIndex-1], &status, &err);
							assert(!err);
						}
				
					}
					else
					{
						LineIndex++;
						lineHandle.resize(LineIndex);
						connect[0] = nodeHandle[nodehandleIndex-2];
						connect[1] = nodeHandle[nodehandleIndex-1];
						iMesh_createEnt(mesh, iMesh_LINE_SEGMENT, &connect[0], 2, &lineHandle[LineIndex-1], &status, &err);
						assert(!err);
						if (k==(numLayers-2))
						{
							LineIndex++;
							lineHandle.resize(LineIndex);
							connect[0] = nodeHandle[nodehandleIndex-1];
							connect[1] = TVertexList[sIndex[j]].gVertexHandle;
							iMesh_createEnt(mesh, iMesh_LINE_SEGMENT, &connect[0], 2, &lineHandle[LineIndex-1], &status, &err);
							assert(!err);
						}					
					}
					//create the inner line segments on the linking surface in the horizontal direction
					if (j==0)
					{
						LineIndex++;
						lineHandle.resize(LineIndex);
						connect[0] = linkVertexList[k][leftIndex].gVertexHandle;
						connect[1] = nodeHandle[nodehandleIndex-1];
						iMesh_createEnt(mesh, iMesh_LINE_SEGMENT, &connect[0], 2, &lineHandle[LineIndex-1], &status, &err);
						assert(!err);
						if (j == (mNodes.size()-1))
						{
							LineIndex++;
							lineHandle.resize(LineIndex);
							connect[0] = nodeHandle[nodehandleIndex-1];
							connect[1] = linkVertexList[k][RightIndex].gVertexHandle;
							iMesh_createEnt(mesh, iMesh_LINE_SEGMENT, &connect[0], 2, &lineHandle[LineIndex-1], &status, &err);
							assert(!err);
						}
					}
					else
					{
						//create the line segments between the adjacent vertical 
						LineIndex++;
						lineHandle.resize(LineIndex);
						connect[0] = linkVertexList[k][sIndex[j-1]].gVertexHandle;
						connect[1] = linkVertexList[k][sIndex[j]].gVertexHandle;
						iMesh_createEnt(mesh, iMesh_LINE_SEGMENT, &connect[0], 2, &lineHandle[LineIndex-1], &status, &err);
						assert(!err);
						if (j == (mNodes.size()-1))
						{
							LineIndex++;
							lineHandle.resize(LineIndex);
							connect[0] = nodeHandle[nodehandleIndex-1];
							connect[1] = linkVertexList[k][RightIndex].gVertexHandle;
							iMesh_createEnt(mesh, iMesh_LINE_SEGMENT, &connect[0], 2, &lineHandle[LineIndex-1], &status, &err);
							assert(!err);
						}
					}
			
			
					//create the face elements for the linking surface
					connect.resize(4);
					int sense_FaceRegion;
					iGeom_getEntNrmlSense(geom, gLinkFaceList[i].gFaceHandle, volEntity, &sense_FaceRegion, &err);
					assert(!err);
					int sense_EdgeFace;
					iGeom_getEgFcSense(geom, gsEdgeList[sEdgeIndex].gEdgeHandle, gLinkFaceList[i].gFaceHandle, &sense_EdgeFace, &err);
					assert(!err);
					int sense_vtxEdge;
					iGeom_getEgVtxSense(geom, gsEdgeList[sEdgeIndex].gEdgeHandle, gsEdgeList[sEdgeIndex].connect[0]->gVertexHandle, gsEdgeList[sEdgeIndex].connect[1]->gVertexHandle, &sense_vtxEdge, &err);
					assert(!err);
							
			
					if ((j==0))
					{
						if (k==0)
						{
							faceIndex++;
							faceHandle.resize(faceIndex);
							if (sense_EdgeFace*sense_vtxEdge > 0)
							{					
								connect[0] = NodeList[leftIndex].gVertexHandle;
								connect[1] = mNodes[j];
								connect[2] = linkVertexList[k][sIndex[j]].gVertexHandle;
								connect[3] = linkVertexList[k][leftIndex].gVertexHandle;
							}
							else
							{
								connect[0] = NodeList[leftIndex].gVertexHandle;
								connect[1] = linkVertexList[k][leftIndex].gVertexHandle;
								connect[2] = linkVertexList[k][sIndex[j]].gVertexHandle;
								connect[3] = mNodes[j];
							}
							iMesh_createEnt(mesh, iMesh_QUADRILATERAL, &connect[0], 4, &faceHandle[faceIndex-1], &status, &err);
							assert(!err);
					
					
					
					
							if (k==(numLayers-2))
							{
								faceIndex++;
								faceHandle.resize(faceIndex);
								if (sense_EdgeFace*sense_vtxEdge > 0)
								{	
									connect[0] = linkVertexList[k][leftIndex].gVertexHandle;
									connect[1] = linkVertexList[k][sIndex[j]].gVertexHandle;			
									connect[2] = TVertexList[sIndex[j]].gVertexHandle;
									connect[3] = TVertexList[leftIndex].gVertexHandle;
								}
								else
								{
									connect[0] = linkVertexList[k][leftIndex].gVertexHandle;
									connect[1] = TVertexList[leftIndex].gVertexHandle;
									connect[2] = TVertexList[sIndex[j]].gVertexHandle;
									connect[3] = linkVertexList[k][sIndex[j]].gVertexHandle;
								}
								iMesh_createEnt(mesh, iMesh_QUADRILATERAL, &connect[0], 4, &faceHandle[faceIndex-1], &status, &err);
								assert(!err);
							}
						}
						else
						{
							faceIndex++;
							faceHandle.resize(faceIndex);
							if (sense_EdgeFace*sense_vtxEdge > 0)
							{
								connect[0] = linkVertexList[k-1][leftIndex].gVertexHandle;
								connect[1] = linkVertexList[k-1][sIndex[j]].gVertexHandle;
								connect[2] = linkVertexList[k][sIndex[j]].gVertexHandle;
								connect[3] = linkVertexList[k][leftIndex].gVertexHandle;
							}
							else
							{
								connect[0] = linkVertexList[k-1][leftIndex].gVertexHandle;
								connect[1] = linkVertexList[k][leftIndex].gVertexHandle;				
								connect[2] = linkVertexList[k][sIndex[j]].gVertexHandle;
								connect[3] = linkVertexList[k-1][sIndex[j]].gVertexHandle;
							}
							iMesh_createEnt(mesh, iMesh_QUADRILATERAL, &connect[0], 4, &faceHandle[faceIndex-1], &status, &err);
							assert(!err);
							if (k==(numLayers-2))
							{
								faceIndex++;
								faceHandle.resize(faceIndex);
								if (sense_EdgeFace*sense_vtxEdge > 0)
								{
									connect[0] = linkVertexList[k][leftIndex].gVertexHandle;
									connect[1] = linkVertexList[k][sIndex[j]].gVertexHandle;	
									connect[2] = TVertexList[sIndex[j]].gVertexHandle;
									connect[3] = TVertexList[leftIndex].gVertexHandle;
								}
								else
								{
									connect[0] = linkVertexList[k][leftIndex].gVertexHandle;
									connect[1] = TVertexList[leftIndex].gVertexHandle;			
									connect[2] = TVertexList[sIndex[j]].gVertexHandle;
									connect[3] = linkVertexList[k][sIndex[j]].gVertexHandle;
								}
								iMesh_createEnt(mesh, iMesh_QUADRILATERAL, &connect[0], 4, &faceHandle[faceIndex-1], &status, &err);
								assert(!err);
							}
						}
						if (j== (mNodes.size()-1))
						{
							if (k==0)
							{
								faceIndex++;
								faceHandle.resize(faceIndex);
								if (sense_EdgeFace*sense_vtxEdge > 0)
								{
									connect[0] = NodeList[sIndex[j]].gVertexHandle;
									connect[1] = NodeList[RightIndex].gVertexHandle;	
									connect[2] = linkVertexList[k][RightIndex].gVertexHandle;
									connect[3] = linkVertexList[k][sIndex[j]].gVertexHandle;
								}
								else
								{
									connect[0] = NodeList[sIndex[j]].gVertexHandle;
									connect[1] = linkVertexList[k][sIndex[j]].gVertexHandle;		
									connect[2] = linkVertexList[k][RightIndex].gVertexHandle;
									connect[3] = NodeList[RightIndex].gVertexHandle;
								}
					
								iMesh_createEnt(mesh, iMesh_QUADRILATERAL, &connect[0], 4, &faceHandle[faceIndex-1], &status, &err);
								assert(!err);
						
								if (k==(numLayers-2))
								{
									faceIndex++;
									faceHandle.resize(faceIndex);
									if (sense_EdgeFace*sense_vtxEdge > 0)
									{
										connect[0] = linkVertexList[k][sIndex[j]].gVertexHandle;
										connect[1] = linkVertexList[k][RightIndex].gVertexHandle;
						
										connect[2] = TVertexList[RightIndex].gVertexHandle;
										connect[3] = TVertexList[sIndex[j]].gVertexHandle;
									}
									else
									{
										connect[0] = linkVertexList[k][sIndex[j]].gVertexHandle;
										connect[1] = TVertexList[sIndex[j]].gVertexHandle;
										connect[2] = TVertexList[RightIndex].gVertexHandle;
										connect[3] = linkVertexList[k][RightIndex].gVertexHandle;
									}
									iMesh_createEnt(mesh, iMesh_QUADRILATERAL, &connect[0], 4, &faceHandle[faceIndex-1], &status, &err);
									assert(!err);
								}
							}
							else
							{
								faceIndex++;
								faceHandle.resize(faceIndex);
								if (sense_EdgeFace*sense_vtxEdge > 0)
								{
									connect[0] = linkVertexList[k-1][sIndex[j]].gVertexHandle;
									connect[1] = linkVertexList[k-1][RightIndex].gVertexHandle;
									connect[2] = linkVertexList[k][RightIndex].gVertexHandle;
									connect[3] = linkVertexList[k][sIndex[j]].gVertexHandle;
								}
								else
								{
									connect[0] = linkVertexList[k-1][sIndex[j]].gVertexHandle;
									connect[1] = linkVertexList[k][sIndex[j]].gVertexHandle;			
									connect[2] = linkVertexList[k][RightIndex].gVertexHandle;
									connect[3] = linkVertexList[k-1][RightIndex].gVertexHandle;
								}
						
								iMesh_createEnt(mesh, iMesh_QUADRILATERAL, &connect[0], 4, &faceHandle[faceIndex-1], &status, &err);
								assert(!err);
								if (k==(numLayers-2))
								{
									faceIndex++;
									faceHandle.resize(faceIndex);
									if (sense_EdgeFace*sense_vtxEdge > 0)
									{
										connect[0] = linkVertexList[k][sIndex[j]].gVertexHandle;
										connect[1] = linkVertexList[k][RightIndex].gVertexHandle;
						
										connect[2] = TVertexList[RightIndex].gVertexHandle;
										connect[3] = TVertexList[sIndex[j]].gVertexHandle;
									}
									else
									{
										connect[0] = linkVertexList[k][sIndex[j]].gVertexHandle;
										connect[1] = TVertexList[sIndex[j]].gVertexHandle;			
										connect[2] = TVertexList[RightIndex].gVertexHandle;
										connect[3] = linkVertexList[k][RightIndex].gVertexHandle;
									}
									iMesh_createEnt(mesh, iMesh_QUADRILATERAL, &connect[0], 4, &faceHandle[faceIndex-1], &status, &err);
									assert(!err);
								}
					
							}
					
						}	
					}
					else
					{
						if (k==0)
						{
							faceIndex++;
							faceHandle.resize(faceIndex);
							if (sense_EdgeFace*sense_vtxEdge > 0)
							{
								connect[0] = NodeList[sIndex[j-1]].gVertexHandle;
								connect[1] = NodeList[sIndex[j]].gVertexHandle;						
								connect[2] = linkVertexList[k][sIndex[j]].gVertexHandle;
								connect[3] = linkVertexList[k][sIndex[j-1]].gVertexHandle;
							}
							else
							{
								connect[0] = NodeList[sIndex[j-1]].gVertexHandle;
								connect[1] = linkVertexList[k][sIndex[j-1]].gVertexHandle;
												
								connect[2] = linkVertexList[k][sIndex[j]].gVertexHandle;
								connect[3] = NodeList[sIndex[j]].gVertexHandle;
							}
					
							iMesh_createEnt(mesh, iMesh_QUADRILATERAL, &connect[0], 4, &faceHandle[faceIndex-1], &status, &err);
							assert(!err);
					
							if (k==(numLayers-2))
							{
								faceIndex++;
								faceHandle.resize(faceIndex);
								if (sense_EdgeFace*sense_vtxEdge > 0)
								{
									connect[0] = linkVertexList[k][sIndex[j-1]].gVertexHandle;
									connect[1] = linkVertexList[k][sIndex[j]].gVertexHandle;
						
									connect[2] = TVertexList[sIndex[j]].gVertexHandle;
									connect[3] = TVertexList[sIndex[j-1]].gVertexHandle;
								}
								else
								{
									connect[0] = linkVertexList[k][sIndex[j-1]].gVertexHandle;
									connect[1] = TVertexList[sIndex[j-1]].gVertexHandle;				
									connect[2] = TVertexList[sIndex[j]].gVertexHandle;
									connect[3] = linkVertexList[k][sIndex[j]].gVertexHandle;
								}
								iMesh_createEnt(mesh, iMesh_QUADRILATERAL, &connect[0], 4, &faceHandle[faceIndex-1], &status, &err);
								assert(!err);
							}
						}
						else
						{
							faceIndex++;
							faceHandle.resize(faceIndex);
							if (sense_EdgeFace*sense_vtxEdge > 0)
							{
								connect[0] = linkVertexList[k-1][sIndex[j-1]].gVertexHandle;
								connect[1] = linkVertexList[k-1][sIndex[j]].gVertexHandle;
								connect[2] = linkVertexList[k][sIndex[j]].gVertexHandle;
								connect[3] = linkVertexList[k][sIndex[j-1]].gVertexHandle;
							}
							else
							{
								connect[0] = linkVertexList[k-1][sIndex[j-1]].gVertexHandle;
								connect[1] = linkVertexList[k][sIndex[j-1]].gVertexHandle;			
								connect[2] = linkVertexList[k][sIndex[j]].gVertexHandle;
								connect[3] = linkVertexList[k-1][sIndex[j]].gVertexHandle;
							}
					
							iMesh_createEnt(mesh, iMesh_QUADRILATERAL, &connect[0], 4, &faceHandle[faceIndex-1], &status, &err);
							assert(!err);
							if (k==(numLayers-2))
							{
								faceIndex++;
								faceHandle.resize(faceIndex);
								if (sense_EdgeFace*sense_vtxEdge > 0)
								{
									connect[0] = linkVertexList[k][sIndex[j-1]].gVertexHandle;
									connect[1] = linkVertexList[k][sIndex[j]].gVertexHandle;	
									connect[2] = TVertexList[sIndex[j]].gVertexHandle;
									connect[3] = TVertexList[sIndex[j-1]].gVertexHandle;
								}
								else
								{
									connect[0] = linkVertexList[k][sIndex[j-1]].gVertexHandle;
									connect[1] = TVertexList[sIndex[j-1]].gVertexHandle;	
									connect[2] = TVertexList[sIndex[j]].gVertexHandle;
									connect[3] = linkVertexList[k][sIndex[j]].gVertexHandle;
								}
								iMesh_createEnt(mesh, iMesh_QUADRILATERAL, &connect[0], 4, &faceHandle[faceIndex-1], &status, &err);
								assert(!err);
							}
					
						}
						if (j== (mNodes.size()-1))
						{
							if (k==0)
							{
								faceIndex++;
								faceHandle.resize(faceIndex);
								if (sense_EdgeFace*sense_vtxEdge > 0)
								{
									connect[0] = NodeList[sIndex[j]].gVertexHandle;
									connect[1] = NodeList[RightIndex].gVertexHandle;	
									connect[2] = linkVertexList[k][RightIndex].gVertexHandle;
									connect[3] = linkVertexList[k][sIndex[j]].gVertexHandle;
								}
								else
								{
									connect[0] = NodeList[sIndex[j]].gVertexHandle;
									connect[1] = linkVertexList[k][sIndex[j]].gVertexHandle;
								
									connect[2] = linkVertexList[k][RightIndex].gVertexHandle;
									connect[3] = NodeList[RightIndex].gVertexHandle;
								}
					
								iMesh_createEnt(mesh, iMesh_QUADRILATERAL, &connect[0], 4, &faceHandle[faceIndex-1], &status, &err);
								assert(!err);
						
								if (k==(numLayers-2))
								{
									faceIndex++;
									faceHandle.resize(faceIndex);
									if (sense_EdgeFace*sense_vtxEdge > 0)
									{
										connect[0] = linkVertexList[k][sIndex[j]].gVertexHandle;
										connect[1] = linkVertexList[k][RightIndex].gVertexHandle;
						
										connect[2] = TVertexList[RightIndex].gVertexHandle;
										connect[3] = TVertexList[sIndex[j]].gVertexHandle;
									}
									else
									{
										connect[0] = linkVertexList[k][sIndex[j]].gVertexHandle;
										connect[1] = TVertexList[sIndex[j]].gVertexHandle;		
										connect[2] = TVertexList[RightIndex].gVertexHandle;
										connect[3] = linkVertexList[k][RightIndex].gVertexHandle;
									}
									iMesh_createEnt(mesh, iMesh_QUADRILATERAL, &connect[0], 4, &faceHandle[faceIndex-1], &status, &err);
									assert(!err);
								}
							}
							else
							{
								faceIndex++;
								faceHandle.resize(faceIndex);
								if (sense_EdgeFace*sense_vtxEdge > 0)
								{
									connect[0] = linkVertexList[k-1][sIndex[j]].gVertexHandle;
									connect[1] = linkVertexList[k-1][RightIndex].gVertexHandle;
									connect[2] = linkVertexList[k][RightIndex].gVertexHandle;
									connect[3] = linkVertexList[k][sIndex[j]].gVertexHandle;
								}
								else
								{
									connect[0] = linkVertexList[k-1][sIndex[j]].gVertexHandle;
									connect[1] = linkVertexList[k][sIndex[j]].gVertexHandle;
							
									connect[2] = linkVertexList[k][RightIndex].gVertexHandle;
									connect[3] = linkVertexList[k-1][RightIndex].gVertexHandle;
								}
					
								iMesh_createEnt(mesh, iMesh_QUADRILATERAL, &connect[0], 4, &faceHandle[faceIndex-1], &status, &err);
								assert(!err);
								if (k==(numLayers-2))
								{
									faceIndex++;
									faceHandle.resize(faceIndex);
									if (sense_EdgeFace*sense_vtxEdge > 0)
									{
										connect[0] = linkVertexList[k][sIndex[j]].gVertexHandle;
										connect[1] = linkVertexList[k][RightIndex].gVertexHandle;
						
										connect[2] = TVertexList[RightIndex].gVertexHandle;
										connect[3] = TVertexList[sIndex[j]].gVertexHandle;
									}
									else
									{
										connect[0] = linkVertexList[k][sIndex[j]].gVertexHandle;
										connect[1] = TVertexList[sIndex[j]].gVertexHandle;			
										connect[2] = TVertexList[RightIndex].gVertexHandle;
										connect[3] = linkVertexList[k][RightIndex].gVertexHandle;
									}
									iMesh_createEnt(mesh, iMesh_QUADRILATERAL, &connect[0], 4, &faceHandle[faceIndex-1], &status, &err);
									assert(!err);
								}					
							}					
						}
					}			
				}
			}
	
			//create the entityset for linking surface[i]
			iBase_EntitySetHandle entityset;
			iRel_getEntSetRelation(assoc, rel, gLinkFaceList[i].gFaceHandle, 0, &entityset, &err);
			if (err)	//there is no entityset associated with gLinkFaceList[i].gFaceHandle
			{
				iMesh_createEntSet(mesh, 1, &entityset, &err);
				assert(!err);
			}
			iMesh_addEntArrToSet(mesh, &nodeHandle[0], nodehandleIndex, entityset, &err);
			assert(!err);
			iMesh_addEntArrToSet(mesh, &lineHandle[0], LineIndex, entityset, &err);
			assert(!err);
			iMesh_addEntArrToSet(mesh, &faceHandle[0], faceIndex, entityset, &err);
			assert(!err);
	
			iRel_getEntSetRelation(assoc, rel, gLinkFaceList[i].gFaceHandle, 0, &entityset, &err);
			if (err)//there is no entityset associated with gLinkFaceList[i].gFaceHandle
			{
				//create the irel between the linking sides and entityset
				iRel_setEntSetRelation(assoc, rel, gLinkFaceList[i].gFaceHandle, entityset, &err);
				assert(!err);
			}
		}
		
	}
	
	return 1;
	
}

double OneToOneSwept::parametricTFI2D(double r, double s, double pt_0s, double pt_1s, double pt_r0, double pt_r1)
{
	assert(r>= 0 && r <= 1.0);
	assert(s>= 0 && s <= 1.0);
	double pt_rs;

	//interpolate the pt_rs based on pt_r0, pt_r1, pt_0s and pt_1s
	pt_rs = 0.5*((1-s)*pt_r0 + s*pt_r1 + (1-r)*pt_0s + r*pt_1s);

	return pt_rs;
}

double OneToOneSwept::linear_interpolation(double r, double x0, double x1)
{
	assert(r >=0 && r <= 1.0);

	double pt= (1-r)*x0 + r*x1;

	return pt;
}

