#include "EdgeMesher.hpp"
#include <iostream>
#include <math.h>

EdgeMesher::EdgeMesher(iGeom_Instance &geom, iBase_EntityHandle *EdgeHandle, iMesh_Instance &Mesh, iRel_Instance &association, iRel_PairHandle *irel)
{
	geometry = geom;
	gEdgeHandle = *EdgeHandle;

	int err;
	iGeom_getEntURange(geometry, gEdgeHandle, &umin, &umax, &err);
	assert(!err);

	int in_u, in_v;
	iGeom_isEntPeriodic(geometry, gEdgeHandle, &in_u, &in_v, &err);
	assert(!err);

	periodic = in_u;
	    
	double x0,y0,z0,x1,y1,z1;
	iGeom_getEntUtoXYZ(geometry, gEdgeHandle, umin, &x0, &y0, &z0, &err);
	iGeom_getEntUtoXYZ(geometry, gEdgeHandle, umax, &x1, &y1, &z1, &err);

	double dx = x0 - x1;
	double dy = y0 - y1;
	double dz = z0 - z1;
	double eps = 1.0E-08;
	if (dx * dx + dy * dy + dz * dz < eps * eps) periodic = 1;    

	mesh  = Mesh;
	assoc = association;
	rel   = *irel;
	
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
EdgeMesher::~EdgeMesher()
{
	std::cout << "it is over now" << std::endl;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double EdgeMesher::getLength() const
{
	int err;
	SimpleArray<double> measure;
	iGeom_measure(geometry, &gEdgeHandle, 1, ARRAY_INOUT(measure), &err);
    	return measure[0];
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double EdgeMesher::getLength(double ustart, double uend) const
{
	int err;
	double arclen;
    	
	SimpleArray<double> measure;
	iGeom_measure(geometry, &gEdgeHandle, 1, ARRAY_INOUT(measure), &err);
	arclen=measure[0]*(uend-ustart)/(umax-umin);
	
    	return arclen;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void EdgeMesher::SetScheme(int SetOption)
{
	switch (SetOption)
	{
		case 0:
			SchemeOption=equalMesh;
			SelectedShemeName="Equal Meshing";
			break;
		case 1:
			SchemeOption=biasMesh;
			SelectedShemeName="Bias Meshing";
			break;
		case 2:
			SchemeOption=dualMesh;
			SelectedShemeName="Dual Meshing";
			break;
		case 3:
			SchemeOption=curvatureMesh;
			SelectedShemeName="Curvature Meshing";
			break;
		default:
			std::cout << "There is an error with scheme option" << std::endl;
			break;
	}		
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
const char *EdgeMesher::getSchemeName()
{
	return SelectedShemeName;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double EdgeMesher::getStepSize()
{
	return stepSize;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void EdgeMesher::SetStepSize(double StepSize)
{
	double len=getLength();
	assert(StepSize);	
	int num=int(len/StepSize+0.5);
	NumEdges=num;
	stepSize=StepSize;	
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int EdgeMesher::Execute()
{
	int Result=0;
	if ((SchemeOption < 4) && (SchemeOption > -1))
	{
		Result=1;		
		EdgeMesh();
		
	}	
	return Result;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void EdgeMesher::EdgeMesh()
{
	int err, NumEntitySet;
	SimpleArray<iBase_EntityHandle> nodeHandles;
	SimpleArray<iBase_EntityHandle> lineHandles;
	SimpleArray<iBase_EntityHandle> edgeHandles;
	SimpleArray<int> status;
	vector<double> NodeCoordinates;
			
	assert(NumEdges);//make sure NumNodes is non-zero integer
	NodeCoordinates.resize( 3*(NumEdges+1));
		
	SimpleArray<iBase_EntityHandle> gNodes;
	iGeom_getEntAdj(geometry, gEdgeHandle, iBase_VERTEX, ARRAY_INOUT(gNodes), &err);	
	
	switch(SchemeOption)
	{
		case equalMesh:
			NodeCoordinates = EqualMeshing();
			break;			
		case biasMesh:
			NodeCoordinates = BiasMeshing();
			break;
		case dualMesh:
			NodeCoordinates = DualBiasMeshing();
			break;
		case curvatureMesh:
			NodeCoordinates = CurvatureMeshing();
			break;
		default:
			std::cout << "The scheme option is set incorrectly!" << std::endl;
			break;
	}
		
	iMesh_createVtxArr(mesh, (NumEdges+1), iBase_INTERLEAVED, &NodeCoordinates[0], 3*(NumEdges+1), ARRAY_INOUT(nodeHandles), &err);
	assert(!err);
	
	edgeHandles.resize(2*NumEdges);
	edgeHandles[0]=nodeHandles[0];
	for (int k=1; k < (NumEdges); k++)
	{
		edgeHandles[2*k-1] = nodeHandles[k];
		edgeHandles[2*k] = nodeHandles[k];		
	}
	edgeHandles[2*NumEdges-1]=nodeHandles[NumEdges];
	
	iMesh_createEntArr(mesh, iMesh_LINE_SEGMENT, &edgeHandles[0], (2*NumEdges), ARRAY_INOUT(lineHandles), ARRAY_INOUT(status), &err);
	assert(!err);
	
	
	iBase_EntitySetHandle mesh_entityset;
	
	get_related_entityset(mesh_entityset);
	//set an array of entities to EntitySet	
	iMesh_addEntArrToSet(mesh, &nodeHandles[0], (NumEdges+1), mesh_entityset, &err);
	iMesh_addEntArrToSet(mesh, &lineHandles[0], (NumEdges), mesh_entityset, &err);
	assert(!err);
			
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
vector<double> EdgeMesher::EqualMeshing()
{
	int err;
	double x, y, z, u, du;
	vector<double> NodeCoordinates;  //temporarily store the coordinates

	assert(NumEdges);//make sure NumEdges is non-zero integer
	du = (umax-umin)/(double)NumEdges;	

	u=umin;
	NodeCoordinates.resize(3*(NumEdges+1));
	for(int i = 0; i < (NumEdges+1); i++)
	{
		u = umin + i*du;
		iGeom_getEntUtoXYZ(geometry, gEdgeHandle, u, &x, &y, &z, &err);
		NodeCoordinates[3*i] = x;
		NodeCoordinates[3*i+1] = y;
		NodeCoordinates[3*i+2] = z;	
	}
	
	return NodeCoordinates;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
vector<double> EdgeMesher::BiasMeshing()
{	
	int err, i;
	double x, y, z, ustart, u, du, len, q, L0, dist;
	vector<double> NodeCoordinates;  //temporarily store the coordinates
	Point3D tempCoordinate;	

	assert(NumEdges);
	NodeCoordinates.resize(3*(NumEdges+1));
	std::cout << "Please input the value for q " << std::endl; 	
	std::cin >> q;
	
	assert(q-1);//make sure q is not equal to 1
	assert(q>0);//make sure q is positive
	len = getLength();
	L0 = len*(1-q)/(1-pow(q,NumEdges));
	
	//discretizing the edge
	tempCoordinate=getXYZCoords(umin);
	NodeCoordinates[3*0] = tempCoordinate.px;
	NodeCoordinates[3*0+1] = tempCoordinate.py;
	NodeCoordinates[3*0+2] = tempCoordinate.pz;
	ustart = umin;
	u=ustart+(umax-umin)*L0/len;
	i=1;
	
	while (i < (NumEdges))
	{
		
		dist = L0*pow(q,i);
		u=getUCoord(ustart, dist, u);
		tempCoordinate=getXYZCoords(u);
		NodeCoordinates[3*i] = tempCoordinate.px;
		NodeCoordinates[3*i+1] = tempCoordinate.py;
		NodeCoordinates[3*i+2] = tempCoordinate.pz;
		ustart = u;
		u=ustart+(umax-umin)*dist/len;
		i++;
	}
	tempCoordinate = getXYZCoords(umax);
	NodeCoordinates[3*NumEdges] = tempCoordinate.px;
	NodeCoordinates[3*NumEdges+1] = tempCoordinate.py;
	NodeCoordinates[3*NumEdges+2] = tempCoordinate.pz;

	return NodeCoordinates;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
vector<double> EdgeMesher::DualBiasMeshing()
{	
	int err, i, NumNodes;
	double x, y, z, ustart, u, du, len, q, L0, dist, ustart0, ustart1;
	vector<double> NodeCoordinates;  //temporarily store the coordinates
	Point3D tempCoordinate;

	assert(NumEdges);
	if ((NumEdges%2)!=0)
	{
		NumNodes = NumEdges + 1;
	}	
	else
	{	
		NumNodes = NumEdges;
	}	
	NodeCoordinates.resize(3*(NumNodes+1));
	std::cout << "Please input the value for q " << std::endl; 	
	std::cin >> q;

	assert(q-1);//make sure q is not equal to 1
	assert(q>0);//make sure q is positive
		
	len = getLength();
	L0 = 0.5*len*(1-q)/(1-pow(q,NumNodes/2));
	
	//discretizing the edge
	tempCoordinate = getXYZCoords(umin);
	NodeCoordinates[0] = tempCoordinate.px;
	NodeCoordinates[1] = tempCoordinate.py;
	NodeCoordinates[2] = tempCoordinate.pz;

	ustart0 = umin;
	ustart1 = umax;
	tempCoordinate = getXYZCoords(umax);
	NodeCoordinates[3*NumNodes] = tempCoordinate.px;
	NodeCoordinates[3*NumNodes+1] = tempCoordinate.py;
	NodeCoordinates[3*NumNodes+2] = tempCoordinate.pz;
	i=1;	
	while (i < ((NumNodes/2)+1))
	{
		dist = L0*pow(q,i-1);
		u=ustart0+(umax-umin)*dist/len;
		u=getUCoord(ustart0, dist, u);
		tempCoordinate = getXYZCoords(u);
		NodeCoordinates[3*i] = tempCoordinate.px;
		NodeCoordinates[3*i+1] = tempCoordinate.py;
		NodeCoordinates[3*i+2] = tempCoordinate.pz;
		ustart0 = u;
		//NodeCoordinates[NumNodes-i] = NodeCoordinates[i];
		
		dist = L0*pow(q,i-1);
		u=ustart1-(umax-umin)*dist/len;
		u=getUCoord(ustart1, dist, u);
		tempCoordinate = getXYZCoords(u);
		NodeCoordinates[3*(NumNodes-i)] = tempCoordinate.px;
		NodeCoordinates[3*(NumNodes-i)+1] = tempCoordinate.py;
		NodeCoordinates[3*(NumNodes-i)+2] = tempCoordinate.pz;		
		
		ustart1 = u;	
		i++;
	}

	return NodeCoordinates;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
vector<double> EdgeMesher::CurvatureMeshing()
{
	int err,index=0;
	double x, y, z, u, du, uMid;
	//temporarily store the coordinates
	vector<double> NodeCoordinates;	
	vector<double> TempNode;
	vector<double> URecord;//record the value of U
		
	Point3D pts0, pts1, ptsMid;

	assert(NumEdges);//make sure NumEdges is a non-zero integer
	du = (umax-umin)/(double)NumEdges;	

	NodeCoordinates.resize(3*(NumEdges+1));

	TempNode.resize(3*1);
	URecord.resize(1);
	
	iGeom_getEntUtoXYZ(geometry, gEdgeHandle, umin, &x, &y, &z, &err);

	NodeCoordinates[3*0]=x;
	NodeCoordinates[3*0+1]=y;
	NodeCoordinates[3*0+2]=z;
	TempNode[3*0] = x;
	TempNode[3*0+1] = y;
	TempNode[3*0+2] = z;
	URecord[0] = umin;
		
	for(int i = 1; i < (NumEdges+1); i++)
	{
		u = umin + i*du;
		iGeom_getEntUtoXYZ(geometry, gEdgeHandle, u, &x, &y, &z, &err);
		
		NodeCoordinates[3*i] = x;
		NodeCoordinates[3*i+1] = y;
		NodeCoordinates[3*i+2] = z;
		
		pts0.px = NodeCoordinates[3*(i-1)];
		pts0.py = NodeCoordinates[3*(i-1)+1];
		pts0.pz = NodeCoordinates[3*(i-1)+2];

		pts1.px = NodeCoordinates[3*i];
		pts1.py = NodeCoordinates[3*i+1];
		pts1.pz = NodeCoordinates[3*i+2];
		uMid = (u-du+u)/2;
		iGeom_getEntUtoXYZ(geometry, gEdgeHandle, uMid, &ptsMid.px, &ptsMid.py, &ptsMid.pz, &err);
		if(!ErrorCalculate(pts0, pts1, ptsMid))
		{
			DivideIntoMore(pts0, ptsMid, pts1, u-du, u, uMid, index, TempNode, URecord);
		}
		// add the other end node to the array		
		{
			index++;
			TempNode.resize(3*(index+1));
			URecord.resize(index+1);
			TempNode[3*index]=pts1.px;
			TempNode[3*index+1]=pts1.py;
			TempNode[3*index+2]=pts1.pz;
			URecord[index]=u;
		}				
	}
	
	//sorting the coordinate data based on the value of u
	assert(TempNode.size()== (3*URecord.size()));
	
	QuickSorting(TempNode, URecord, URecord.size());

	//std::cout << "The maximal size of URecord and TempNode are " << TempNode.size() << "  " << URecord.size() <<std::endl;
	NumEdges=URecord.size()-1;
	return TempNode;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void EdgeMesher::DivideIntoMore(Point3D p0, Point3D pMid, Point3D p1, double u0, double u1, double uMid, int &index, vector<double> &nodes, vector<double> &URecord)
{
	//this is a recursive process
	double uu0, uu1, uumid;
	Point3D pts0, pts1, ptsMid;
	int err;
	
	index++;
	nodes.resize(3*(index+1));
	URecord.resize(index+1);
	nodes[3*index] = pMid.px;
	nodes[3*index+1] = pMid.py;
	nodes[3*index+2] = pMid.pz;
	URecord[index] = uMid;
	
	//left side
	uu0=u0;
	uu1=uMid;
	uumid=(uu0+uu1)/2;
	pts0=p0;
	pts1=pMid;
	iGeom_getEntUtoXYZ(geometry, gEdgeHandle, uumid, &ptsMid.px, &ptsMid.py, &ptsMid.pz, &err);
	if(!ErrorCalculate(pts0, pts1, ptsMid))
	{
		DivideIntoMore(pts0, ptsMid, pts1, uu0, uu1, uumid, index, nodes, URecord);
	}
	
	//right side
	uu0 = uMid;
	uu1=u1;
	uumid=(uu0+uu1)/2;
	pts0=pMid;
	pts1=p1;
	iGeom_getEntUtoXYZ(geometry, gEdgeHandle, uumid, &ptsMid.px, &ptsMid.py, &ptsMid.pz, &err);
	if(!ErrorCalculate(pts0, pts1, ptsMid))
	{
		DivideIntoMore(pts0, ptsMid, pts1, uu0, uu1, uumid, index, nodes, URecord);
	}
			
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//rapid sorting
void EdgeMesher::RapidSorting(vector<double> &nodes, vector<double> &URecord, int left, int right)
{
	int i, j;
	double middle, iTemp, x, y, z;
	Point3D TempData;
	
	middle=URecord[(left+right)/2];
	i=left;
	j=right;
	
	do
	{
		//search the values which are greater than the middle value from the left side		
		while((URecord[i] < middle)&&(i<right))
		{
			i++;
		}
		//search the values which are greater than the middle value from the right side
		while((URecord[j] > middle)&&(j > left))
		{
			j--;
		}
		if (i<=j)//find a pair of values
		{
			iTemp = URecord[i];
			URecord[i] = URecord[j];
			URecord[j]=iTemp;
			
			
			TempData.px = nodes[3*i];
			TempData.py = nodes[3*i+1];
			TempData.pz = nodes[3*i+2];

			nodes[3*i] = nodes[3*j];
			nodes[3*i+1] = nodes[3*j+1];
			nodes[3*i+2] = nodes[3*j+2];
			nodes[3*j] = TempData.px;
			nodes[3*j+1] = TempData.py;
			nodes[3*j+2] = TempData.pz;			
			
			i++;
			j--;
		}
	}while(i<=j);
	if (left < j)
		RapidSorting(nodes, URecord, left, j);
	if (right > i)
		RapidSorting(nodes, URecord, i, right);	
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void EdgeMesher::QuickSorting(vector<double> &nodes, vector<double> &URecord, int count)
{
	RapidSorting(nodes, URecord, 0, count-1);
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Point3D EdgeMesher::getXYZCoords(double u) const
{
	Point3D pts3D;
	double x, y, z;
	
	int err;
	iGeom_getEntUtoXYZ(geometry, gEdgeHandle, u, &x, &y, &z, &err);
	assert(!err);

	pts3D.px = x;
	pts3D.py = y;
	pts3D.pz = z;
	return pts3D;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//give a distance and starting point ustart, determine the next point
double EdgeMesher::getUCoord(double ustart, double dist, double uguess) const
{

	Point3D p0 = getXYZCoords(ustart);
	Point3D p1 = getXYZCoords(uguess);


	double dx, dy, dz, dl, u=uguess;
	double tol = 1.0E-7;
	int test=0;

	int ntrials=0;
	while(1)
	{
		dx = p1.px - p0.px;
        	dy = p1.py - p0.py;
        	dz = p1.pz - p0.pz;
        	dl = sqrt(dx * dx + dy * dy + dz * dz);
        	if ( fabs(dl-dist) < tol) break;
		
        	u = ustart + (u - ustart) * (dist/dl);
		if (u > umax)
		{
			u=umax;
			test++;
			if (test>10) break;
		}		
		if (u < umin)
		{		
			u=umin;
			test++;
			if (test>10) break;
		}        	
		p1 = getXYZCoords(u);
		

        if (ntrials++ == 100000)
        {
            cout << " Warning: Searching for U failed " << endl;
        }
	}
	uguess = u;
	return uguess;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool EdgeMesher::ErrorCalculate(Point3D p0, Point3D p1, Point3D pMid)
{
	double lengtha, lengthb, lengthc;
	double deltax, deltay, deltaz;
	double angle, error, tol=1.0E-3, H;
	double cvtr_i, cvtr_j, cvtr_k, curvature;
	bool result;
	int err;
	deltax = pMid.px-p0.px;
	deltay = pMid.py-p0.py;
	deltaz = pMid.pz-p0.pz;	
	lengtha = sqrt(deltax*deltax + deltay*deltay + deltaz*deltaz);

	deltax = p1.px-p0.px;
	deltay = p1.py-p0.py;
	deltaz = p1.pz-p0.pz;	
	lengthb = sqrt(deltax*deltax + deltay*deltay + deltaz*deltaz);

	deltax = pMid.px-p1.px;
	deltay = pMid.py-p1.py;
	deltaz = pMid.pz-p1.pz;	
	lengthc = sqrt(deltax*deltax + deltay*deltay + deltaz*deltaz);

	angle = acos((lengtha*lengtha + lengthb*lengthb - lengthc*lengthc)/(2*lengtha*lengthb));
	H = fabs(lengtha*sin(angle));

	iGeom_getEgCvtrXYZ(geometry, gEdgeHandle, pMid.px, pMid.py, pMid.pz, &cvtr_i, &cvtr_j, &cvtr_k, &err);
	curvature = sqrt(cvtr_i*cvtr_i+cvtr_j*cvtr_j+cvtr_k*cvtr_k);
	error= H*curvature;
	
	
	if (error > tol)
		result = false;
	else
		result = true;
	return result;		
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void EdgeMesher::get_related_entityset(iBase_EntitySetHandle &mesh_entityset)
{
	int switch_order=0;
	int ierr;

	iRel_getEntSetRelation(assoc, rel, gEdgeHandle, switch_order, &(mesh_entityset), &ierr);
	if (ierr)
	{
                iMesh_createEntSet(mesh, 1, &mesh_entityset, &ierr);
                assert(!ierr);

                iBase_TagHandle dim_tag;
                const char *tag2 = "GEOM_DIMENSION";
                int namelen = strlen(tag2);
                iMesh_getTagHandle(mesh, tag2, &dim_tag, &ierr, namelen);
                assert(!ierr);

                int dim = 1;
                iMesh_setEntSetIntData(mesh, mesh_entityset, dim_tag, dim, &ierr);
		assert(!ierr);

                iRel_setEntSetRelation(assoc, rel, gEdgeHandle, mesh_entityset, &ierr);
		assert(!ierr);
	}
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void EdgeMesher::ShowCoorData()
{
	int err;
	iBase_EntitySetHandle meshSet;
	SimpleArray<iBase_EntityHandle> VertexArr;
	SimpleArray<double> CoordinateArr;	
	
	iRel_getEntSetRelation(assoc, rel, gEdgeHandle, 0, &meshSet, &err);
	assert(!err);
	
	iMesh_getEntities(mesh, meshSet, iBase_VERTEX, iMesh_POINT, ARRAY_INOUT(VertexArr), &err);
	assert(!err);
	
	double x, y, z;
	for(int i=0;i<VertexArr.size();i++)
	{
		iMesh_getVtxCoord(mesh, VertexArr[i], &x, &y, &z, &err);		
		std::cout << "coordinates" << (i+1) << " are x: " << x << " y: " <<y << " z:" << z << std::endl;
	}	
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int EdgeMesher::SaveFile(const char *FileName)
{
	int namelen = strlen(FileName), err;
	iBase_EntitySetHandle mesh_root_set;
	iMesh_getRootSet(mesh, &mesh_root_set, &err);
    	assert(!err);
    	iMesh_save(mesh, mesh_root_set, FileName, NULL, &err, namelen, 0);
	if(err)
	{
		return 1;//fail
	}
	else
	{
		return 0;//succeed
	}
	
}

