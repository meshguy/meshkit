#include "meshkit/EdgeMesher.hpp"
#include "meshkit/MKCore.hpp"
#include "meshkit/ModelEnt.hpp"
#include "meshkit/SizingFunction.hpp"
#include "meshkit/RegisterMeshOp.hpp"
#include "moab/ReadUtilIface.hpp"
#include <vector>
#include <math.h>

namespace MeshKit
{

//---------------------------------------------------------------------------//
//Entity Type initilization for edge meshing
moab::EntityType EdgeMesher_tps[] = {moab::MBVERTEX, moab::MBEDGE, moab::MBMAXTYPE};
const moab::EntityType* EdgeMesher::output_types()
  { return EdgeMesher_tps; }

//---------------------------------------------------------------------------//
// Construction Function for Edge Mesher
EdgeMesher::EdgeMesher(MKCore *mk_core, const MEntVector &me_vec) 
        : MeshScheme(mk_core, me_vec), schemeType(EQUAL)
{

}
#if 0
//---------------------------------------------IBERRCHK(g_err, "Trouble get the adjacent geometric nodes on a surface.");------------------------------//
// measure function: compute the distance between the parametric coordinate 
// ustart and the parametric coordinate uend
double EdgeMesher::measure(iGeom::EntityHandle ent, double ustart, double uend) const
{
	double umin, umax;

	//get the minimal and maximal parametrical coordinates of the edge
	iGeom::Error err = mk_core()->igeom_instance()->getEntURange(ent, umin, umax);
	IBERRCHK(err, "Trouble getting parameter range for edge.");

	if (umin == umax) throw Error(MK_BAD_GEOMETRIC_EVALUATION, "Edge evaluated to same parameter umax and umin.");

	//compute the distance for edges
	double measure;
	err = mk_core()->igeom_instance()->measure(&ent, 1, &measure);

	IBERRCHK(err, "Trouble getting edge measure.");

	return measure * (uend - ustart) / (umax - umin);
}
#endif
//---------------------------------------------------------------------------//
// setup function: set up the number of intervals for edge meshing through the 
// sizing function
void EdgeMesher::setup_this()
{
    //compute the number of intervals for the associated ModelEnts, from the size set on them
    //the sizing function they point to, or a default sizing function
  for (MEntSelection::iterator mit = mentSelection.begin(); mit != mentSelection.end(); mit++)
  {
    ModelEnt *me = mit->first;

      //first check to see whether entity is meshed
    if (me->get_meshed_state() >= COMPLETE_MESH || me->mesh_intervals() > 0)
      continue;
    
    SizingFunction *sf = mk_core()->sizing_function(me->sizing_function_index());
    if (!sf && me -> mesh_intervals() < 0 && me -> interval_firmness() == DEFAULT &&
        mk_core()->sizing_function(0))
      sf = mk_core()->sizing_function(0);
    
    if (!sf && me -> mesh_intervals() < 0 && me -> interval_firmness() == DEFAULT)
    {
        //no sizing set, just assume default #intervals as 4
      me->mesh_intervals(4);
      me->interval_firmness(DEFAULT);
    }
    else
    {
        //check # intervals first, then size, and just choose for now
      if (sf->intervals() > 0)
      {
        if (me->constrain_even() && sf->intervals()%2)
          me -> mesh_intervals(sf->intervals()+1);
        else
          me -> mesh_intervals(sf->intervals());
        me -> interval_firmness(HARD);
      }
      else if (sf->size()>0)
      {
        int intervals = me->measure()/sf->size();
        if (!intervals) intervals++;
        if (me->constrain_even() && intervals%2) intervals++;
        me->mesh_intervals(intervals);
        me->interval_firmness(SOFT);
      }
      else
        throw Error(MK_INCOMPLETE_MESH_SPECIFICATION,  "Sizing function for edge had neither positive size nor positive intervals.");
    }
  }

    // now call setup_boundary to treat vertices
  setup_boundary();
}

//---------------------------------------------------------------------------//
// execute function: Generate the mesh for edges
void EdgeMesher::execute_this()
{
  std::vector<double> coords;
  std::vector<moab::EntityHandle> nodes;

  for (MEntSelection::iterator mit = mentSelection.begin(); mit != mentSelection.end(); mit++)
  {
    ModelEnt *me = mit -> first;
    if (me->get_meshed_state() >= COMPLETE_MESH)
	continue;
    //resize the coords based on the interval setting
    int num_edges = me->mesh_intervals();
    coords.resize(3*(num_edges+1));
    nodes.clear();
    nodes.reserve(num_edges + 1);

    //get bounding mesh entities, use 1st 2 entries of nodes list temporarily
    //pick up the boundary end nodes
    me->boundary(0, nodes);

    bool periodic = (nodes.size() == 1);
    
    //get coords in list, then move one tuple to the last position
    moab::ErrorCode rval = mk_core()->moab_instance()->get_coords(&nodes[0], nodes.size(), &coords[0]);
    MBERRCHK(rval, mk_core()->moab_instance());

    //move the second node to the endmost position in the node list
    // if periodic, the last node coordinates are also at index 0 in coords array
    // there is only one node, coords[3+i] are not even initialized
    int index2 = (periodic) ? 0 : 3;
    for (int i = 0; i < 3; i++)
      coords[3*num_edges+i] = coords[index2+i];


    //choose the scheme for edge mesher
    switch(schemeType)
    {
      case EQUAL://equal meshing for edges
          EqualMeshing(me, num_edges, coords);
          break;
      case BIAS://bias meshing for edges
          BiasMeshing(me, num_edges, coords);
          break;
      case DUAL://dual bias meshing for edges
          DualBiasMeshing(me, num_edges, coords);
          break;
      case CURVATURE://curvature-based meshing for edges
          CurvatureMeshing(me, num_edges, coords);
          break;
      default:
          break;			
    }
    //the variable nodes should be resized, node size may be changed in the different scheme for edge meshing
    me->mesh_intervals(num_edges);
    nodes.resize(num_edges+1);
    
    //move the other nodes to the end of nodes' list
    if (periodic) nodes[num_edges] = nodes[0];
    else nodes[num_edges] = nodes[1];

    //create the vertices' entities on the edge
    if (num_edges > 1) {
      rval = mk_core()->moab_instance()->create_vertices(&coords[3], num_edges - 1, mit -> second);
      MBERRCHK(rval, mk_core()->moab_instance());
    }
    
    //distribute the nodes into vector
    int j = 1;
    for (moab::Range::iterator rit = mit -> second.begin(); rit != mit -> second.end(); rit++)
      nodes[j++] = *rit;

      //get the query interface, which we will use to create the edges directly 
    moab::ReadUtilIface *iface;
    rval = mk_core() -> moab_instance() -> query_interface(iface);
    MBERRCHK(rval, mk_core()->moab_instance());		

      //create the edges, get a direct ptr to connectivity
    moab::EntityHandle starth, *connect, *tmp_connect;
    rval = iface -> get_element_connect(num_edges, 2, moab::MBEDGE, 1, starth, connect);
    MBERRCHK(rval, mk_core()->moab_instance());

      //add edges to range for the MESelection
    mit -> second.insert(starth, starth + num_edges - 1);

      //now set the connectivity array from the nodes
    tmp_connect = &nodes[0];
    for (int i = 0; i < num_edges; i++)
    {
      connect[0] = tmp_connect[0];
      connect[1] = tmp_connect[1];

        //increment connectivity ptr by 2 to go to next edge
      connect += 2;
			
        //increment tmp_connect by 1, to go to next node
      tmp_connect++;
    }

      //   ok, we are done, commit to ME
    me->commit_mesh(mit->second, COMPLETE_MESH);
   
	
  }

	
}

//---------------------------------------------------------------------------//
// Deconstruction function
EdgeMesher::~EdgeMesher()
{

}

//---------------------------------------------------------------------------//
// Create the mesh for edges with the equal distances
void EdgeMesher::EqualMeshing(ModelEnt *ent, int num_edges, std::vector<double> &coords)
{
  double umin, umax, measure;

    //get the u range for the edge
  iGeom::Error gerr =  ent->igeom_instance()->getEntURange(ent->geom_handle(), umin, umax);
  IBERRCHK(gerr, "Trouble get parameter range for edge.");

  if (umin == umax) throw Error(MK_BAD_GEOMETRIC_EVALUATION, "Edge evaluated to some parameter umax and umin.");

  //get the arc length
  measure = ent -> measure();

  double u, du;
  if (!num_edges) throw Error(MK_BAD_INPUT, "Trying to mesh edge with zero edges.");
  du = (umax - umin)/(double)num_edges;
	
  u = umin;
  //distribute the nodes with equal distances
  for (int i = 1; i < num_edges; i++)
  {
    u = umin + i*du;
    gerr = ent->igeom_instance()->getEntUtoXYZ(ent->geom_handle(), u, coords[3*i], coords[3*i+1], coords[3*i+2]);
    IBERRCHK(gerr, "Trouble getting U from XYZ along the edge.");
  }

}

//---------------------------------------------------------------------------//
// Create the mesh for edges based on the curvatures
void EdgeMesher::CurvatureMeshing(ModelEnt *ent, int &num_edges, std::vector<double> &coords)
{
  double umin, umax, measure;
  //store the initial edge size, the edge size may be changed during meshing
  int initial_num_edges = num_edges;

  //get the u range for the edge
  iGeom::Error gerr = ent->igeom_instance() ->getEntURange(ent->geom_handle(), umin, umax);
  IBERRCHK(gerr, "Trouble get parameter range for edge.");

  if (umin == umax) throw Error(MK_BAD_GEOMETRIC_EVALUATION, "Edge evaluated to some parameter umax and umin.");

  //get the arc length
  measure = ent -> measure();

  int index = 0;
  double u, du, uMid;
  du = (umax - umin)/(double)num_edges;
	
  std::vector<double> NodeCoordinates;
  std::vector<double> TempNode;
  std::vector<double> URecord;		//record the value of U

  Point3D pts0, pts1, ptsMid;
  double tmp[3];

  NodeCoordinates.resize(3*(num_edges + 1));

  TempNode.resize(3*1);
  URecord.resize(1);	

  gerr = ent -> igeom_instance() -> getEntUtoXYZ(ent->geom_handle(), umin, TempNode[0], TempNode[1], TempNode[2]);
  IBERRCHK(gerr, "Trouble getting U from XYZ along the edge");
	
  NodeCoordinates[3*0] = TempNode[0];
  NodeCoordinates[3*0+1] = TempNode[1];
  NodeCoordinates[3*0+2] = TempNode[2];

  URecord[0] = umin;

  u = umin;
  //distribute the mesh nodes on the edge based on the curvature
  for (int i = 1; i < num_edges; i++)
  {
    //first distribute the nodes evenly
    u = umin + i*du;
    gerr = ent->igeom_instance()->getEntUtoXYZ(ent->geom_handle(), u, NodeCoordinates[3*i], NodeCoordinates[3*i+1], NodeCoordinates[3*i+2]);
    IBERRCHK(gerr, "Trouble getting U from XYZ along the edge.");

    //store the two adjacent mesh nodes on the edge
    pts0.px = NodeCoordinates[3*(i-1)];
    pts0.py = NodeCoordinates[3*(i-1)+1];
    pts0.pz = NodeCoordinates[3*(i-1)+2];

    pts1.px = NodeCoordinates[3*i];
    pts1.py = NodeCoordinates[3*i+1];
    pts1.pz = NodeCoordinates[3*i+2];

    //get the coordinates for mid point between two adjacent mesh nodes on the edge
    uMid = (u-du+u)/2;
    gerr = ent->igeom_instance()->getEntUtoXYZ(ent->geom_handle(), uMid, tmp[0], tmp[1], tmp[2]);
    ptsMid.px = tmp[0];
    ptsMid.py = tmp[1];
    ptsMid.pz = tmp[2];

    //calculate the error and check whether it requires further refinement based on the curvature
    if (!ErrorCalculate(ent, pts0, pts1, ptsMid))
    {
      	DivideIntoMore(ent, pts0, ptsMid, pts1, u-du, u, uMid, index, TempNode, URecord);
    }
    
    // add the other end node to the array, get the next two adjacent mesh nodes
    index++;
    TempNode.resize(3*(index + 1));
    URecord.resize(index + 1);
    TempNode[3*index] = pts1.px;
    TempNode[3*index + 1] = pts1.py;
    TempNode[3*index + 2] = pts1.pz;

    URecord[index] = u;	
  }

  //sorting the parametrical coordinate data based on the value of u
  assert(TempNode.size()== (3*URecord.size()) );
	
  QuickSorting(TempNode, URecord, URecord.size());
  num_edges = URecord.size() - 1;
	
  //resize the variable coords
  coords.resize(3*(num_edges+1));	

  //move the other end node to the endmost of the list
  for (int i = 0; i < 3; i++)
    coords[3*num_edges+i] = coords[3*initial_num_edges+i];

  //return the mesh nodes	
  for (int i = 1; i < num_edges; i++)
  {
    coords[3*i] = TempNode[3*i];
    coords[3*i+1] = TempNode[3*i+1];
    coords[3*i+2] = TempNode[3*i+2];
  }

}

//---------------------------------------------------------------------------//
// Create the mesh for edges with dual bias distances
void EdgeMesher::DualBiasMeshing(ModelEnt *ent, int &num_edges, std::vector<double> &coords)
{
  double umin, umax, measure;

    //get the u range for the edge
  iGeom::Error gerr = ent->igeom_instance()->getEntURange(ent->geom_handle(), umin, umax);
  IBERRCHK(gerr, "Trouble get parameter range for edge.");

  if (umin == umax) throw Error(MK_BAD_GEOMETRIC_EVALUATION, "Edge evaluated to some parameter umax and umin.");

  //get the arc length
  measure = ent -> measure();

  double u, L0, dist, u0, u1;
	
  //if the node # is odd, node # will increase by 1
  if ((num_edges%2)!=0)
  {
    num_edges++;
    coords.resize(3*(num_edges+1));
      //move the other end node's position because the variable coords has been resized.
    for (int k = 0; k < 3; k++)
      coords[3*num_edges + k] = coords[3*num_edges + k - 3];
  }

  //set up the default bias ratio 1.2 temporarily	
  double q = 1.2;
	
  //get the distance between the first two nodes
  L0 = 0.5 * measure * (1-q) / (1 - pow(q, num_edges/2));
		
	
  u0 = umin;
  u1 = umax;
  //distribute the mesh nodes on the edge with dual-bias distances
  for (int i = 1; i < (num_edges/2 + 1); i++)
  {
      //distribute the mesh nodes on the one side		
    dist = L0*pow(q, i-1);
    u = u0 + (umax - umin) * dist/measure;
    u = getUCoord(ent, u0, dist, u, umin, umax);
    u0 = u;
    gerr = ent->igeom_instance()->getEntUtoXYZ(ent->geom_handle(), u, coords[3*i], coords[3*i+1], coords[3*i+2]);
    IBERRCHK(gerr, "Trouble getting U from XYZ along the edge.");

    //distribute the mesh nodes on the other side
    if (i < num_edges/2)
    {
      u = u1 - (umax-umin) * dist / measure;
      u = getUCoord(ent, u1, dist, u, umin, umax);
      gerr = ent->igeom_instance()->getEntUtoXYZ(ent->geom_handle(), u, coords[3*(num_edges-i)], coords[3*(num_edges-i)+1], coords[3*(num_edges-i)+2]);
      IBERRCHK(gerr, "Trouble getting U from XYZ along the edge");
    }
		
  }
	
}

//---------------------------------------------------------------------------//
// Create the mesh for edges with bias distances
void EdgeMesher::BiasMeshing(ModelEnt *ent, int num_edges, std::vector<double> &coords)
{
  double umin, umax, measure;

  //get the u range for the edge
  iGeom::Error gerr = ent->igeom_instance()->getEntURange(ent->geom_handle(), umin, umax);
  IBERRCHK(gerr, "Trouble get parameter range for edge.");

  if (umin == umax) throw Error(MK_BAD_GEOMETRIC_EVALUATION, "Edge evaluated to some parameter umax and umin.");

  //get the arc length
  measure = ent -> measure();

  double u, L0, dist = 0, u0;
	
  //set up the default bias ratio 1.2	
  double q = 1.2;
  L0 = measure * (1-q) / (1 - pow(q, num_edges));
		
	
  u0 = umin;
  //distribute the mesh nodes on the edge with bias distances
  for (int i = 1; i < num_edges; i++)
  {
    dist = L0*pow(q, i-1);
    u = u0 + (umax - umin)*dist/measure;
    u = getUCoord(ent, u0, dist, u, umin, umax);
    u0 = u;
    gerr = ent->igeom_instance()->getEntUtoXYZ(ent->geom_handle(), u, coords[3*i], coords[3*i+1], coords[3*i+2]);
    IBERRCHK(gerr, "Trouble getting U from XYZ along the edge.");
  }
}

//---------------------------------------------------------------------------//
// Mesh Refinement function based on the curvature: if the error is too big, refine the mesh on the edge
void EdgeMesher::DivideIntoMore(ModelEnt *ent, Point3D p0, Point3D pMid, Point3D p1, double u0, double u1, double uMid, int &index, vector<double> &nodes, vector<double> &URecord)
{
  //this is a recursive process, the process continues until the error is smaller than what is required.
  //first get two adjacent mesh nodes on the edge and coordinates for mid point between two adjacent nodes.
  //then check the left side and right side whether the error is too big nor not
  double uu0, uu1, uumid, tmp[3];
  Point3D pts0, pts1, ptsMid;
	
  index++;
  nodes.resize(3*(index+1));
  URecord.resize(index+1);
  nodes[3*index] = pMid.px;
  nodes[3*index+1] = pMid.py;
  nodes[3*index+2] = pMid.pz;
  URecord[index] = uMid;
	
  //check the left side
  uu0=u0;
  uu1=uMid;
  uumid=(uu0+uu1)/2;
  pts0=p0;
  pts1=pMid;


  iGeom::Error gerr = ent->igeom_instance()->getEntUtoXYZ(ent->geom_handle(), uumid, tmp[0], tmp[1], tmp[2]);
  IBERRCHK(gerr, "Trouble getting U from XYZ along the edge.");
  ptsMid.px = tmp[0];
  ptsMid.py = tmp[1];
  ptsMid.pz = tmp[2];

  //check the error
  if(!ErrorCalculate(ent, pts0, pts1, ptsMid))
  {
    //further refinement
    DivideIntoMore(ent, pts0, ptsMid, pts1, uu0, uu1, uumid, index, nodes, URecord);
  }
	
  //check the right side
  uu0 = uMid;
  uu1=u1;
  uumid=(uu0+uu1)/2;
  pts0=pMid;
  pts1=p1;
  //get the coorindinates for mid point 
  gerr = ent->igeom_instance()->getEntUtoXYZ(ent->geom_handle(), uumid, tmp[0], tmp[1], tmp[2]);
  IBERRCHK(gerr, "Trouble getting U from XYZ along the edge.");
  ptsMid.px = tmp[0];
  ptsMid.py = tmp[1];
  ptsMid.pz = tmp[2];
	
  //check the error
  if(!ErrorCalculate(ent, pts0, pts1, ptsMid))
  {
    //further refinement
    DivideIntoMore(ent, pts0, ptsMid, pts1, uu0, uu1, uumid, index, nodes, URecord);
  }
			
}


//---------------------------------------------------------------------------//
// Rapid sorting the mesh nodes on the edge based on the parametric coordinates. This is a recursive
// process
void EdgeMesher::RapidSorting(vector<double> &nodes, vector<double> &URecord, int left, int right)
{
  int i, j;
  double middle, iTemp;
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

//---------------------------------------------------------------------------//
// Quick Sorting: this function comes together with the function RapidSorting
void EdgeMesher::QuickSorting(vector<double> &nodes, vector<double> &URecord, int count)
{
  RapidSorting(nodes, URecord, 0, count-1);
}


//---------------------------------------------------------------------------//
// return the x, y, z coordinates from the parametric coordinates
EdgeMesher::Point3D EdgeMesher::getXYZCoords(ModelEnt *ent, double u) const
{
  Point3D pts3D;
  double xyz[3];


  //get the coordinates in the physical space
  iGeom::Error gerr = ent->igeom_instance()->getEntUtoXYZ(ent->geom_handle(), u, xyz[0], xyz[1], xyz[2]);
  IBERRCHK(gerr, "Trouble getting U from XYZ along the edge.");	

  pts3D.px = xyz[0];
  pts3D.py = xyz[1];
  pts3D.pz = xyz[2];
  return pts3D;
}

//---------------------------------------------------------------------------//
// get the parametric coordinate based on first parametric coordinate ustart and distance in the physical space
double EdgeMesher::getUCoord(ModelEnt *ent, double ustart, double dist, double uguess, double umin, double umax) const
{

  Point3D p0 = getXYZCoords(ent, ustart);
  Point3D p1 = getXYZCoords(ent, uguess);


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
    p1 = getXYZCoords(ent, u);
		

    if (ntrials++ == 100000)
    {
      cout << " Warning: Searching for U failed " << endl;
    }
  }
  uguess = u;
  return uguess;
}

//---------------------------------------------------------------------------//
// calculate the error: the distance between the mid point of two adjacent 
// mesh nodes (on the mesh line segments) and mid point on the edge
bool EdgeMesher::ErrorCalculate(ModelEnt *ent, Point3D p0, Point3D p1, Point3D pMid)
{
  double lengtha, lengthb, lengthc;
  double deltax, deltay, deltaz;
  double angle, error, tol=1.0E-3, H;
  double cvtr_ijk[3], curvature;
  bool result;

  //calculate the distance between the first mesh node and mid point on the edge
  deltax = pMid.px-p0.px;
  deltay = pMid.py-p0.py;
  deltaz = pMid.pz-p0.pz;	
  lengtha = sqrt(deltax*deltax + deltay*deltay + deltaz*deltaz);

  //calculate the distance between two adjacent mesh nodes
  deltax = p1.px-p0.px;
  deltay = p1.py-p0.py;
  deltaz = p1.pz-p0.pz;	
  lengthb = sqrt(deltax*deltax + deltay*deltay + deltaz*deltaz);

  //calculate the distance between the second mesh node and mid point on the edge
  deltax = pMid.px-p1.px;
  deltay = pMid.py-p1.py;
  deltaz = pMid.pz-p1.pz;	
  lengthc = sqrt(deltax*deltax + deltay*deltay + deltaz*deltaz);

  //calculate the angle
  angle = acos((lengtha*lengtha + lengthb*lengthb - lengthc*lengthc)/(2*lengtha*lengthb));
  H = fabs(lengtha*sin(angle));

  //calculate the curvature	
  iGeom::Error gerr = ent->igeom_instance()->getEgCvtrXYZ(ent->geom_handle(), pMid.px, pMid.py, pMid.pz, cvtr_ijk[0], cvtr_ijk[1], cvtr_ijk[2]);
  IBERRCHK(gerr, "Trouble getting U from XYZ along the edge.");		
  curvature = sqrt(cvtr_ijk[0]*cvtr_ijk[0]+cvtr_ijk[1]*cvtr_ijk[1]+cvtr_ijk[2]*cvtr_ijk[2]);
  error= H*curvature;
	
	
  if (error > tol)
    result = false;
  else
    result = true;
  return result;		
}

}

