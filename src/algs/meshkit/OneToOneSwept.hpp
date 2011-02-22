//-----------------------------------C++-------------------------------------//
// File: src/algs/OneToOneSwept.hpp
// Wednesday February 11 10:50 2011
// Brief: OneToOneSwept class definition: generate the all-hexahedral mesh by
//        general sweeping 
//---------------------------------------------------------------------------//


#ifndef MESHKIT_ONETOONESWEPT_HPP
#define MESHKIT_ONETOONESWEPT_HPP

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string>
#include <iostream>
#include <fstream>
#include <string.h>
#include <limits.h>

#include <iGeom.h>
#include "meshkit/MKCore.hpp"
#include "meshkit/ModelEnt.hpp"
#include "meshkit/SizingFunction.hpp"
#include "moab/ReadUtilIface.hpp"
#include <iMesh.h>
#include <set>
#include <iRel.h>
#include <vector>

#include "meshkit/MeshScheme.hpp"

#include "SimpleArray.hpp"

using namespace std;

namespace MeshKit
{
//===========================================================================//
  /*!
   * \class OneToOneSwept
   * \brief Generate the all-hexahedral mesh by sweeping
   * 
   * OneToOneSwept generates the all-hexahedral mesh by sweeping the source mesh
   * to the target surface
   */
//===========================================================================//

struct Point2D
{
	double pu;
	double pv;	
};
struct Point3D
{
	double px;
	double py;
	double pz;
};
struct Vertex
{
   	Vertex() { onBoundary = 0; }
   	int    id;
	int index;
   	bool   onBoundary;
   	bool   onCorner;
   	double uCoord, uvCoords[2], xyzCoords[3];
	iBase_EntityHandle gVertexHandle;
};
struct Edge
{
  	Vertex* connect[2];
	int  EdgeID;
	int index;
	bool onBoundary;
	iBase_EntityHandle gEdgeHandle;
};

struct Face
{
    	int  getNumNodes() const { return connect.size(); }
    	Vertex* getVertex(int i) const { return connect[i]; }
    	vector<Vertex*> connect;
    	int FaceID;
	int index;
	iBase_EntityHandle gFaceHandle;
};

class OneToOneSwept :  public MeshScheme
{	
public:

	OneToOneSwept(MKCore *mk_core, const MEntVector &me_vec);
	~OneToOneSwept();

	virtual void setup_this();
	virtual void execute_this();

 
        /**\brief Get class name */
        static const char* name() 
          { return "OneToOneSwept"; }

        /**\brief Function returning whether this scheme can mesh entities of t
         *        the specified dimension.
         *\param dim entity dimension
         */
        static bool can_mesh(iBase_EntityType dim)
          { return iBase_REGION == dim; }

        /** \brief Function returning whether this scheme can mesh the specified entity
         * 
         * Used by MeshOpFactory to find scheme for an entity.
         * \param me ModelEnt being queried
         * \return If true, this scheme can mesh the specified ModelEnt
         */
        static bool can_mesh(ModelEnt *me)
          { return canmesh_region(me); }

        /**\brief Get list of mesh entity types that can be generated.
         *\return array terminated with \c moab::MBMAXTYPE
         */
        static const moab::EntityType* output_types();

        /** \brief Return the mesh entity types operated on by this scheme
         * \return array terminated with \c moab::MBMAXTYPE
         */
        virtual const moab::EntityType* mesh_types_arr() const
          { return output_types(); }


private:
	//Build Association function: try to build the association between the geometry and mesh
	void buildAssociation();

	//determine whether a mesh edge is on the boundary or not
	int isEdgeBoundary(iBase_EntityHandle gEdgeHandle);	

	//find the corner node list, inner node list and edge node list for the mesh on the source surface
	void GetList();
	
	//specify the source surface for OneToOneSwept class
	void SetSourceSurface();
	
	//specify the target surface for OneToOneSwept class
	void SetTargetSurface();
	
	//function for all-quad meshing on the target surface
	int TargetSurfProjection();

	//function for obtaining the parametric coordinates from x,y,z coordinates
	int getUVCoords(iBase_EntityHandle gFaceHandle, Point3D pts3, Point2D &pts2);

	//function for obtaining the x,y,z coordinates from parametric coordinates
	int getXYZCoords(iBase_EntityHandle gFaceHandle, Point3D &pts3, double uv[2]);

	//interpolate linearly between x0 and x1
	double linear_interpolation(double r, double x0, double x1);

	//implement the transfinite interpolation between (pt_0s, pt_1s) and (pt_r0, pt_r1)
	double parametricTFI2D(double r, double s, double pt_0s, double pt_1s, double pt_r0, double pt_r1);

	//generate the mesh on the linking surface
	int LinkSurfMeshing(vector<vector <Vertex> > &linkVertexList);

	//generate the all-hex mesh by sweeping in the interior
	int InnerLayerMeshing();

	//generate the nodes between the Source surface and target surface.
	int InnerNodesProjection(vector<vector <Vertex> > &linkVertexList);

	//create the hexahedral elements between the source surface and target surface
	int CreateElements(vector<vector <Vertex> > &linkVertexList);
private://private member variable
	iBase_EntityHandle sourceSurface;
	iBase_EntityHandle targetSurface;
	std::vector<Edge> gsEdgeList;  //geometrical edges on the source surfaces
	std::vector<Edge> gtEdgeList;  //geometrical edges on the target surfaces
	std::vector<Vertex> NodeList;  //mesh nodes on the source surface
	map<int, int> edgePairs;  //store the relationship between the  edge id on the source surface and target surface
	map<int, int> cornerPairs; ////store the relationship between the corner vertex id on the source surface and target surface	
	std::vector<Vertex> gVertexList; //geometrical vertex
	std::vector<Face> gLinkFaceList;  //geometrical face list for linking surface
	std::vector<Vertex> TVertexList;
	std::vector<Edge> TEdgeList;
	std::vector<Face> TFaceList;
	std::vector<Edge> EdgeList;
	std::vector<Face> FaceList;
	vector<Edge> gLinkSides;  //geometrical edges for linking sides between source and target
	int numLayers;
	iBase_TagHandle  geom_id_tag, mesh_id_tag;
	iBase_EntityHandle volEntity;
	iBase_EntitySetHandle volumeSet;	

	
	
};

}

#endif


