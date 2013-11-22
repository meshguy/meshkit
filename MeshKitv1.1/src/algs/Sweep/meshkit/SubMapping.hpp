//-----------------------------------C++-------------------------------------//
// File: src/algs/Submapping.hpp
// Wednesday February 11 10:50 2011
// Brief: SubMapping class definition: generate the structural mesh by
//        submapping 
//---------------------------------------------------------------------------//


#ifndef MESHKIT_SUBMAPPING_HPP
#define MESHKIT_SUBMAPPING_HPP

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
#include "Global.hpp"
#include <iMesh.h>
#include <iGeom.h>
#include <set>
#include <iRel.h>
#include <vector>
#include <set>
#include <list>

#include "meshkit/MeshScheme.hpp"

using namespace std;

namespace MeshKit
{
//===========================================================================//
  /*!
   * \class SubMapping
   * \brief Generate the structural mesh for linking surface
   * 
   * SubMapping generates the structural mesh
   * 3 steps: vertex classification, assign each node on the bounary i-j coordinates
              subdivide the surface
   */         
//===========================================================================//
enum VertexTypes{REVERSAL=-2, CORNER=-1, SIDE = 0, END = 1};
enum EdgeTypes{POSI_I = 0, NEG_I= 1, POSI_J = 2, NEG_J = 3};
class SubMapping :  public MeshScheme
{	
public:

	SubMapping(MKCore *mk_core, const MEntVector &me_vec);
	~SubMapping();

	virtual void setup_this();
	virtual void execute_this();

 
        /**\brief Get class name */
        static const char* name() 
          { return "SubMapping"; }

        /**\brief Function returning whether this scheme can mesh entities of t
         *        the specified dimension.
         *\param dim entity dimension
         */
        static bool can_mesh(iBase_EntityType dim)
          { return iBase_FACE == dim; }

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
          { return output_types();}

		//setup the mesh size
		void SetupMeshSize(double size = 1.0);

private:
		//classify the vertices as side, end, reversal and corner
		void VertexClassification(ModelEnt *ent);

		//classify the boundary edges as -I, +I, -J, +J
		void EdgeClassification();

		//assign the i,j coordinates for each node on the boundary
		void EdgeDiscretization(ModelEnt *me);

		//ReOrganize the vertices and edges on the boundary
		void VerEdgOrganize(std::set<iBase_EntityHandle> edge_set, std::vector<iBase_EntityHandle> g_edge, iBase_EntityHandle surf);

		//subdivide the surface 
		void InteriorNodeInterpolation(ModelEnt *me);

		//calculate the angle for each vertex
		void GetAngle(iBase_EntityHandle surf, vector<double> &angle);

		//assign the i-j coordinates for boundary nodes
		void AssignIJCoords(double &u, double &v, EdgeTypes type, int index);

		//call the linear programming library to solve the linear programm
		void VtxClassificationLP();

		//check whether a point is outside the surface or not
		bool isOutSideSurf(vector<Vertex> corner, int i, int j);
		double Angle2D(double x1, double y1, double x2, double y2);

		//Use the Winslow to smooth the structured mesh
		void MeshSmoothing(ModelEnt *ent);

		void buildAssociation();

		bool isCurved(int vtx_index, vector<double> u1, vector<double> u2, vector<double> u3, vector<double> u4, vector<vector<double> > tang_pre, vector<vector<double> > tang_next);

		void LinearProgramming();

		
private://private member variable
	vector<VertexTypes> vertices_types;
	vector<Vertex> nodes, vertices;
	vector<Edge> edges;
	
	vector<double> interior_angle;
	vector<int> sorted_vertex_list, sorted_node_list, sorted_edge_list;
	vector<EdgeTypes> edges_types;	
	
	iBase_TagHandle g_taghandle, m_taghandle;

	vector<vector<int> > coordinate_i_j;

	double size_low_bound;

	unsigned int start_index;

	vector<int> edge_size;
	map<int, int> geom_mesh_node;
	map<int, int> mesh_geom_vertex;

	vector<Face> quads;

	
};

}

#endif


