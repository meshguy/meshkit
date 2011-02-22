//-----------------------------------C++-------------------------------------//
// File: src/algs/EdgeMesher.hpp
// Wednesday February 11 10:50 2011
// Brief: EdgeMesher class definition: four schemes are provided: equal meshing,
//        Bias Meshing, Dual Bias Meshing, Curvature-based meshing 
//---------------------------------------------------------------------------//

#ifndef MESHKIT_EDGE_MESHER_HPP
#define MESHKIT_EDGE_MESHER_HPP

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string>
#include <iostream>
#include <fstream>
#include <string.h>
#include <limits.h>

#include "meshkit/iGeom.hpp"
#include <set>
#include <vector>

#include "meshkit/MeshScheme.hpp"

namespace MeshKit
{
//===========================================================================//
  /*!
   * \class EdgeMesher
   * \brief Generate the mesh for edges
   * 
   * EdgeMesher generates a mesh for edges, creating the nodes and line segments
   * on edges. There are four schemes for edge mesher: equal meshing, Bias Meshing
   * Dual Bias Meshing, Curvature-Based Meshing
   */
//===========================================================================//

using namespace std;

class EdgeMesher : public MeshScheme
{
public:
	//Four schemes for edge meshing
	enum EdgeSchemeType {EQUAL=0, BIAS, DUAL, CURVATURE};

	
public:
	//construction function for edge mesher
	EdgeMesher(MKCore *mk_core, const MEntVector &me_vec);

	//set up the parameters for edge meshing, e.g. compute the number of intervals
	virtual void setup_this();

	//Generate the edge mesh
	virtual void execute_this();

	//make an instance of the EdgeMesher class
	EdgeSchemeType get_edge_scheme() const;
	 void set_edge_scheme(EdgeSchemeType scheme);

       /**\brief Get class name */
       static const char* name() 
         { return "EdgeMesher"; }

       /**\brief Function returning whether this scheme can mesh entities of t
        *        the specified dimension.
        *\param dim entity dimension
        */
       static bool can_mesh(iBase_EntityType dim)
          { return iBase_EDGE == dim; }

       /** \brief Function returning whether this scheme can mesh the specified entity
        * 
        * Used by MeshOpFactory to find scheme for an entity.
        * \param me ModelEnt being queried
        * \return If true, this scheme can mesh the specified ModelEnt
        */
       static bool can_mesh(ModelEnt *me)
          { return canmesh_edge(me); }


       /**\brief Get list of mesh entity types that can be generated.
        *\return array terminated with \c moab::MBMAXTYPE
        */
       static const moab::EntityType* output_types();

       /** \brief Return the mesh entity types operated on by this scheme
        * \return array terminated with \c moab::MBMAXTYPE
        */
       virtual const moab::EntityType* mesh_types_arr() const
         { return output_types(); }

	~EdgeMesher();
	
private:

        struct Point3D
        {
	        double px;
	        double py;
	        double pz;	
        };

 	EdgeSchemeType schemeType;

	//return x, y, z coordinates based on the parametric coordinate u on the edge
	Point3D getXYZCoords(ModelEnt *ent, double u) const;
	
	//return the parametrical coordinate based the starting parametric coordinate ustart and distance in physical space
	double getUCoord(ModelEnt *ent, double ustart, double dist, double uguess, double umin, double umax) const;

	//Create more nodes if there is a high curvature on the edge
	void DivideIntoMore(ModelEnt *ent, Point3D p0, Point3D pMid, Point3D p1, double u0, double u1, double uMid, int &index, vector<double> &nodes, vector<double> &URecord);

	//calculate the error between the line segments and the edge in the physical space
	bool ErrorCalculate(ModelEnt *ent, Point3D p0, Point3D p1, Point3D pMid);

	//Sort the nodes on the edge
	void RapidSorting(vector<double> &nodes, vector<double> &URecord, int left, int right);
	void QuickSorting(vector<double> &nodes, vector<double> &URecord, int count);

	//four schemes for edge meshing
	//create the mesh for edges with equal distances
	void EqualMeshing(ModelEnt *ent, int num_edges, std::vector<double> &coords);

	//create the mesh for edges with bias distances
	void BiasMeshing(ModelEnt *ent, int num_edges, std::vector<double> &coords);

	//create the mesh for edges with dual bias distances
	void DualBiasMeshing(ModelEnt *ent, int &num_edges, std::vector<double> &coords);

	//create the mesh for edges based on curvatures 
	void CurvatureMeshing(ModelEnt *ent, int &num_edges, std::vector<double> &coords);

	//compute the distance between the parametric coordinate ustart and parametric coordinate uend.
	double measure(iGeom::EntityHandle ent, double ustart, double uend) const;

};

//set up the scheme type for edge meshing
inline void EdgeMesher::set_edge_scheme(EdgeMesher::EdgeSchemeType scheme)
{
	schemeType = scheme;
}

//return the scheme type for edge meshing
inline EdgeMesher::EdgeSchemeType EdgeMesher::get_edge_scheme() const
{
	return schemeType;
}
}

#endif
