//-----------------------------------C++-------------------------------------//
// File: src/algs/TFIMapping.hpp
// Wednesday February 11 10:50 2011
// Brief: TFIMapping class definition: generate the all-quad mesh by
//       TFI mapping. Before using this class, four edges should have already
//	 been meshed by EdgeMesher
//---------------------------------------------------------------------------//


#ifndef MESHKIT_TFIMAPPING_HPP
#define MESHKIT_TFIMAPPING_HPP

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
#include "meshkit/SizingFunction.hpp"
#include "moab/ReadUtilIface.hpp"
#include <iMesh.h>
#include <set>
#include <iRel.h>
#include <vector>

#include "meshkit/MeshScheme.hpp"

using namespace std;

namespace MeshKit
{
//===========================================================================//
  /*!
   * \class TFIMapping
   * \brief Generate the all-Quad mesh by TFI mapping
   * 
   * TFIMapping generates the all-quad mesh by trans-finite interpolation with
   * i, j parameters
   */
//===========================================================================//
class TFIMapping :  public MeshScheme
{	
public:

	TFIMapping(MKCore *mk_core, const MEntVector &me_vec);
	~TFIMapping();

	virtual void setup_this();
	virtual void execute_this();

 
        /**\brief Get class name */
        static const char* name() 
          { return "TFIMapping"; }

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
          { return output_types(); }


private:

	//interpolate linearly between x0 and x1
	double linear_interpolation(double r, double x0, double x1);

	//implement the transfinite interpolation between (pt_0s, pt_1s) and (pt_r0, pt_r1)
	double parametricTFI2D(double r, double s, double pt_0s, double pt_1s, double pt_r0, double pt_r1);

	//generate the mesh on the linking surface
	int SurfMapping(ModelEnt *ent);

	
private://private member variable
	
	
};

}

#endif


