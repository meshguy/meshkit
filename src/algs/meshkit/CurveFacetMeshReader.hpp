

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string>
#include <iostream>
#include <fstream>
#include <string.h>
#include <limits.h>

#include <iGeom.h>
#include <iMesh.h>
#include <set>
#include <iRel.h>
#include <vector>
#include "meshkit/MKCore.hpp"
#include "meshkit/MeshScheme.hpp"

namespace MeshKit
{
 
using namespace std;

  /** \class SolidCurveMesher SolidCurveMesher.hpp "meshkit/SolidCurveMesher.hpp"
   * \brief A class for generating facet-based mesh of geometric curves.
   * INPUT: one or more ModelEnts representing geometric curves
   * MESH TYPE(S): MBEDGE, MBVERTEX
   * OUTPUT: a set of edges and vertices for each ModelEnt
   * DEPENDENCIES: meshkit build with iGeom interface
   * 
   * This class uses the facets generated for the visualization of solid models in typical CAD 
   * software to represent a geometric curve as mesh. Upon execution, this class will call for
   * the facet data and store it as part of the ModelEnt's mesh.
   */

class CurveFacetMeshReader : public MeshScheme
{
public: 
       CurveFacetMeshReader(MKCore *mk, const MEntVector &ments);

       ~CurveFacetMeshReader();
private:
        double facet_tol;
        double geom_res;
        MKCore *mk;

        /** \brief Returns the distance between an iGeom vertex and an iMesh vertex.
         * \param vtx1 iGeom vertex handle
         * \param vtx2 iMesh vertex handle
         */
        virtual double vtx2vtx_dist(iGeom::EntityHandle vtx1, iMesh::EntityHandle vtx2);

        /** \brief Returns the distance between an iMesh vertex and an iMesh vertex.
         * \param vtx1 iMesh vertex handle
         * \param vtx2 iMesh vertex handle
         */
        virtual double mvtx2mvtx_dist(iMesh::EntityHandle vtx1, iMesh::EntityHandle vtx2);


        virtual void facet(ModelEnt *curve);

        /** \brief Sets the senses wrt all surfaces adjacent to the curve
         * \param curve Pointer to the ModelEnt to be meshed
         */
        virtual void set_senses( ModelEnt *curve);

public:  
       virtual void setup_this();
       virtual void execute_this();

       /** \brief Sets the faceting tolerance and geom_reabs values. If
        *    this function is not run before mk->setup(). Default values
        *    for these parameters will be used. 
        * \param faceting_tolerance value to be set for the faceting tolerance
        * \param geom_resabs value to be set for the geom_resabs (used for vertex proximity checks)
        */
       void set_mesh_params(double faceting_tolerance = 0, double geom_resabs = 0);

       /** \brief Function returning whether this scheme can mesh entities of the 
        *  specified dimension
        *  \param dim entity dimension 
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

      /** \brief Get the class name */
      static const char* name()
      {return "CurveFacetMeshReader";}

      /** \brief Function that names the output types of this meshing class
       *  \return array terminated with \c moab::MBMAXTYPE
       */
      static const moab::EntityType* output_types();

      virtual const moab::EntityType* mesh_types_arr() const
      { return output_types(); }

};


}
