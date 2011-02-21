//-----------------------------------C++-------------------------------------//
// File: algs/SCDMesh.hpp
// Author: Stuart R. Slattery
// Wednesday February 2 16:15:8 2011
// Brief: SCDMesh class definition 
//---------------------------------------------------------------------------//

#ifndef __SCDMESH_HPP
#define __SCDMESH_HPP

#include <vector>

#include "meshkit/Types.hpp"
#include "meshkit/Error.hpp"
#include "meshkit/MeshScheme.hpp"
#include "meshkit/ModelEnt.hpp"
#include "meshkit/iGeom.hpp"

#include "moab/Core.hpp"
#include "moab/Range.hpp"
#include "moab/Interface.hpp"
#include "moab/ScdInterface.hpp"

namespace MeshKit
{

//===========================================================================//
  /*!
   * \class SCDMesh
   * \brief Generate a structured Cartesian grid over a geometry
   * 
   * SCDMesh generates a structured Cartesian grid with either a lightweight
   * representation using ScdInterface in MOAB or a full representation
   * of all the entities. 
   */
//===========================================================================//


  class SCDMesh : public MeshScheme
  {
  public:

    // 2 different interface types
    // full - full element and vertex representation in moab
    // scd - minimal structured mesh represenation using ScdInterface in MOAB
    enum InterfaceSchemeType { full, scd};

    // 2 different grid scheme types
    // cfMesh - provide a coarse number of i, j, k cells and a number of fine
    //          sub-cells in a given coarse division
    // vtxMesh - provide an exact list of the vertex coordinates to be used
    //           in each direction
    enum GridSchemeType { cfMesh, vtxMesh};

    // 2 different axis scheme types
    // cartesian - i, j, k directions are parallel to the x, y, z axes
    // oriented - i, j, k directions oriented with the primary geometric axes
    enum AxisSchemeType { cartesian, oriented};

  private:

    // enum types defining the mesh generation
    InterfaceSchemeType interfaceType;
    GridSchemeType gridType;
    AxisSchemeType axisType;

    // number of coarse divisions in each direction
    int coarse_i;
    int coarse_j;
    int coarse_k;

    // array of fine divisions in each division
    std::vector<int> fine_i;
    std::vector<int> fine_j;
    std::vector<int> fine_k;

    // vertex arrays in each direction
    std::vector<double> i_arr;
    std::vector<double> j_arr;
    std::vector<double> k_arr;

    // error codes
    iGeom::Error gerr;
    moab::ErrorCode rval;

    // bounding box min, max coordiantes
    double minCoord[3], maxCoord[3];


  public:

    //! Bare Constructor
    SCDMesh(MKCore *mkcore, const MEntVector &me_vec);

    //! Destructor
    virtual ~SCDMesh();

    /*! 
     * \brief Setup
     */
    virtual void setup_this();

    //! Generates the SCD Mesh
    virtual void execute_this();

   /**\brief Get class name */
    static const char* name() 
      { return "SCDMesh"; }
      
    /**\brief Function returning whether this scheme can mesh entities of
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
    static bool can_mesh(ModelEnt* ent) 
      { return canmesh_region(ent); }
    
    /**\brief Get list of mesh entity types that can be generated.
     *\return array terminated with \c moab::MBMAXTYPE
     */
    static const moab::EntityType* output_types();
  
    /** \brief Return the mesh entity types operated on by this scheme
     * \return array terminated with \c moab::MBMAXTYPE
     */
    virtual const moab::EntityType* mesh_types_arr() const
      { return output_types(); }

    /*!
     * \brief Set the interface type to be used
     * \param scheme The type of grid scheme to be used
     */
    void set_interface_scheme(InterfaceSchemeType scheme);

    /*!
     * \brief Get the interface type
     */
    InterfaceSchemeType get_interface_scheme() const;

    /*!
     * \brief Set the grid generation scheme to be used
     * \param scheme The type of grid scheme to be used
     */
    void set_grid_scheme(GridSchemeType scheme);

    /*!
     * \brief Get the grid scheme assigned to the mesh
     */
    GridSchemeType get_grid_scheme() const;

    /*!
     * \brief Set the axis type to be used
     * \param scheme The type of axis scheme to use
     */
    void set_axis_scheme(AxisSchemeType scheme);

    /*!
     * \brief Get the axis scheme assigned ot the mesh
     */
    AxisSchemeType get_axis_scheme() const;

    /*!
     * \brief Set the i direction number of coarse grid divisions for the cfmesh case
     */
    void set_coarse_i_grid(int coarse_i_);

    /*!
     * \brief Set the j direction number of coarse grid divisions for the cfmesh case
     */
    void set_coarse_j_grid(int coarse_j_);

    /*!
     * \brief Set the k direction number of coarse grid divisions for the cfmesh case
     */
    void set_coarse_k_grid(int coarse_k_);

    /*!
     * \brief Set the i direction fine grid divisions
     */
    void set_fine_i_grid(std::vector<int> fine_i_);

    /*!
     * \brief Set the j direction fine grid divisions
     */
    void set_fine_j_grid(std::vector<int> fine_j_);

    /*!
     * \brief Set the k direction fine grid divisions
     */
    void set_fine_k_grid(std::vector<int> fine_k_);

     /*!
     * \brief Set the geometry encompasing box dimensions
     */
    void set_box_dimension();

    /*!
     * \brief Get the geometry encompasing box dimensions
     */
    void get_box_dimension(double* min, double* max);


  private:
    
    //! create cartesian bounding box
    void create_cart_edges();

    //! create full mesh representation
    void create_full_mesh();
    
  }; // end class SCDMesh


  /* Inline Functions */

  inline SCDMesh::~SCDMesh()
  {
  }

  inline void SCDMesh::set_interface_scheme(SCDMesh::InterfaceSchemeType scheme)
  {
    interfaceType = scheme;
  }

  inline SCDMesh::InterfaceSchemeType SCDMesh::get_interface_scheme() const
  {
    return interfaceType;
  }

  inline void SCDMesh::set_grid_scheme(SCDMesh::GridSchemeType scheme)
  {
    gridType = scheme;
  }

  inline SCDMesh::GridSchemeType SCDMesh::get_grid_scheme() const
  {
    return gridType;
  }

  inline void SCDMesh::set_axis_scheme(SCDMesh::AxisSchemeType scheme)
  {
    axisType = scheme;
  }

  inline SCDMesh::AxisSchemeType SCDMesh::get_axis_scheme() const
  {
    return axisType;
  }

  inline void SCDMesh::set_coarse_i_grid(int coarse_i_)
  {
    coarse_i = coarse_i_;
  }

  inline void SCDMesh::set_coarse_j_grid(int coarse_j_)
  {
    coarse_j = coarse_j_;
  }

  inline void SCDMesh::set_coarse_k_grid(int coarse_k_)
  {
    coarse_k = coarse_k_;
  }

  inline void SCDMesh::set_fine_i_grid(std::vector<int> fine_i_)
  {
    fine_i = fine_i_;
  }

  inline void SCDMesh::set_fine_j_grid(std::vector<int> fine_j_)
  {
    fine_j = fine_j_;
  }

  inline void SCDMesh::set_fine_k_grid(std::vector<int> fine_k_)
  {
    fine_k = fine_k_;
  }

  inline void SCDMesh::get_box_dimension(double* min, double* max)
  {
    for (int i = 0; i < 3; i++) {
      min[i] = minCoord[i];
      max[i] = maxCoord[i];
    }
  }

} // end namespace MeshKit

#endif // end if __SCDMESH_HPP

//---------------------------------------------------------------------------//
// end algs/SCDMesh.hpp
//---------------------------------------------------------------------------//
