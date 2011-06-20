//-----------------------------------C++-------------------------------------//
// File: src/algs/meshkit/AssyGen.hpp
//
// Brief: AssyGen class definition:
//        Creates reactor assembly geometry (ACIS or OCC format) and cubit mesh script as specified in a user defined input
//         class, AssyGen
//---------------------------------------------------------------------------// 

#ifndef MESHKIT_ASSYGEN_HPP
#define MESHKIT_ASSYGEN_HPP

#include <cassert>
#include <string>
#include <vector>
#include <set>

#include "meshkit/Types.hpp"
#include "meshkit/Error.hpp"
#include "meshkit/MeshScheme.hpp"
#include "meshkit/ModelEnt.hpp"
#include "meshkit/MKCore.hpp"
#include "meshkit/SizingFunction.hpp"
#include "meshkit/RegisterMeshOp.hpp"

#include "meshkit/LocalSet.hpp"
#include "meshkit/LocalTag.hpp"
#include "meshkit/Matrix.hpp"
#include "../CopyUtils.hpp"

#include "meshkit/iMesh.hpp"
#include "meshkit/iGeom.hpp"
#include "MBCN.h"

#include "../AssyGen/vectortemplate.hpp"
#include "../AssyGen/matrixtemplate.hpp"
#include "../AssyGen/pincell.hpp"
#include "../AssyGen/parser.hpp"
#include "SimpleArray.hpp"
#include "../AssyGen/clock.hpp"


namespace MeshKit {

#define DEFAULT_TEST_FILE  "assygen_default"
#define TEST_FILE_NAME "assygen_default"

  enum ErrorStates {PINCELLS, INVALIDINPUT, EMAT, EGEOMTYPE, EGEOMENGINE, ENEGATIVE, EALIAS, EPIN};

  class MKCore;

  class AssyGen : public MeshScheme
  {
  public:
    /* \brief Constructor
     *
     * Create a new AssyGen instance
     * \param impl the iGeom instance handle for the Geom
     */
    AssyGen(MKCore *mk, const MEntVector &me_vec);

    /* \brief Destructor
     */
    virtual ~AssyGen();

    enum ErrorStates {PINCELLS, INVALIDINPUT, EMAT, EGEOMTYPE, EGEOMENGINE, ENEGATIVE, EALIAS, EPIN};

    /**\brief Get class name */
    static const char* name();

    /**\brief Function returning whether this scheme can mesh entities of t
     *        the specified dimension.
     *\param dim entity dimension
     */
    static bool can_mesh(iBase_EntityType dim);

    /** \brief Function returning whether this scheme can mesh the specified entity
     *
     * Used by MeshOpFactory to find scheme for an entity.
     * \param me ModelEnt being queried
     * \return If true, this scheme can mesh the specified ModelEnt
     */
    static bool can_mesh(ModelEnt *me);

    /**\brief Get list of mesh entity types that can be generated.
     *\return array terminated with \c moab::MBMAXTYPE
     */
    static const moab::EntityType* output_types();

    /** \brief Return the mesh entity types operated on by this scheme
     * \return array terminated with \c moab::MBMAXTYPE
     */
    virtual const moab::EntityType* mesh_types_arr() const;

    /** \brief Re-implemented here so we can check topological dimension of model_ent
     * \param model_ent ModelEnt being added
     */
    virtual bool add_modelent(ModelEnt *model_ent);

    //! Setup is a no-op, but must be provided since it's pure virtual
    virtual void setup_this();

    //! The only setup/execute function we need, since meshing vertices is trivial
    virtual void execute_this();

    /** \brief Prepare input/output files for reading/writing
     * \param command line args and testdir for default test case
     */
    void PrepareIO (int argc, char *argv[], std::string TestDir);

    /** \brief Read the command based text input file
     * \param input file
     */
    void ReadInputPhase1 ();

    /** \brief Keep reading input file and create
     * \param input file
     */
    void ReadAndCreate ();

    /** \brief Name the surface created
     * \param material name from input file, surface entity, name tag
     */
    void Name_Faces( const std::string sMatName, const iBase_EntityHandle body,
		     iBase_TagHandle this_tag);

    /** \brief Move the assembly to the center
     * \param direction
     */
    void Center_Assm( char&);

    /** \brief Section assembly
     * \param direction, offset, reverse/forward
     */
    void Section_Assm ( char&, double&, const std::string);

    /** \brief Rotate assembly
     * \param direction, angle
     */
    void Rotate_Assm ( char&, double&);

    /** \brief Move assembly
     * \param X, Y, Z distance
     */
    void Move_Assm ( double&, double&, double&);

    /** \brief Create hexagonal assembly
     * \param data from input file
     */
    void Create_HexAssm( std::string &);

    /** \brief Create cartesian or rectangular assembly
     * \param data from input file 
     */
    void Create_CartAssm( std::string &);

    /** \brief Create outermost ducts
     * \param data from input file 
     */
    void CreateOuterCovering();

    /** \brief Merge and impring the geometry creaed
     * \param geometry created
     */
    void Imprint_Merge ();

    /** \brief Subtract the pins from innermost duct
     * \param geometry entities
     */
    void Subtract_Pins ();

    /** \brief Get the top surface from 3D assembly geometry
     * \param pins, ducts
     */
    void Create2DSurf();

    /** \brief Read pincell data
     * \param input file
     */
    void ReadPinCellData( int i);

    /** \brief Create pincell i, pincell intersects the assembly
     * \param i and location
     */
    void CreatePinCell_Intersect( int i, double dX,
				  double dY, double dZ);

    /** \brief Create pincell i
     * \param i and location
     */
    void CreatePinCell( int i, double dX,
			double dY, double dZ);  

    /** \brief Write cubit journal file 
     * \param information read from text based input file
     */
    void CreateCubitJournal();

    /** \brief Computes the location of the pincells in the assembly
     * \param pin-number and location of the pincell
     */
    void ComputePinCentroid( int, CMatrix<std::string>, int, int,
			     double&, double&, double&);

  private:
    // iGeom Impl for calling geometry creation/manipulation operations
    iGeom *igeomImpl;

    // number of sides in the geometry
    int m_nSides;
  
    // file Input
    std::ifstream m_FileInput;  
    
    // journal file Output
    std::ofstream m_FileOutput, m_SchemesFile;

    // string for file names
    std::string m_szFile, m_szInFile, m_szGeomFile,m_szJouFile, m_szSchFile;    

    // matrix for holding pincell arrangement
    CMatrix<std::string> m_Assembly; 

    // matrix for holding verts coordinates used in tet-meshing
    CMatrix<double> m_dMTopSurfCoords; 

    // vector for duct specification 
    CMatrix<double> m_dMAssmPitch, m_dMAssmPitchX, m_dMAssmPitchY, m_dMXYAssm, m_dMZAssm;
  
    // vector for material names
    CVector<std::string> m_szAssmMat, m_szAssmMatAlias;
    CMatrix<std::string> m_szMMAlias;  

    // vector holding a pincell
    CVector<CPincell> m_Pincell; 

    // string for geomtype, engine, meshtype
    std::string m_szEngine;
    std::string m_szGeomType;       
    std::string m_szMeshType;
    std::string m_szSideset; 

    // integers for vectors sizes, err etc
    int m_nAssemblyMat, m_nDimensions, m_nPincells , m_nAssmVol, m_nPin, m_nPinX, m_nPinY, err, m_nLineNumber, m_nPlanar, 
      m_nNeumannSetId, m_nMaterialSetId, m_nDuct, m_nDuctNum, m_nJouFlag; 

    // doubles for pincell pitch, pi and mesh sizes resp.
    double m_dPitch, pi, m_dRadialSize, m_dAxialSize, m_dTetMeshSize, m_dMergeTol;      
 
    // igeom related
    std::vector<iBase_EntityHandle> assms, in_pins;
    //iGeom_Instance geom;
    iBase_EntitySetHandle root_set;

 
    // error handlers
    void IOErrorHandler (ErrorStates) const;
    friend class CPincell;

    // parsing related
    std::string szInputString;
    std::string szComment;
    int MAXCHARS;

  };

  inline const char* AssyGen::name()
  {
    return "AssyGen";
  }

  inline bool AssyGen::can_mesh(iBase_EntityType)
  {
    return false;
  }

  inline bool AssyGen::can_mesh(ModelEnt *)
  {
    return true;
  }

  inline const moab::EntityType* AssyGen::mesh_types_arr() const
  {
    return output_types();
  }

} // namespace MeshKit
#endif
