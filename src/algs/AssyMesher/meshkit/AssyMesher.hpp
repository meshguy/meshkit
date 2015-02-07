#ifndef MESHKIT_ASSYMESHER_HPP
#define MESHKIT_ASSYMESHER_HPP

#include <cassert>
#include <string>
#include <vector>
#include <set>
#include <iomanip>

#include "meshkit/Types.hpp"
#include "meshkit/Error.hpp"
#include "meshkit/MeshScheme.hpp"
#include "meshkit/ModelEnt.hpp"

#include "meshkit/LocalTag.hpp"
#include "meshkit/Matrix.hpp"

#include "meshkit/iMesh.hpp"
#include "meshkit/iGeom.hpp"
#include "meshkit/mstream.hpp"

#include "meshkit/SimpleArray.hpp"

#include "meshkit/vectortemplate.hpp"
#include "meshkit/matrixtemplate.hpp"
#include "meshkit/pincell.hpp"
#include "meshkit/parser.hpp"
#include "meshkit/clock.hpp"

namespace MeshKit {
#define DEFAULT_TEST_AM  "assygen_default"


class MKCore;

class AssyMesher : public MeshScheme
{
public:
  /* \brief Constructor
   *
   * Create a new AssyMesher instance
   * \param impl the iGeom instance handle for the Geom
   */
  AssyMesher(MKCore *mkcore, const MEntVector &me_vec);

  /* \brief Destructor
   */
  virtual ~AssyMesher();

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
   *  command line args and testdir for default test case
   */
  void PrepareIO (int argc, char *argv[], std::string TestDir);

  /** \brief Read pincell data
   *  input file
   */
  void ReadPinCellData( int i);


private:
  //! iGeom Impl for calling geometry creation/manipulation operations
  iGeom *igeom;

  //! iMesh Impl for calling mesh creation/manipulation operations
  iMesh *imesh;

  //! MOAB Impl for calling mesh creation/manipulation operations
  moab::Interface *mb;

  // igeom related
  std::vector<iBase_EntityHandle> assms, in_pins;

  // number of sides in the geometry
  int m_nSides;

  // !! file Input
  std::ifstream m_FileInput;
  mstream m_LogFile;
  std::string szInputString;
  std::string szComment;
  int MAXCHARS;

  // ! variables to parse
  std::string m_InputFile, m_GeomFile, m_MeshFile, m_OutFile, m_LogName, m_MeshType;
  std::string m_Card;

  // ! error handlers
  enum ErrorStates {PINCELLS, INVALIDINPUT, EMAT, EGEOMTYPE, EGEOMENGINE, ENEGATIVE, EALIAS, EPIN};
  void IOErrorHandler (ErrorStates) const;


  //// from assygen
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
  int m_nAssemblyMat, m_nDimensions, m_nPincells, m_nLineNumber, m_nPlanar,
    m_nNeumannSetId, m_nMaterialSetId, m_nDuct, m_nDuctNum, m_nJouFlag;

  // doubles for pincell pitch, pi and mesh sizes resp.
  double m_dPitch, pi, m_dRadialSize, m_dAxialSize, m_dTetMeshSize, m_dMergeTol;



};

inline const char* AssyMesher::name()
{
  return "AssyMesher";
}

inline bool AssyMesher::can_mesh(iBase_EntityType)
{
  // Given just a dimension, AssyMesher can't do anything since it doesn't know
  // what to copy.
  return false;
}

inline bool AssyMesher::can_mesh(ModelEnt *)
{
  return true;
}

inline const moab::EntityType* AssyMesher::mesh_types_arr() const
{
  return output_types();
}

} // namespace MeshKit
#endif
