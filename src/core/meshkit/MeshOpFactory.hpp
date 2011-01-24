#ifndef MESHOPFACTORY_HPP
#define MESHOPFACTORY_HPP

/** \file MeshOpFactory.hpp
 */
#include "meshkit/MeshOp.hpp"
#include "meshkit/ModelEnt.hpp"
#include <vector>

namespace MeshKit {

class MKCore;
    
/** \class MeshOpFactory MeshOpFactory.hpp "meshkit/MeshOpFactory.hpp"
 * \brief A factory for constructing meshing operation (MeshOp) instances
 *
 * MeshOpFactory is one of the only singleton classes in MeshKit.  This factory registers mesh operations
 * and manages their construction.  Valid MeshOp classes define the following:
 * - a name
 * - a list of MOAB entity types or topological dimensions they apply to or generate
 * - a factory function for instantiating objects of this type
 * - a function for evaluating whether the scheme can mesh a given ModelEnt
 *
 * All MeshKit-provided classes derived from MeshOp register themselves with the MeshOpFactory singleton
 * during the static initialization phase of the code.  Applications wanting to define their own MeshOp
 * classes simply register themselves with the MeshOpFactory, either during their static initialization
 * or sometime later.
 *
 * Once a MeshOp has been registered, instances of that type of object can be requested from MeshOpFactory
 * using MeshOpFactory::instance()->get_meshop.  Most of the MeshKit-provided classes derived from MeshOp 
 * cannot be constructed without going through MeshOpFactory, though this is not strictly enforced.  
 * MeshOpFactory keeps track of which MeshOps can generate or operate on a given type of mesh
 * entity, and this information can be queried from applications.
 *
 */

class MeshOpFactory
{
public:

    /** \brief Singleton access function
     *
     * If an MKCore object is passed and the MeshOpFactory has already been instantiated,
     * that argument is compared against the registered one, and an Error thrown if they
     * don't match.
     * \param mk_core MeshKit instance
     * \param create_if_missing If true, an instance will always be returned from this function
     * \return A pointer to the singleton instance
     */
  static MeshOpFactory *instance(MKCore *mk_core = NULL, bool create_if_missing = false);

    /** \brief Destroy the factory instance, and the core if I created it
     * \param dont_destroy_core If true, don't destroy the MKCore, even if I created it
     */
  void destroy_instance(bool dont_destroy_core = false);
  
    /** \brief Change the MKCore instance pointed to by the current factory
     *
     * This function should be used rarely, but is needed to allow multiple MKCore objects
     * to interact with the (singleton) MeshOpFactory.
     * \param mk_core New MKCore for this MeshOpFactory
     */
  void change_core(MKCore *mk_core);
  
    //! Function for constructing MeshOp instances, provided by registrants
  typedef MeshOp* (*meshop_factory_t)(MKCore *, const MEVector &vec);

    //! Function returning whether a given MeshOp can mesh the specified ModelEnt
  typedef bool (*meshop_canmesh_t)(ModelEnt *);

    /** \brief Register a new MeshOp factory
     * \param op_name The name by which this type of MeshOp can be requested
     * \param tp The MOAB entity type operated on by this MeshOp
     * \param meshop The (static) factory function producing instances of this MeshOp type
     * \param canmesh If provided, a static function that returns whether the MeshOp can mesh a given ModelEnt
     * \return Returns true if successful, false otherwise
     */
  bool register_meshop(const char *op_name, moab::EntityType tp, 
                       meshop_factory_t meshop, meshop_canmesh_t canmesh = NULL);
  
    /** \brief Register a new MeshOp factory
     * \param op_name The name by which this type of MeshOp can be requested
     * \param tps The MOAB entity types operated on by this MeshOp
     * \param num_tps Number of entity types in tps
     * \param meshop The (static) factory function producing instances of this MeshOp type
     * \param canmesh If provided, a static function that returns whether the MeshOp can mesh a given ModelEnt
     * \return Returns true if successful, false otherwise
     */
  bool register_meshop(const char *op_name, moab::EntityType *tps, int num_tps,
                       meshop_factory_t meshop, meshop_canmesh_t canmesh = NULL);
  
    /** \struct OpInfo MeshOpFactory.hpp "meshkit/MeshOpFactory.hpp"
     * \brief Struct used to store information about each MeshOp type registered with MeshOpFactory
     */
  struct OpInfo
  {
    std::string opName;
    unsigned int opIndex;
    std::vector<moab::EntityType> opEntTypes;
    meshop_factory_t opFactory;
    meshop_canmesh_t opCanMesh;
  };

    /** \brief A map from MeshOp name to index in the OpInfos vector
     */
  typedef std::map<std::string, unsigned int> OpNameMap;
  
    /** \brief Return the MeshOp type with the given name
     * \param op_name Operation name requested
     * \return OpInfo for the corresponding MeshOp type
     */
  OpInfo meshop_info(const char *op_name);
  
    /** \brief Return the MeshOp index given the name
     * \param op_name Operation name requested
     * \return OpInfo index for the corresponding MeshOp type
     */
  unsigned int meshop_index(const char *op_name);

    /** \brief Make the specified MeshOp name the default for the given dimension(s)
     * 
     * If the specified MeshOp cannot produce entities of the specified dimension, an error is
     * thrown with type MK_BAD_INPUT.
     * \param op_name MeshOp name being set
     * \param dims Bitmask, where 2^x indicates that this MeshOp should be the default for dimension x 
     */
  void set_default_meshop(const char *op_name, unsigned short dims);
  
    /** \brief Make the specified MeshOp name the default for the given dimension(s)
     * 
     * If the specified MeshOp cannot produce entities of the specified dimension, an error is
     * thrown with type MK_BAD_INPUT.
     * \param op_index MeshOp index being set
     * \param dims Bitmask, where 2^x indicates that this MeshOp should be the default for dimension x 
     */
  void set_default_meshop(unsigned short op_index, unsigned short dims);
  
    /** \brief Return MeshOp types that can operate on the specified entity type
     * \param tp Entity type requested
     * \param ops MeshOp types returned
     */
  void meshop_by_type(moab::EntityType tp, std::vector<OpInfo> &ops);
    
    /** \brief Return MeshOp types that can operate on mesh of specified dimension
     * \param dim Entity dimension requested
     * \param ops MeshOp types returned
     */
  void meshop_by_dimension(int dim, std::vector<OpInfo> &ops);
    
    /** \brief Return MeshOp types that can mesh the specified ModelEnt
     * \param ent ModelEnt* requested
     * \param ops MeshOp types returned
     */
  void meshop_by_modelent(ModelEnt * const ent, std::vector<OpInfo> &ops);
    
    /** \brief Construct a new MeshOp of the specified name
     * \param op_name MeshOp name being requested
     * \param me_vec MEVector of entities the operation applies to
     * \return Pointer to new MeshOp constructed
     */
  MeshOp *construct_meshop(std::string op_name, const MEVector &me_vec = MEVector());
    
    /** \brief Construct the default type of MeshOp for the specified dimension
     * \param dim Dimension requested
     * \param me_vec MEVector of entities the operation applies to
     * \return Pointer to new MeshOp constructed
     */
  MeshOp *construct_meshop(int dim, const MEVector &me_vec = MEVector());
    
    //! Return the MKCore
  MKCore *mk_core();
    
private:
    /** \brief Constructor, private because this class is a singleton
     * \param mk_core MKCore instance this factory works with
     * \param create_if_missing If true and mk_core is NULL, will construct a core instance
     */
  MeshOpFactory(MKCore *mk_core, bool create_if_missing);
  
    //! Destructor is private; applications should use destroy_instance() member function;
    //! virtual because applications may want to substitute their own factory
  virtual ~MeshOpFactory();

    //! Static instance
  static MeshOpFactory *instance_;

    //! MKCore instance
  MKCore *mkCore;

    //! Keep track of whether I created the core or not
  bool iCreatedCore;

    /** \brief Master vector of OpInfo's
     */
  std::vector<OpInfo> registeredOps;

    //! Map from MeshOp name to OpInfo
  OpNameMap opNameMap;

    //! OpInfo's by dimension they can produce; first for a given dimension is the default
  std::vector<unsigned short> opsByDim[4];
  
};

inline MeshOpFactory *MeshOpFactory::instance(MKCore *mk_core, bool create_if_missing) 
{
  if (!instance_ && (create_if_missing || !mk_core)) 
    instance_ = new MeshOpFactory(mk_core, create_if_missing);

  else if (instance_ && mk_core && mk_core != instance_->mkCore)
    throw Error(MK_BAD_INPUT, "MKCore passed to MeshOpFactory::instance doesn't match instance's mkCore.");
  
  return instance_;
}

inline void MeshOpFactory::change_core(MKCore *mk_core) 
{
  mkCore = mkCore;
}
    
inline MKCore *MeshOpFactory::mk_core() 
{
  return mkCore;
}

} // namespace MeshKit

#endif

  
