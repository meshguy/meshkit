#ifndef MESHKIT_MESQUITE_OPT_HPP
#define MESHKIT_MESQUITE_OPT_HPP

#include "IQInterface.hpp"
#include "meshkit/MeshOp.hpp"
#include "iBase.h"

namespace MeshKit {

/**\brief Base class for creating MeshOps that call Mesquite smoother
 *
 * This class can be used to either create a MeshOp
 * from any custom Mesquite algorithm or Mesquite wrapper,
 * or as a base class for MeshOp classes that implement
 * a specific Mesquite algorithm.
 */
class MesquiteOpt : public MeshOp
{
  public:
  
    MesquiteOpt( MKCore* core, const MEntVector &me_vec = MEntVector());
    
    /**\brief specify mesquite optimization algorithm to execute */
    void set_mesquite_op( Mesquite::IQInterface* msq_algo ) { msqAlgo = msq_algo; }
    
    /**\brief Get handle of tag used to designate vertices as fixed */
    iBase_TagHandle get_fixed_tag();

    /**\brief Set handle of tag used to designate vertices as fixed */
    void set_fixed_tag( iBase_TagHandle tag );

    /**\brief Set handle of tag used to designate vertices as fixed 
     *
     * Tag will be created if it does not already exists.  If
     * tag already exists then it must be a single integer value.
     */
    void set_fixed_tag( const char* name );
    
    /**\brief smooth only interior of mesh 
     * 
     * This is the default behavior
     */
    void smooth_with_fixed_boundary();
    
    /**\brief smooth interior and exiterior of mesh
     *
     * This method will throw an exception if Mesquite is
     * configured w/out iRel support or the boundary of the
     * mesh is not constrained to geometric entities.
     */
    void smooth_with_free_boundary();
    
    /* virtual functions from MeshOp */
    virtual void setup_this();
    virtual void execute_this();
    virtual const moab::EntityType* mesh_types_arr() const
      { return output_types(); }
    
    /* pseudo-virtual functions from MeshOpProxy */
    static const moab::EntityType* output_types();
    static bool can_mesh( iBase_EntityType dimension );
    static bool can_mesh( ModelEnt* entity );
  
    void print_quality( bool yesno )
      { verboseOutput = yesno; }
  
  protected:
  
    /** Create Mesquite's internal-use tag for it so that we can
     *  create it as a dense tag.
     */
    void create_byte_tag();
  
    /**\brief Set fixed tag on vertices of ModelEnt's mesh entities 
     *\param ent Entity for which to set 'fixed' tag on vertices of
     *           all owned elements.
     *\param value Value to set 'fixed' tag to.  Should be either zero
     *           or one, corresponding to free or fixed, respectively. 
     */
    void set_fixed_tag( ModelEnt* ent, int value );
  
    /**\brief Set fixed tag on boundary of ModelEnt's mesh 
     *\param ent Entity for which to set 'fixed' tag on vertices of boundary
     *\param value Value to set 'fixed' tag to.  Should be either zero
     *           or one, corresponding to free or fixed, respectively. 
     */
    void set_fixed_tag_on_skin( ModelEnt* ent, int value );
  
    /**\brief Verify that all model entities have associated geometry
     */
    bool all_model_ents_have_geom() const;
    
    /**\brief Helper function for free smooth of multiple entities 
     *
     * For proper results when doing a free smooth of multiple model
     * entities, the mesh for all model entities connected by a curve
     * or surface should be smoothed together.  This finds a single
     * connected group of ModelEnts in \c from_this, removes them from
     * \c from_this, and creates a single \c iBase_EntitySetHandle 
     * containing the mesh of all of the ModelEnts in the connected set.
     *
     *\param from_this  Set of ModelEnts in which to search for a group
     *                  of connected ModelEnts.  Connected ModelEnts are
     *                  _removed_ from the passed MEntSet.
     *\param result     An iBase_EntitySetHandle corresponding to a set
     *                  that contains all of the mesh entities of all
     *                  of the connected ModelEnts removed from \c from_this
     *\param dimension  If all connected ModelEnts have the same dimension,
     *                  that dimension. Otherwise iBase_ALL_TYPES.
     *\param created_result_set  True if the function created a new
     *                  iBase_EntitySetHandle.  Will be false if the
     *                  set corresponds to a single ModelEnt that was
     *                  connected to no others, in which case the returned
     *                  set is just that of the ModelEnt.
     */
    void get_adjacent_entity_set( MEntSet& from_this,
                                  iBase_EntitySetHandle& result,
                                  iBase_EntityType& dimension,
                                  bool& created_result_set );   
    
  private:
   
    Mesquite::IQInterface* msqAlgo;
    bool fixedBoundary;
    bool haveFixedTag;
    bool createdByteTag;
    bool verboseOutput;
    iBase_TagHandle fixedTag;
};

}


#endif // MESHKIT_MESQUITE_OPT_HPP
