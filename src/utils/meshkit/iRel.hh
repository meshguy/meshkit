#ifndef ITAPS_REL_HH
#define ITAPS_REL_HH

/** \file iRel.hh
 */

#include "iRel.h"
#include "meshkit/iGeom.hh"
#include "meshkit/iMesh.hh"
#include <string.h>
#include <stdlib.h>
#include <vector>
#include <limits>
#include <string>
#include <algorithm>
#include "assert.h"

#define IRELI irelInstance->mInstance

typedef void* iBase_Instance;

/** \class iRel
 * \brief C++ interface to ITAPS iRel interface
 *
 * This class is a simple wrapper for the ITAPS %iRel interface.  The primary benefit to using this class
 * instead of %iRel directly is that lists of handles are passed as std::vectors instead of pointers to
 * handle arrays.  This file includes both declaration and definition of all iRel class functions, i.e.
 * all functions are inlined.  The class can be constructed and destructed in the standard C++ way; the
 * implementation of those functions call into the standard %iRel C functions newRel and dtor.
 *
 * For complete documentation of these functions, see the 
 * [http://trac.mcs.anl.gov/projects/ITAPS/browser/Lasso/trunk/iRel.h iRel header] in the Lasso source
 * (for now).
 */
class iRel {
public:
    
  typedef iBase_ErrorType Error;

  iRel( const char* options = 0 );
  iRel( iRel_Instance instance );
  
  ~iRel();
    
  enum IfaceType 
  {IGEOM_IFACE = 0, 
   IMESH_IFACE, 
   IFIELD_IFACE, 
   IREL_IFACE};

  enum RelationType 
  {ENTITY = 0, 
   SET, 
   BOTH};

/** \class PairHandle iRel.hh "iRel.hh"
 * \brief Class for storing, querying and modifying relations data.
 *
 * This class encapsulates most of the functions for querying relations in the C++ version of %iRel.
 */
  class PairHandle 
  {
  public:
    PairHandle(iRel *instance,
               iBase_Instance iface1 = NULL, RelationType rtype1 = ENTITY, IfaceType itype1 = IREL_IFACE,
               iBase_Instance iface2 = NULL, RelationType rtype2 = ENTITY, IfaceType itype2 = IREL_IFACE);
    
    PairHandle(iRel *instance, iRel_PairHandle ph);
    
    ~PairHandle();

    friend class iRel;

    iBase_Instance get_iface(int i);

    IfaceType get_iface_type(int i);
    RelationType get_rel_type(int i);
    
    Error changeType (RelationType ent_or_set1, 
                             RelationType ent_or_set2);

    Error setEntEntRelation (
        /*in*/ iBase_EntityHandle ent1,
        /*in*/ iBase_EntityHandle ent2);
    Error setEntSetRelation (
        /*in*/ iBase_EntityHandle ent1,
        /*in*/ iBase_EntitySetHandle entset2);
    Error setSetEntRelation (
        /*in*/ iBase_EntitySetHandle entset1,
        /*in*/ iBase_EntityHandle ent2);
    Error setSetSetRelation (
        /*in*/ iBase_EntitySetHandle entset1,
        /*in*/ iBase_EntitySetHandle entset2);

    Error setEntEntArrRelation (
        /*in*/ iBase_EntityHandle ent1,
        /*in*/ bool switch_order,
        /*in*/ iBase_EntityHandle *ent_array_2,
        /*in*/ int num_entities);
    Error setSetEntArrRelation (
        /*in*/ iBase_EntitySetHandle entset1,
        /*in*/ bool switch_order,
        /*in*/ iBase_EntityHandle *ent_array_2,
        /*in*/ int num_entities);
    Error setEntSetArrRelation (
        /*in*/ iBase_EntityHandle ent1,
        /*in*/ bool switch_order,
        /*in*/ iBase_EntitySetHandle *entset_array_2,
        /*in*/ int num_entities);
    Error setSetSetArrRelation (
        /*in*/ iBase_EntitySetHandle entset1,
        /*in*/ bool switch_order,
        /*in*/ iBase_EntitySetHandle *entset_array_2,
        /*in*/ int num_entities);

    Error setEntArrEntArrRelation (
        /*in*/ iBase_EntityHandle *ent_array_1,
        /*in*/ int num_ent1,
        /*in*/ iBase_EntityHandle *ent_array_2,
        /*in*/ int num_ent2);
    Error setSetArrEntArrRelation (
        /*in*/ iBase_EntitySetHandle *entset_array_1,
        /*in*/ int num_ent1,
        /*in*/ iBase_EntityHandle *ent_array_2,
        /*in*/ int num_ent2);
    Error setEntArrSetArrRelation (
        /*in*/ iBase_EntityHandle *ent_array_1,
        /*in*/ int num_ent1,
        /*in*/ iBase_EntitySetHandle *entset_array_2,
        /*in*/ int num_ent2);
    Error setSetArrSetArrRelation (
        /*in*/ iBase_EntitySetHandle *entset_array_1,
        /*in*/ int num_ent1,
        /*in*/ iBase_EntitySetHandle *entset_array_2,
        /*in*/ int num_ent2);

    Error getEntEntRelation (
        /*in*/ iBase_EntityHandle ent1,
        /*in*/ bool switch_order,
        /*out*/ iBase_EntityHandle &ent2);
    Error getEntSetRelation (
        /*in*/ iBase_EntityHandle ent1,
        /*in*/ bool switch_order,
        /*out*/ iBase_EntitySetHandle &entset2);
    Error getSetEntRelation (
        /*in*/ iBase_EntitySetHandle entset1,
        /*in*/ bool switch_order,
        /*out*/ iBase_EntityHandle &ent2);
    Error getSetSetRelation (
        /*in*/ iBase_EntitySetHandle entset1,
        /*in*/ bool switch_order,
        /*out*/ iBase_EntitySetHandle &entset2);
    Error getEntSetIterRelation (
        /*in*/ iBase_EntityHandle ent1,
        /*in*/ bool switch_order,
        /*out*/ iBase_EntityIterator &entset2);

    Error getEntEntArrRelation (
        /*in*/ iBase_EntityHandle ent1,
        /*in*/ bool switch_order,
        /*inout*/ std::vector<iBase_EntityHandle> &ent_array_2);
    Error getSetEntArrRelation (
        /*in*/ iBase_EntitySetHandle entset1,
        /*in*/ bool switch_order,
        /*inout*/ std::vector<iBase_EntityHandle> &ent_array_2);

    Error getEntArrEntArrRelation (
        /*in*/ iBase_EntityHandle *ent_array_1,
        /*in*/ int ent_array_1_size,
        /*in*/ bool switch_order,
        /*inout*/ std::vector<iBase_EntityHandle> &ent_array_2,
        /*inout*/ std::vector<int> &offset);
    Error getEntArrSetArrRelation (
        /*in*/ iBase_EntityHandle *ent_array_1,
        /*in*/ int ent_array_1_size,
        /*in*/ bool switch_order,
        /*inout*/ std::vector<iBase_EntitySetHandle> &entset_array_2);
    Error getSetArrEntArrRelation (
        /*in*/ iBase_EntitySetHandle *entset_array_1,
        /*in*/ int entset_array_1_size,
        /*in*/ bool switch_order,
        /*inout*/ std::vector<iBase_EntityHandle> &ent_array_2,
        /*inout*/ std::vector<int> &offset);
    Error getSetArrSetArrRelation (
        /*in*/ iBase_EntitySetHandle *entset_array_1,
        /*in*/ int entset_array_1_size,
        /*in*/ bool switch_order,
        /*inout*/ std::vector<iBase_EntitySetHandle> &entset_array_2);
    Error getEntArrSetIterArrRelation (
        /*in*/ iBase_EntityHandle *ent_array_1,
        /*in*/ int ent_array_1_size,
        /*in*/ bool switch_order,
        /*inout*/ std::vector<iBase_EntityIterator> &entiter);

    Error inferAllRelations ();

    Error inferAllRelationsAndType ();

    Error inferEntRelations (
        /*in*/ iBase_EntityHandle entity,
        /*in*/ int iface_no);
    Error inferSetRelations (
        /*in*/ iBase_EntitySetHandle entity_set,
        /*in*/ int iface_no);

    Error inferEntArrRelations (
        /*in*/ iBase_EntityHandle *entities,
        /*in*/ int entities_size,
        /*in*/ int iface_no);
    Error inferSetArrRelations (
        /*in*/ iBase_EntitySetHandle *entity_sets,
        /*in*/ int entities_size,
        /*in*/ int iface_no);

  private:
    iRel *irelInstance;
    iBase_Instance iFaces[2];
    IfaceType iFaceTypes[2];
    RelationType relTypes[2];
    
    iRel_PairHandle irelPair;

    bool pairOwner;
  };
  
  Error getErrorType();
    
  std::string getDescription();
    
  Error createPair (
    /*in*/ iBase_Instance iface1,
    /*in*/ const RelationType ent_or_set1,
    /*in*/ const IfaceType iface_type1,
    /*in*/ iBase_Instance iface2,
    /*in*/ const RelationType ent_or_set2,
    /*in*/ const IfaceType iface_type2,
    /*out*/ PairHandle *&pair);

  Error findPairs (
    /*in*/ iBase_Instance iface,
    /*inout*/ std::vector<PairHandle *> &pairs);

protected:

  Error add_pair(PairHandle *pair);

  Error remove_pair(PairHandle *pair);
  
private:
  iRel_Instance mInstance;
  
  bool iRelInstanceOwner;

  std::vector<PairHandle*> relPairs;
  
      // prohibit copying
  iRel( const iRel& ) {}
  void operator=(const iRel&) {}
};

inline 
iRel::iRel( const char* options )
  : iRelInstanceOwner(true)
{
  int err, len = options ? strlen(options) : 0;
  iRel_newRel( options, &mInstance, &err, len );
  if (iBase_SUCCESS != err) {
    mInstance = 0;
    iRelInstanceOwner = false;
  }
}

inline 
iRel::iRel( iRel_Instance instance )
  : iRelInstanceOwner(false)
{
  mInstance = instance;
}

inline iRel::~iRel()
{
  if (iRelInstanceOwner) {
    int err;
    iRel_dtor( mInstance, &err );
  }
}

inline iRel::Error iRel::getErrorType()
{
  int err1, err2;
  iRel_getErrorType(mInstance, &err1, &err2 );
  if (iBase_SUCCESS != err2)
    return (Error)err2;
  else
    return (Error)err1;
}
    
inline std::string iRel::getDescription() 
{
  std::vector<char> buffer(1024);
  int err;
  iRel_getDescription(mInstance, &buffer[0], &err, buffer.size() );
  if (iBase_SUCCESS != err)
    return std::string();
  else
    return std::string(&buffer[0]);
}
    
inline iRel::Error iRel::createPair (
    /*in*/ iBase_Instance iface1,
    /*in*/ const RelationType ent_or_set1,
    /*in*/ const IfaceType iface_type1,
    /*in*/ iBase_Instance iface2,
    /*in*/ const RelationType ent_or_set2,
    /*in*/ const IfaceType iface_type2,
    /*out*/ PairHandle *&pair) 
{
  pair = new PairHandle(this, iface1, ent_or_set1, iface_type1, iface2, ent_or_set2, iface_type2);
  return (Error)iBase_SUCCESS;
}

inline iRel::Error iRel::findPairs (
    /*in*/ iBase_Instance iface,
    /*inout*/ std::vector<PairHandle *> &pairs) 
{
  for (std::vector<PairHandle*>::iterator vit = relPairs.begin(); vit != relPairs.end(); vit++) {
    if ((*vit)->iFaces[0] == iface || (*vit)->iFaces[1] == iface)
      pairs.push_back(*vit);
  }
  
  return (Error)iBase_SUCCESS;
}

inline iRel::Error iRel::add_pair(PairHandle *pair) 
{
  assert(std::find(relPairs.begin(), relPairs.end(), pair) == relPairs.end());
  relPairs.push_back(pair);
  return (Error)iBase_SUCCESS;
}

inline iRel::Error iRel::remove_pair(PairHandle *pair) 
{
  std::vector<PairHandle*>::iterator vit = std::find(relPairs.begin(), relPairs.end(), pair);
  if (vit == relPairs.end()) return (Error)iBase_FAILURE;
  else relPairs.erase(vit);
  return (Error)iBase_SUCCESS;
}
  
inline iRel::PairHandle::PairHandle(iRel *instance,
                                    iBase_Instance iface1, RelationType rtype1, IfaceType itype1,
                                    iBase_Instance iface2, RelationType rtype2, IfaceType itype2)
  : irelInstance(instance), 
    pairOwner(true)
{
  iFaces[0] = iface1;
  iFaces[1] = iface2;
  iFaceTypes[0] = itype1;
  iFaceTypes[1] = itype2;
  relTypes[0] = rtype1;
  relTypes[1] = rtype2;
    
  int err;
  iRel_createPair(IRELI, iface1, rtype1, itype1, iface2, rtype2, itype2, &irelPair, &err);
  instance->add_pair(this);
}
    
inline iRel::PairHandle::PairHandle(iRel *instance, iRel_PairHandle ph)
        : irelInstance(instance), pairOwner(false)
{
  int err;
  int rtype1, rtype2, itype1, itype2;
  iRel_getPairInfo(IRELI, ph, iFaces, &rtype1, &itype1,
                   iFaces+1, &rtype2, &itype2, &err);
  if (iBase_SUCCESS == err) {
    iFaceTypes[0] = (IfaceType)itype1;
    iFaceTypes[1] = (IfaceType)itype2;
    relTypes[0] = (RelationType)rtype1;
    relTypes[1] = (RelationType)rtype2;
  }

  instance->add_pair(this);
}
    
inline iRel::PairHandle::~PairHandle() 
{
  irelInstance->remove_pair(this);

  if (pairOwner) {
    int err;
    iRel_destroyPair(IRELI, irelPair, &err);
  }
}

inline iBase_Instance iRel::PairHandle::get_iface(int i) 
{
  assert(i == 0 || i == 1);
  return iFaces[i];
}

inline iRel::IfaceType iRel::PairHandle::get_iface_type(int i) 
{
  assert(i == 0 || i == 1);
  return iFaceTypes[i];
}

inline iRel::RelationType iRel::PairHandle::get_rel_type(int i)
{
  assert(i == 0 || i == 1);
  return relTypes[i];
}
    
inline iRel::Error iRel::PairHandle::changeType (RelationType ent_or_set1, 
                                           RelationType ent_or_set2) 
{
  int err;
  iRel_changePairType(IRELI, irelPair, ent_or_set1, ent_or_set2, &err);
  return (Error)err;
}

inline iRel::Error iRel::PairHandle::setEntEntRelation (
    /*in*/ iBase_EntityHandle ent1,
    /*in*/ iBase_EntityHandle ent2)
{
  int err;
  iRel_setEntEntRelation (IRELI, irelPair, ent1, ent2, &err);
  return (Error)err;
}

inline iRel::Error iRel::PairHandle::setEntSetRelation (
    /*in*/ iBase_EntityHandle ent1,
    /*in*/ iBase_EntitySetHandle entset2)
{
  int err;
  iRel_setEntSetRelation (IRELI, irelPair, ent1, entset2, &err);
  return (Error)err;

}
inline iRel::Error iRel::PairHandle::setSetEntRelation (
    /*in*/ iBase_EntitySetHandle entset1,
    /*in*/ iBase_EntityHandle ent2)
{
  int err;
  iRel_setSetEntRelation (IRELI, irelPair, entset1, ent2, &err);
  return (Error)err;
}

inline iRel::Error iRel::PairHandle::setSetSetRelation (
    /*in*/ iBase_EntitySetHandle entset1,
    /*in*/ iBase_EntitySetHandle entset2)
{
  int err;
  iRel_setSetSetRelation (IRELI, irelPair, entset1, entset2, &err);
  return (Error)err;
}


inline iRel::Error iRel::PairHandle::setEntEntArrRelation (
    /*in*/ iBase_EntityHandle ent1,
    /*in*/ bool switch_order,
    /*in*/ iBase_EntityHandle *ent_array_2,
    /*in*/ int num_entities)
{
  int err;
  iRel_setEntEntArrRelation (IRELI, irelPair, ent1, switch_order, ent_array_2, num_entities, &err);
  return (Error)err;
}

inline iRel::Error iRel::PairHandle::setSetEntArrRelation (
    /*in*/ iBase_EntitySetHandle entset1,
    /*in*/ bool switch_order,
    /*in*/ iBase_EntityHandle *ent_array_2,
    /*in*/ int num_entities)
{
  int err;
  iRel_setSetEntArrRelation (IRELI, irelPair, entset1, switch_order, ent_array_2, num_entities, &err);
  return (Error)err;
}

inline iRel::Error iRel::PairHandle::setEntSetArrRelation (
    /*in*/ iBase_EntityHandle ent1,
    /*in*/ bool switch_order,
    /*in*/ iBase_EntitySetHandle *entset_array_2,
    /*in*/ int num_entities)
{
  int err;
  iRel_setEntSetArrRelation (IRELI, irelPair, ent1, switch_order, entset_array_2, num_entities, &err);
  return (Error)err;
}

inline iRel::Error iRel::PairHandle::setSetSetArrRelation (
    /*in*/ iBase_EntitySetHandle entset1,
    /*in*/ bool switch_order,
    /*in*/ iBase_EntitySetHandle *entset_array_2,
    /*in*/ int num_entities)
{
  int err;
  iRel_setSetSetArrRelation (IRELI, irelPair, entset1, switch_order, entset_array_2, num_entities, &err);
  return (Error)err;
}


inline iRel::Error iRel::PairHandle::setEntArrEntArrRelation (
    /*in*/ iBase_EntityHandle *ent_array_1,
    /*in*/ int num_ent1,
    /*in*/ iBase_EntityHandle *ent_array_2,
    /*in*/ int num_ent2)
{
  int err;
  iRel_setEntArrEntArrRelation (IRELI, irelPair, ent_array_1, num_ent1, ent_array_2, num_ent2, &err);
  return (Error)err;
}

inline iRel::Error iRel::PairHandle::setSetArrEntArrRelation (
    /*in*/ iBase_EntitySetHandle *entset_array_1,
    /*in*/ int num_ent1,
    /*in*/ iBase_EntityHandle *ent_array_2,
    /*in*/ int num_ent2)
{
  int err;
  iRel_setSetArrEntArrRelation (IRELI, irelPair, entset_array_1, num_ent1, ent_array_2, num_ent2, &err);
  return (Error)err;
}

inline iRel::Error iRel::PairHandle::setEntArrSetArrRelation (
    /*in*/ iBase_EntityHandle *ent_array_1,
    /*in*/ int num_ent1,
    /*in*/ iBase_EntitySetHandle *entset_array_2,
    /*in*/ int num_ent2)
{
  int err;
  iRel_setEntArrSetArrRelation (IRELI, irelPair, ent_array_1, num_ent1, entset_array_2, num_ent2, &err);
  return (Error)err;
}

inline iRel::Error iRel::PairHandle::setSetArrSetArrRelation (
    /*in*/ iBase_EntitySetHandle *entset_array_1,
    /*in*/ int num_ent1,
    /*in*/ iBase_EntitySetHandle *entset_array_2,
    /*in*/ int num_ent2)
{
  int err;
  iRel_setSetArrSetArrRelation (IRELI, irelPair, 
                                /*in*/ entset_array_1, num_ent1, entset_array_2, num_ent2, &err);
  return (Error)err;
}


inline iRel::Error iRel::PairHandle::getEntEntRelation (
    /*in*/ iBase_EntityHandle ent1,
    /*in*/ bool switch_order,
    /*out*/ iBase_EntityHandle &ent2)
{
  int err;
  iRel_getEntEntRelation (IRELI, irelPair, ent1, switch_order, &ent2, &err);
  return (Error)err;
}

inline iRel::Error iRel::PairHandle::getEntSetRelation (
    /*in*/ iBase_EntityHandle ent1,
    /*in*/ bool switch_order,
    /*out*/ iBase_EntitySetHandle &entset2)
{
  int err;
  iRel_getEntSetRelation (IRELI, irelPair, ent1, switch_order, &entset2, &err);
  return (Error)err;
}

inline iRel::Error iRel::PairHandle::getSetEntRelation (
    /*in*/ iBase_EntitySetHandle entset1,
    /*in*/ bool switch_order,
    /*out*/ iBase_EntityHandle &ent2)
{
  int err;
  iRel_getSetEntRelation (IRELI, irelPair, entset1, switch_order, &ent2, &err);
  return (Error)err;
}

inline iRel::Error iRel::PairHandle::getSetSetRelation (
    /*in*/ iBase_EntitySetHandle entset1,
    /*in*/ bool switch_order,
    /*out*/ iBase_EntitySetHandle &entset2)
{
  int err;
  iRel_getSetSetRelation (IRELI, irelPair, entset1, switch_order, &entset2, &err);
  return (Error)err;
}

inline iRel::Error iRel::PairHandle::getEntSetIterRelation (
    /*in*/ iBase_EntityHandle ent1,
    /*in*/ bool switch_order,
    /*out*/ iBase_EntityIterator &entset2)
{
  int err;
  iRel_getEntSetIterRelation (IRELI, irelPair, ent1, switch_order, &entset2, &err);
  return (Error)err;
}

inline iRel::Error iRel::PairHandle::getEntEntArrRelation (
    /*in*/ iBase_EntityHandle ent1,
    /*in*/ bool switch_order,
    /*inout*/ std::vector<iBase_EntityHandle> &ent_array_2)
{
  int err;
  int dum1 = 0, dum2;
  iBase_EntityHandle *dum_ptr = NULL;
  
  iRel_getEntEntArrRelation (IRELI, irelPair, ent1, switch_order,
                             &dum_ptr, &dum1, &dum2, &err);
  if (iBase_SUCCESS == err) {
    ent_array_2.resize(dum1);
    memcpy(&ent_array_2[0], dum_ptr, dum1*sizeof(iBase_EntityHandle));
  }

  return (Error)err;

}
inline iRel::Error iRel::PairHandle::getSetEntArrRelation (
    /*in*/ iBase_EntitySetHandle entset1,
    /*in*/ bool switch_order,
    /*inout*/ std::vector<iBase_EntityHandle> &ent_array_2)
{
  int err;
  int dum1 = 0, dum2;
  iBase_EntityHandle *dum_ptr = NULL;
  iRel_getSetEntArrRelation (IRELI, irelPair, entset1, switch_order, &dum_ptr, &dum1, &dum2, &err);
  if (iBase_SUCCESS == err) {
    ent_array_2.resize(dum1);
    memcpy(&ent_array_2[0], dum_ptr, dum1*sizeof(iBase_EntityHandle));
  }
  return (Error)err;

}

inline iRel::Error iRel::PairHandle::getEntArrEntArrRelation (
    /*in*/ iBase_EntityHandle *ent_array_1,
    /*in*/ int ent_array_1_size,
    /*in*/ bool switch_order,
    /*inout*/ std::vector<iBase_EntityHandle> &ent_array_2,
    /*inout*/ std::vector<int> &offset)
{
  int err;
  int dum1 = 0, dum2, *dum_ptr2 = NULL, dum3, dum4;
  iBase_EntityHandle *dum_ptr = NULL;
  iRel_getEntArrEntArrRelation (IRELI, irelPair, ent_array_1, ent_array_1_size, switch_order,
                                &dum_ptr, &dum1, &dum2, 
                                &dum_ptr2, &dum3, &dum4, &err);
  if (iBase_SUCCESS == err) {
    ent_array_2.resize(dum1);
    memcpy(&ent_array_2[0], dum_ptr, dum1*sizeof(iBase_EntityHandle));
    offset.resize(dum3);
    memcpy(&offset[0], dum_ptr2, dum3*sizeof(int));
  }
  return (Error)err;

}
inline iRel::Error iRel::PairHandle::getEntArrSetArrRelation (
    /*in*/ iBase_EntityHandle *ent_array_1,
    /*in*/ int ent_array_1_size,
    /*in*/ bool switch_order,
    /*inout*/ std::vector<iBase_EntitySetHandle> &entset_array_2)
{
  int err;
  int dum1 = 0, dum2;
  iBase_EntitySetHandle *dum_ptr = NULL;
  iRel_getEntArrSetArrRelation (IRELI, irelPair, ent_array_1, ent_array_1_size, switch_order,
                                &dum_ptr, &dum1, &dum2, &err);
  if (iBase_SUCCESS == err) {
    entset_array_2.resize(dum1);
    memcpy(&entset_array_2[0], dum_ptr, dum1*sizeof(iBase_EntitySetHandle));
  }
  return (Error)err;

}
inline iRel::Error iRel::PairHandle::getSetArrEntArrRelation (
    /*in*/ iBase_EntitySetHandle *entset_array_1,
    /*in*/ int entset_array_1_size,
    /*in*/ bool switch_order,
    /*inout*/ std::vector<iBase_EntityHandle> &ent_array_2,
    /*inout*/ std::vector<int> &offset)
{
  int err;
  int dum1 = 0, dum2, *dum_ptr2 = NULL, dum3, dum4;
  iBase_EntityHandle *dum_ptr = NULL;
  iRel_getSetArrEntArrRelation (IRELI, irelPair, entset_array_1, entset_array_1_size, switch_order,
                                &dum_ptr, &dum1, &dum2, 
                                &dum_ptr2, &dum3, &dum4, &err);
  if (iBase_SUCCESS == err) {
    ent_array_2.resize(dum1);
    memcpy(&ent_array_2[0], dum_ptr, dum1*sizeof(iBase_EntityHandle));
    offset.resize(dum3);
    memcpy(&offset[0], dum_ptr2, dum3*sizeof(int));
  }
  return (Error)err;

}
inline iRel::Error iRel::PairHandle::getSetArrSetArrRelation (
    /*in*/ iBase_EntitySetHandle *entset_array_1,
    /*in*/ int entset_array_1_size,
    /*in*/ bool switch_order,
    /*inout*/ std::vector<iBase_EntitySetHandle> &entset_array_2)
{
  int err;
  int dum1 = 0, dum2;
  iBase_EntitySetHandle *dum_ptr = NULL;
  iRel_getSetArrSetArrRelation (IRELI, irelPair, entset_array_1, entset_array_1_size, switch_order,
                                &dum_ptr, &dum1, &dum2, &err);
  if (iBase_SUCCESS == err) {
    entset_array_2.resize(dum1);
    memcpy(&entset_array_2[0], dum_ptr, dum1*sizeof(iBase_EntitySetHandle));
  }
  return (Error)err;

}
inline iRel::Error iRel::PairHandle::getEntArrSetIterArrRelation (
    /*in*/ iBase_EntityHandle *ent_array_1,
    /*in*/ int ent_array_1_size,
    /*in*/ bool switch_order,
    /*inout*/ std::vector<iBase_EntityIterator> &entiter)
{
  int err;
  int dum1 = 0, dum2;
  iBase_EntityIterator *dum_ptr = NULL;
  iRel_getEntArrSetIterArrRelation (IRELI, irelPair, ent_array_1, ent_array_1_size, switch_order,
                                    &dum_ptr, &dum1, &dum2, &err);
  if (iBase_SUCCESS == err) {
    entiter.resize(dum1);
    memcpy(&entiter[0], dum_ptr, dum1*sizeof(iBase_EntityIterator));
  }
  return (Error)err;

}

inline iRel::Error iRel::PairHandle::inferAllRelations ()
{
  int err;
  iRel_inferAllRelations (IRELI, irelPair, &err);
  return (Error)err;
}

inline iRel::Error iRel::PairHandle::inferAllRelationsAndType ()
{
  int err;
  iRel_inferAllRelationsAndType (IRELI, &irelPair, &err);
  return (Error)err;
}


inline iRel::Error iRel::PairHandle::inferEntRelations (
    /*in*/ iBase_EntityHandle entity,
    /*in*/ int iface_no)
{
  int err;
  iRel_inferEntRelations (IRELI, irelPair, entity, iface_no, &err);
  return (Error)err;
}

inline iRel::Error iRel::PairHandle::inferSetRelations (
    /*in*/ iBase_EntitySetHandle entity_set,
    /*in*/ int iface_no)
{
  int err;
  iRel_inferSetRelations (IRELI, irelPair, entity_set, iface_no, &err);
  return (Error)err;
}

inline iRel::Error iRel::PairHandle::inferEntArrRelations (
    /*in*/ iBase_EntityHandle *entities,
    /*in*/ int entities_size,
    /*in*/ int iface_no)
{
  int err;
  iRel_inferEntArrRelations (IRELI, irelPair, entities, entities_size, iface_no, &err);
  return (Error)err;
}

inline iRel::Error iRel::PairHandle::inferSetArrRelations (
        /*in*/ iBase_EntitySetHandle *entity_sets,
        /*in*/ int entities_size,
        /*in*/ int iface_no)
{
  int err;
  iRel_inferSetArrRelations (IRELI, irelPair, entity_sets, entities_size, iface_no, &err);
  return (Error)err;
}

#endif /* #ifndef ITAPS_REL_HH */

