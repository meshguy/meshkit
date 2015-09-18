#include "meshkit/MesquiteOpt.hpp"
#include "meshkit/ModelEnt.hpp"
#include "meshkit/iMesh.hpp"
#include "moab/Skinner.hpp"
#include "MsqIMesh.hpp"
#include "QualityAssessor.hpp"
#include "MsqError.hpp"
#ifdef MSQIGEOM
# include "MsqIGeom.hpp"
#endif
#include "meshkit/FreeSmoothDomain.hpp"

#include "mesquite_version.h"
#if MSQ_VERSION_MAJOR > 2 || (MSQ_VERSION_MAJOR == 2 && MSQ_VERSION_MINOR > 1)
#include "ShapeImprover.hpp"
typedef Mesquite::ShapeImprover DefaultAlgorithm;
#else
#include "ShapeImprovementWrapper.hpp"
typedef Mesquite::ShapeImprovementWrapper DefaultAlgorithm;
#endif

    
#define MSQERRCHK(err)                                 \
  do {                                                 \
    if ((err)) {                                       \
      throw Error(MK_FAILURE,"%s:%d: Mesquite error: %s",     \
          __FILE__, __LINE__, (err).error_message() ); \
    }                                                  \
  } while(false)


using namespace Mesquite;

namespace MeshKit { 

MesquiteOpt::MesquiteOpt( MKCore* core, const MEntVector &me_vec )
  : MeshOp(core, me_vec),
    msqAlgo(0),
    fixedBoundary(true),
    haveFixedTag(false),
    createdByteTag(false),
    verboseOutput(false)
{}

const moab::EntityType* MesquiteOpt::output_types() 
{
  const static moab::EntityType end = moab::MBMAXTYPE;
  return &end;
}

bool MesquiteOpt::can_mesh( iBase_EntityType dim )
{
  return dim >= 2;
}

bool MesquiteOpt::can_mesh( ModelEnt* ent )
{
  return can_mesh( (iBase_EntityType)ent->dimension() );
}

iBase_TagHandle MesquiteOpt::get_fixed_tag()
{
  if (!haveFixedTag)
    set_fixed_tag( "fixed" );
  return fixedTag;
}

void MesquiteOpt::set_fixed_tag( iBase_TagHandle tag )
{
  int size;
  iMesh::TagValueType type;
  mk_core()->imesh_instance()->getTagSizeValues( tag, size );
  mk_core()->imesh_instance()->getTagType( tag, type );
  if (size != 1 || type != iBase_INTEGER) 
    throw Error(MK_WRONG_DIMENSION,"Tag as invalid type or size");
  haveFixedTag = true;
  fixedTag = tag;
}

void MesquiteOpt::set_fixed_tag( const char* name )
{
  iBase_TagHandle tag;
  iBase_ErrorType err = mk_core()->imesh_instance()->getTagHandle( name, tag );
  if (iBase_SUCCESS == err) {
    set_fixed_tag(tag);
  }
  else {
      // create though moab so we can make it a dense tag
    moab::Tag t;
    moab::ErrorCode rval;
    int zero = 0;
    rval = mk_core()->moab_instance()->tag_get_handle( name,
                                                   1,
                                                   moab::MB_TYPE_INTEGER,
                                                   t,
                                                   moab::MB_TAG_DENSE|moab::MB_TAG_CREAT,
                                                   &zero );
    MBERRCHK(rval,mk_core()->moab_instance());
    
    err = mk_core()->imesh_instance()->getTagHandle( name, tag );
    IBERRCHK(err,"Tag created via Moab not accessable via iMesh");
    set_fixed_tag( tag );
  }
}

void MesquiteOpt::create_byte_tag()
{
  if (!createdByteTag) {
        // Create though moab so we can make it a dense tag
        // Don't check error code.  We don't care if it already exists
        // and futher, we really aren't responsiblefor creating it anyway.
        // We do it here only so we can make it a dense tag.
    moab::Tag t;
    unsigned char zero = 0;
    mk_core()->moab_instance()->tag_get_handle( Mesquite::VERTEX_BYTE_TAG_NAME,
                                            1,
                                            moab::MB_TYPE_INTEGER,
                                            t,
                                            moab::MB_TAG_DENSE|moab::MB_TAG_CREAT,
                                            &zero );
    createdByteTag = true;
  }
}

void MesquiteOpt::smooth_with_fixed_boundary()
{ 
  fixedBoundary = true;
}

void MesquiteOpt::smooth_with_free_boundary()
{
#ifndef MSQIGEOM
  throw Error(MK_NOT_IMPLEMENTED,"Cannot do free smooth w/out Mesquite iRel support");
#else
  fixedBoundary = false;
#endif
}

void MesquiteOpt::setup_this()
{
  get_fixed_tag();
  create_byte_tag();
}

void MesquiteOpt::set_fixed_tag( ModelEnt* ent, int value )
{
  std::vector<int> fixed_vals;
  iBase_TagHandle fixed = get_fixed_tag();
  iMesh* mesh = mk_core()->imesh_instance();
  std::vector<iBase_EntityHandle> ents, verts;
  iMesh::Error err;

  ents.clear();
  err = mesh->getEntities( (iBase_EntitySetHandle)ent->mesh_handle(), 
                           (iMesh::EntityType)ent->dimension(),
                           iMesh_ALL_TOPOLOGIES,
                           ents );
  IBERRCHK(err,"Error getting ModelEnt mesh entities via iMesh");

  if (ent->dimension() > 0) {
    verts.clear();
    fixed_vals.resize(ents.size()+1);
    err = mesh->getEntArrAdj( &ents[0], ents.size(), iBase_VERTEX, verts, &fixed_vals[0] );
    IBERRCHK(err,"Error mesh vertices from mesh entities via iMesh");
  }
  else {
    verts.swap(ents);
  }
  
  fixed_vals.clear();
  fixed_vals.resize( verts.size(), value );
  err = mesh->setIntArrData( &verts[0], verts.size(), fixed, &fixed_vals[0] );
  IBERRCHK(err,"Clearing fixed tag on mesh vertices");
}

void MesquiteOpt::set_fixed_tag_on_skin( ModelEnt* ent, int value )
{
  if (ent->dimension() > 0) {
      // if geometric entity, skin is defined by bounding geometry
    if (ent->geom_handle()) {
      MEntVector children;
      ent->get_adjacencies( ent->dimension() - 1, children );
      for (MEntVector::iterator i = children.begin(); i != children.end(); ++i)
        set_fixed_tag( *i, value );
    }
      // otherwise find the actual skin
    else {
      moab::Interface* mb = mk_core()->moab_instance();
      moab::EntityHandle set = static_cast<moab::EntityHandle>( ent->mesh_handle() );
      moab::ErrorCode rval;
      moab::Range elems;
      rval = mb->get_entities_by_dimension( set, ent->dimension(), elems );
      MBERRCHK(rval, mk_core()->moab_instance());
      
      moab::Range verts;
      moab::Skinner tool(mb);
      moab::ErrorCode err = tool.find_skin( 0, elems, 0, verts );
      MBERRCHK(err, mk_core()->moab_instance());
     
      moab::Tag t = reinterpret_cast<moab::Tag>( get_fixed_tag() );
      err = mb->tag_clear_data( t, verts, &value );
      MBERRCHK(err, mk_core()->moab_instance());
    }
  }
}

bool MesquiteOpt::all_model_ents_have_geom() const
{
  MEntSelection::const_iterator i;
  for (i = mentSelection.begin(); i != mentSelection.end(); ++i) {
    ModelEnt* ent = i->first;
    if (!ent->geom_handle())
      return false;
  }
  return true;
}

void MesquiteOpt::get_adjacent_entity_set( MEntSet& from_this,
                                           MEntVector& ents,
                                           iBase_EntitySetHandle& result,
                                           iBase_EntityType& dimension,
                                           bool& created_result_set )
{
  if (from_this.empty())
    throw Error(MK_BAD_INPUT,"Input set is empty");
  
  iMesh* imesh = mk_core()->imesh_instance();
  bool first = true;
  iBase_ErrorType err;
  created_result_set = false;
  MEntVector front, adj, adj2;
  front.push_back(*from_this.begin());
  from_this.erase(from_this.begin());
  while (!front.empty()) {
    ModelEnt* ent = front.back();
    front.pop_back();
    ents.push_back(ent);
    
    try {
      if (first) {
        dimension = (iBase_EntityType)ent->dimension();
        result = (iBase_EntitySetHandle)ent->mesh_handle();
        first = false;
      }
      else {
        if (ent->dimension() != dimension)
          dimension = iBase_ALL_TYPES;

        if (!created_result_set) {
          err = imesh->unite( result, (iBase_EntitySetHandle)ent->mesh_handle(), result );
          IBERRCHK(err,"WTF?");
          created_result_set = true;
        }
        else {
          iBase_EntitySetHandle tmp;
          err = imesh->unite( result, (iBase_EntitySetHandle)ent->mesh_handle(), tmp );
          IBERRCHK(err,"WTF?");
          imesh->destroyEntSet( result );
          result = tmp;
        }
      }
    }
    catch (...) {
      if (created_result_set)
        imesh->destroyEntSet( result );
      throw;
    }
    
    if (ent->dimension() == 3) {
      adj.clear();
      ent->get_adjacencies( 2, adj );
      for (MEntVector::iterator i = adj.begin(); i != adj.end(); ++i) {
        MEntSet::iterator j = from_this.find(*i);
        if (j != from_this.end()) {
          front.push_back( *j );
          from_this.erase( j );
        }
        adj2.clear();
        (*i)->get_adjacencies( 3, adj2 );
        for (MEntVector::iterator k = adj2.begin(); k != adj2.end(); ++k) {
          j = from_this.find(*k);
          if (j != from_this.end()) {
            front.push_back( *j );
            from_this.erase( j );
          }
        }
      }
    }
    
    adj.clear();
    ent->get_adjacencies( 1, adj );
    for (MEntVector::iterator i = adj.begin(); i != adj.end(); ++i) {
      for (int dim = 2; dim <= 3; ++dim) {
        adj2.clear();
        (*i)->get_adjacencies( dim, adj2 );
        for (MEntVector::iterator k = adj2.begin(); k != adj2.end(); ++k) {
          MEntSet::iterator j = from_this.find(*k);
          if (j != from_this.end()) {
            front.push_back( *j );
            from_this.erase( j );
          }
        }
      }
    }
  }
}


void MesquiteOpt::execute_this()
{
  MEntSelection::iterator i;
  MsqError err;
  iMesh* imesh = mk_core()->imesh_instance();
  iGeom* igeom = mk_core()->igeom_instance();

  IQInterface* smoother;
  DefaultAlgorithm defaultAlgo;
  if (msqAlgo == 0) {
    smoother = &defaultAlgo;
    if (!verboseOutput)
      defaultAlgo.quality_assessor().disable_printing_results();
  }
  else
    smoother = msqAlgo;
    
  create_byte_tag();

  MsqIMesh msqmesh( imesh->instance(), iBase_ALL_TYPES, err, &fixedTag );
  MSQERRCHK(err);

  if (fixedBoundary) {
    for (int dim = 2; dim <= 3; ++dim) {
      for (i = mentSelection.begin(); i != mentSelection.end(); ++i) {
        ModelEnt* ent = i->first;
        if (ent->dimension() != dim)
          continue;
        
        set_fixed_tag( ent, 0 );
        set_fixed_tag_on_skin( ent, 1 );
        MsqIMesh msqmesh( imesh->instance(), 
                          (iBase_EntitySetHandle)ent->mesh_handle(),
                          (iBase_EntityType)dim,
                          err, &fixedTag );
        MSQERRCHK(err);
        
        if (ent->dimension() == 2 && ent->geom_handle()) {
#ifndef MSQIGEOM
          throw Error(MK_BAD_GEOMETRIC_EVALUATION,
              "Mesquite not configured with iGeom support.  "
              "Cannot optimize surface meshes.");
#else
          MsqIGeom msqgeom( igeom->instance(), ent->geom_handle() );
          MeshDomainAssoc mesh_and_domain(&msqmesh, &msqgeom);
          smoother->run_instructions( &mesh_and_domain, err );
#endif
        }
        else {
          MeshDomainAssoc mesh_and_domain(&msqmesh,0);
          smoother->run_instructions( &msqmesh, err );
        }
        MSQERRCHK(err);
      }
    }
  }
  else {
    // shouuld be checked by FreeSmoothDomain
    //if (!all_model_ents_have_geom())
    //  throw Error(MK_BAD_GEOMETRIC_EVALUATION,
    //        "Cannot do free smooth of entity that doesn't have bounding geometry.");

    MEntSet remaining;
    for (i = mentSelection.begin(); i != mentSelection.end(); ++i) {
      ModelEnt* ent = i->first;
      MEntSet::iterator hint = remaining.begin();
      if (ent->dimension() == 2 || ent->dimension() == 3) 
        hint = remaining.insert( hint, ent );
    }
    
    MEntVector ents;
    while (!remaining.empty()) {
      iBase_EntitySetHandle set;
      bool free_set = false;
      iBase_EntityType dim;
      ents.clear();
      get_adjacent_entity_set( remaining, ents, set, dim, free_set );
      
      try {
        FreeSmoothDomain msqgeom( mk_core(), ents );
        MsqIMesh msqmesh( imesh->instance(), set, dim, err, &fixedTag );
        MSQERRCHK(err);
        MeshDomainAssoc mesh_and_domain(&msqmesh, &msqgeom);
        smoother->run_instructions( &mesh_and_domain, err );
        MSQERRCHK(err);
      }
      catch (...) {
        if (free_set)
          imesh->destroyEntSet( set );
        throw;
      }
    }
  }
}

} // namespace MeshKit
  
