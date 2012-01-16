#include "meshkit/ParExchangeMesh.hpp"
#include "meshkit/ParPostRecv.hpp"
#include "moab/Range.hpp"
#include "moab/Types.hpp"
#include "RefEntity.hpp"
#include "TDParallel.hpp"

const bool debug_par_exchange_mesh = false;

namespace MeshKit 
{
// static registration of this mesh scheme
  moab::EntityType ParExchangeMesh_tps[] = {moab::MBVERTEX, moab::MBEDGE, moab::MBTRI, moab::MBTET, moab::MBMAXTYPE};
const moab::EntityType* ParExchangeMesh::output_types()
{ return ParExchangeMesh_tps; }
  
ParExchangeMesh::ParExchangeMesh(MKCore *mkcore, const MEntVector &me_vec)
  : MeshScheme(mkcore, me_vec)
{
  // get information related to MOAB parallel communication
  m_mpcomm = moab::ParallelComm::get_pcomm(mk_core()->moab_instance(), 0);
  if (NULL == m_mpcomm) {
    throw Error(iBase_FAILURE, "Parallel communication should be already created.");
  }
  m_rank = m_mpcomm->proc_config().proc_rank();
  m_proc_size = m_mpcomm->proc_config().proc_size();

  // create tag
  iMesh::Error err = mk_core()->imesh_instance()->createTag("PARALLEL_UNIQUE_ID", 1, iBase_INTEGER, m_mPuniqueIDTag);
  if (err != iBase_TAG_ALREADY_EXISTS) {
    IBERRCHK(err, "Trouble create a parallel unique id tag handle.");
  }
}

ParExchangeMesh::~ParExchangeMesh()
{
  //delete m_mpcomm;
  std::vector<Range*>::iterator vit;
  for (vit = m_shared_entities.begin(); vit != m_shared_entities.end(); vit++)
    delete (*vit);
  m_shared_entities.clear();
}

void ParExchangeMesh::setup_this()
{
}

void ParExchangeMesh::execute_this()
{
  iMesh::Error err;
  MEntVector meshed_mes;

  for (MEntSelection::iterator mit = mentSelection.begin(); mit != mentSelection.end(); mit++) {
    ModelEnt *me = mit->first;
    RefEntity* entity = reinterpret_cast<RefEntity*> (me->geom_handle());
    iBase_EntitySetHandle entityset = reinterpret_cast<iBase_EntitySetHandle> (me->mesh_handle());
    TDParallel *td_par = (TDParallel *) entity->get_TD(&TDParallel::is_parallel);
    if (td_par == NULL) ECERRCHK(MK_FAILURE, "Exchange entity should have partitioned information.");
    int charge_p = td_par->get_charge_proc();

    if (m_rank == charge_p) { // send proc
      // get all child meshes
      std::vector<iBase_EntityHandle> entities;
      err = mk_core()->imesh_instance()->getEntities(entityset, iBase_ALL_TYPES,
                                                     iMesh_ALL_TOPOLOGIES, entities);
      IBERRCHK(err, "Couldn't get entities of surface entityset.");
      
      // put processors sharing this interface vertex to buffer
      DLIList<int>* proc_list = td_par->get_shared_proc_list();
      int n_proc_list = proc_list->size();
      proc_list->reset();
      
      for (int i = 0; i < n_proc_list; i++) {
        int proc = proc_list->get_and_step();
        if (proc != m_rank) {
          int ind = get_shared_list(proc);
          int n_entity = entities.size(); 
          for (int j = 0; j < n_entity; j++) {
            m_shared_entities[ind]->insert(reinterpret_cast<moab::EntityHandle> (entities[j]));
          }
          m_shared_entities[ind]->insert(reinterpret_cast<moab::EntityHandle> (entityset));
        }
      }
    }
    else {
      get_shared_list(charge_p); // receive proc: just create shared entity/proc list
      meshed_mes.push_back(me);
    }
  }

  // exchange shared entities
  if (debug_par_exchange_mesh) {
    std::cout << "m_shared_procs_size=" << m_shared_procs.size() << std::endl;
    for (int i = 0; i < m_shared_procs.size(); i++) {
      std::cout << "m_shared_procs[" << i << "]=" << m_shared_procs[i] << std::endl;
    }
    std::cout << "m_shared_entities_size=" << m_shared_entities.size() << std::endl;
    for (int i = 0; i < m_shared_entities.size(); i++) {
      std::cout << "m_shared_entities_range[" << i << "]_size=" << m_shared_entities[i]->size() << std::endl;
    }
  }

  int dim = mentSelection.begin()->first->dimension();
  moab::ErrorCode rval = m_mpcomm->exchange_owned_meshs(m_shared_procs,
                                                        m_shared_entities,
                                                        m_recv_reqs,
                                                        m_recv_remoteh_reqs,
                                                        true, false, false, dim);
  MBERRCHK(rval, mk_core()->moab_instance());
  
  // set the model ent mesh with received mesh
  int n_meshed = meshed_mes.size();
  for (int i = 0; i < n_meshed; i++) {
    meshed_mes[i]->set_meshed_state(COMPLETE_MESH);
  }

  if (debug_par_exchange_mesh) {
    int dim = -1;
    if (mentSelection.size() > 0) dim = mentSelection.begin()->first->dimension();
    std::cout << "Parallel_exchange_execution_time(dim:" << dim << ")=" << t2 - t1 << std::endl;
  }
}

int ParExchangeMesh::get_shared_list(const int proc)
{
  int ind = -1;
  std::vector<unsigned int>::iterator vit = 
    std::find(m_shared_procs.begin(), m_shared_procs.end(), proc);
  if (vit == m_shared_procs.end()) {
    ind = m_shared_procs.size();
    m_shared_procs.push_back(proc);
    m_shared_entities.push_back(new Range);
  }
  else ind = vit - m_shared_procs.begin();

  return ind;
}

void ParExchangeMesh::print_mesh()
{
  // test
  int tmp_procs[MAX_SHARING_PROCS];
  moab::EntityHandle tmp_hs[MAX_SHARING_PROCS];
  unsigned char pstat;
  int num_ps;
  moab::Range entities;
  moab::ErrorCode rval;

  for (moab::EntityType type = MBVERTEX; type != MBMAXTYPE; type++) {
    entities.clear();
    rval = mk_core()->moab_instance()->get_entities_by_type(NULL, type, entities);
    MBERRCHK(rval, mk_core()->moab_instance());

    for (moab::Range::iterator rit = entities.begin(); rit != entities.end(); rit++) {
      rval = m_mpcomm->get_sharing_data(*rit, tmp_procs, tmp_hs, pstat, num_ps);
      MBERRCHK(rval, mk_core()->moab_instance());

      //std::cout << "ParExchangeMesh::entity=" << *rit << ", type=" << type;
      if (type == MBVERTEX) {
        double coord[3];
        rval = mk_core()->moab_instance()->get_coords(&(*rit), 1, coord);
        MBERRCHK(rval, mk_core()->moab_instance());
        //std::cout << ", coords=" << coord[0] << " " << coord[1] << " " << coord[2];
      }
      else if (type != MBENTITYSET) {
        std::vector<moab::EntityHandle> connect;
        rval = mk_core()->moab_instance()->get_connectivity(&(*rit), 1, connect);
        MBERRCHK(rval, mk_core()->moab_instance());
        int n_conn = connect.size();
        //std::cout << ", connect=";
        //for (int j = 0; j < n_conn; j++) {
        //std::cout << connect[j] << " ";
        //}
      }
      /*
      std::cout << ", shared_info=";
      for (int ii = 0; ii < num_ps; ii++) {
        std::cout << tmp_procs[ii] << ":" << tmp_hs[ii] << ", ";
      }
      std::cout << std::endl;
      */
    }
  }    
}
} // namespace MeshKit
