#include "meshkit/ParRecvSurfMesh.hpp"
#include "meshkit/ParSendPostSurfMesh.hpp"
#include "moab/Types.hpp"
#include "RefEntity.hpp"
#include "TDParallel.hpp"

namespace MeshKit 
{
// static registration of this mesh scheme
  moab::EntityType ParRecvSurfMesh_tps[] = {moab::MBVERTEX, moab::MBEDGE, moab::MBTRI, moab::MBTET, moab::MBMAXTYPE};
const moab::EntityType* ParRecvSurfMesh::output_types()
{ return ParRecvSurfMesh_tps; }
  
ParRecvSurfMesh::ParRecvSurfMesh(MKCore *mkcore, const MEntVector &me_vec)
  : MeshScheme(mkcore, me_vec)
{
  // get information related to MOAB parallel communication
  m_mpcomm = moab::ParallelComm::get_pcomm(mk_core()->moab_instance(), 0);
  if (NULL == m_mpcomm) {
    throw Error(iBase_FAILURE, "Parallel communication should be already created.");
  }
  m_rank = m_mpcomm->proc_config().proc_rank();
}

ParRecvSurfMesh::~ParRecvSurfMesh()
{
  std::vector<Range*>::iterator vit;
  for (vit = m_shared_entities.begin(); vit != m_shared_entities.end(); vit++)
    delete (*vit);
  m_shared_entities.clear();
}

void ParRecvSurfMesh::setup_this()
{
}

void ParRecvSurfMesh::execute_this()
{
  iMesh::Error err;
  std::set<unsigned int> recv_procs;
  MEntVector recv_ents;

  // get send procs and mesh
  MEntSelection::iterator mit = mentSelection.begin();
  MEntSelection::iterator emit = mentSelection.end();
  for (; mit != emit; mit++) {
    ModelEnt *me = mit->first;
    RefEntity* entity = reinterpret_cast<RefEntity*> (me->geom_handle());
    iBase_EntitySetHandle entityset = reinterpret_cast<iBase_EntitySetHandle> (me->mesh_handle());
    TDParallel *td_par = (TDParallel *) entity->get_TD(&TDParallel::is_parallel);
    if (td_par == NULL) ECERRCHK(MK_FAILURE, "Recv surface should have partitioned information.");

    unsigned int charge_p = td_par->get_charge_proc();
    DLIList<int>* proc_list = td_par->get_shared_proc_list();
    int n_proc_list = proc_list->size();

    if (m_rank == charge_p) { // send proc
      // get all child meshes
      std::vector<iBase_EntityHandle> entities;
      err = mk_core()->imesh_instance()->getEntities(entityset, iBase_ALL_TYPES,
                                                     iMesh_ALL_TOPOLOGIES, entities);
      IBERRCHK(err, "Couldn't get entities of surface entityset.");
      
      // put processors sharing this interface vertex to buffer
      proc_list->reset();
      for (int i = 0; i < n_proc_list; i++) { // to procs
        unsigned int to_proc = proc_list->get_and_step();
        if (to_proc != m_rank) {
          int ind = get_shared_list(to_proc);
          int n_entity = entities.size(); 
          for (int j = 0; j < n_entity; j++) {
            m_shared_entities[ind]->insert(reinterpret_cast<moab::EntityHandle> (entities[j]));
          }
          m_shared_entities[ind]->insert(reinterpret_cast<moab::EntityHandle> (entityset));
        }
      }
    }
    else { // recv procs
      proc_list->reset();
      proc_list->step();
      for (int i = 1; i < n_proc_list; i++) {
        unsigned int proc = proc_list->get_and_step();
        if (proc == m_rank) { // recv proc
          get_shared_list(charge_p);
          recv_procs.insert(charge_p);
          recv_ents.push_back(mit->first);
          break;
        }
      }
    }
  }

  // get send sets
  std::vector<Range*> send_sets;
  int n_proc = m_shared_procs.size();
  for (int i = 0; i < n_proc; i++) {
    Range set_range = m_shared_entities[i]->subset_by_type(MBENTITYSET);
    *m_shared_entities[i] = subtract(*m_shared_entities[i], set_range);
    Range* tmp_range = new Range(set_range);
    send_sets.push_back(tmp_range);
  }

  // recv mesh
  int incoming1 = 0, incoming2 = 0;
  GraphNode* pspsm_node = get_graph()->find_node("ParSendPostSurfMesh");
  if (pspsm_node != NULL) {
    static_cast<ParSendPostSurfMesh*> (pspsm_node)->get_incoming(incoming1, incoming2);
  }
  else ECERRCHK(MK_FAILURE, "ParSendPostSurfMesh should be inserted before.");

  moab::ErrorCode rval = m_mpcomm->recv_entities(recv_procs, incoming1, incoming2, true);
  MBERRCHK(rval, mk_core()->moab_instance());

  // post irecv for surface meshset
  rval = m_mpcomm->post_irecv(m_shared_procs, recv_procs);
  MBERRCHK(rval, mk_core()->moab_instance());

  // send meshset
  incoming1 = recv_procs.size();
  incoming2 = 0;
  rval = m_mpcomm->send_entities(m_shared_procs,
                                 send_sets, incoming1,
                                 incoming2, true);
  MBERRCHK(rval, mk_core()->moab_instance());

  // recv mesh 
  rval = m_mpcomm->recv_entities(recv_procs, incoming1, incoming2, true);
  MBERRCHK(rval, mk_core()->moab_instance());

  // set the model ent mesh with received mesh
  int n_recv = recv_ents.size();
  for (int i = 0; i < n_recv; i++) {
    recv_ents[i]->set_meshed_state(COMPLETE_MESH);
  }
}

int ParRecvSurfMesh::get_shared_list(const int proc)
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
} // namespace MeshKit
