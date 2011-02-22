#ifndef MESHKIT_COPY_UTILS_HPP
#define MESHKIT_COPY_UTILS_HPP

#include "meshkit/iMesh.hpp"
#include <vector>

/**\brief Get the entities and unique adjacent vertices of a set
 *
 * Return an array of all entities in a set, sorted by topology, and all
 * unique adjacent vertices. The adjacent vertices can be accessed with an
 * offset array indexing into an index buffer.
 *
 * \param instance the iMesh instance handle
 * \param set the set of entities from which to query
 * \param ents an array of entity handles in the set
 * \param unique_adj an array of unique vertices adjacent to |ents|
 * \param indices index buffer into |unique_adj|
 * \param offsets offset array indicating start and end of |indices| for each
                  entity in |ents|
 */
iMesh::Error iMesh_getStructure(iMesh_Instance instance,
                                iBase_EntitySetHandle set,
                                std::vector<iBase_EntityHandle> &ents,
                                std::vector<iBase_EntityHandle> &unique_adj,
                                std::vector<int> &indices,
                                std::vector<int> &offsets);

#endif
