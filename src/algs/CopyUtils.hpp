#ifndef COPYUTILS_HPP
#define COPYUTILS_HPP

#include <iMesh.h>

/**\brief Get the entities and unique adjacent vertices of a set
 *
 * Return an array of all entities in a set, sorted by topology, and all
 * unique adjacent vertices. The adjacent vertices can be accessed with an
 * offset array indexing into an index buffer.
 *
 * \param instance the iMesh instance handle
 * \param set the set of entities from which to query
 * \param ents pointer to an array of entity handles in the set
 * \param ents_allocated allocated space of |ents|
 * \param ents_size number of elements in |ents|
 * \param unique_adj pointer to an array of unique vertices adjacent to |ents|
 * \param unique_adj_allocated allocated space of |unique_adj|
 * \param unique_adj_size number of elements in |unique_adj|
 * \param indices index buffer into |unique_adj|
 * \param indices_allocated allocated space of |indices|
 * \param indices_size number of elements in |indices|
 * \param offsets offset array indicating start and end of |indices| for each
                  entity in |ents|
 * \param offsets_allocated allocated space of |offsets|
 * \param offsets_size number of elements in |offsets|
 * \param err pointer to an ITAPS error code 
 */
void iMesh_getStructure(iMesh_Instance instance, iBase_EntitySetHandle set,
                        iBase_EntityHandle **ents,
                        int *ents_allocated,
                        int *ents_size,
                        iBase_EntityHandle **unique_adj,
                        int *unique_adj_allocated,
                        int *unique_adj_size,
                        int **indices,
                        int *indices_allocated,
                        int *indices_size,
                        int **offsets,
                        int *offsets_allocated,
                        int *offsets_size,
                        int *err);

/**\brief Create higher-dimension entities on a set of vertices
 *
 * Given an array of vertices and a template to build from, recreate the
 * higher-dimension entities on the array of vertices. This is useful when
 * making copies of a selection.
 *
 * \param instance the iMesh instance handle
 * \param ents the entities for the template
 * \param size the number of entities in |entities|
 * \param local_tag the local copy tag
 * \param indices index buffer into |entities|
 * \param offsets offset array indicating start and end of |indices| for each
                  entity in |entities|
 * \param verts the array of vertices to build upon
 */
void connect_the_dots(iMesh_Instance instance, iBase_EntityHandle *ents,
                      int size, iBase_TagHandle local_tag, int *indices,
                      int *offsets, iBase_EntityHandle *verts);

#endif
