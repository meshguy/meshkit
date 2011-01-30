#ifndef COPYGEOM_HPP
#define COPYGEOM_HPP

#include "Transform.hpp"
#include "LocalTag.hpp"
#include "CESets.hpp"
#include "MKException.hpp"

#include <cassert>
#include <string>
#include <vector>
#include <set>

class CopyGeom
{
public:
  /* \brief Constructor
   *
   * Create a new CopyGeom instance
   * \param impl the iGeom instance handle for the Geom
   */
  CopyGeom(iGeom_Instance impl);

  /* \brief Destructor
   */
  virtual ~CopyGeom();

  /* \brief copy/move igeom entities
   *
   * to location dx
   * \param entities, entity handle size and the location
   */
  void copy(iBase_EntityHandle *entities, const int entities_ehsize,
	    const double *dx,
	    iBase_EntityHandle **new_ents = NULL,
	    int *new_ents_allocated = 0,
	    int *new_ents_size = 0,
	    bool do_merge = false);

  void tag_copied_sets(const char **tag_names, const char **tag_vals,
                       const int num_tags);

  iGeom_Instance igeomImpl;
};
#endif
