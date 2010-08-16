#ifndef EXTRUDEMESH_HPP
#define EXTRUDEMESH_HPP

#include "iMesh_extensions.h"
#include "Transform.hpp"
#include "LocalTag.hpp"
#include "CESets.hpp"

class ExtrudeMesh
{
public:
  explicit ExtrudeMesh(iMesh_Instance impl);
  virtual ~ExtrudeMesh();

  iMesh_Instance impl() const
  {
    return impl_;
  }

  iBase_TagHandle extrude_tag();
  iBase_TagHandle copy_tag();

  CESets & extrude_sets();
  CESets & copy_sets();
  CESets & expand_sets();

  /*void add_unique_tag(const std::string &tag_name);
  void add_unique_tag(iBase_TagHandle tag_handle);

  std::set<iBase_EntitySetHandle> & unique_sets();*/

  void update_sets();
  void reset_sets();

  void translate(iBase_EntityHandle *src, iBase_EntityHandle *dest, int size,
                 int steps);
  void translate(iBase_EntitySetHandle src, iBase_EntitySetHandle dest,
                 int steps);

  void extrude(iBase_EntityHandle *src, int size,
               const extrude::Transform &trans, bool copy_faces = false);
  void extrude(iBase_EntitySetHandle src,
               const extrude::Transform &trans, bool copy_faces = false);
  void extrude(iBase_EntityHandle *src, iBase_EntityHandle *dest,
               int steps, const extrude::Transform &trans);
  void extrude(iBase_EntitySetHandle src, iBase_EntitySetHandle dest,
               const extrude::Transform &trans);
private:
  std::vector<int> get_normals(iBase_EntityHandle *verts, int *indices,
                               int *offsets, int size, double *dv);

  void connect_up_dots(iBase_EntityHandle *src, int size,
                       iBase_TagHandle local_tag, int *norms, int *inds,
                       int *offs, iBase_EntityHandle *pre,
                       iBase_EntityHandle *post);
  void connect_up_dots(
    iBase_EntityHandle *src, int size, iBase_TagHandle local_tag,
    int *pre_norms,  int *pre_inds,  int *pre_offs,  iBase_EntityHandle *pre,
    int *post_norms, int *post_inds, int *post_offs, iBase_EntityHandle *post);

  iMesh_Instance impl_;

  LocalTag copy_tag_;
  LocalTag extrude_tag_;

  CESets copy_sets_;
  CESets expand_sets_;
  CESets extrude_sets_;
};

inline iBase_TagHandle
ExtrudeMesh::copy_tag()
{
  return copy_tag_;
}

inline iBase_TagHandle
ExtrudeMesh::extrude_tag()
{
  return extrude_tag_;
}

inline CESets &
ExtrudeMesh::copy_sets()
{
  return copy_sets_;
}

inline CESets &
ExtrudeMesh::expand_sets()
{
  return expand_sets_;
}

inline CESets &
ExtrudeMesh::extrude_sets()
{
  return extrude_sets_;
}

/*inline void
ExtrudeMesh::add_unique_tag(const std::string &tag_name)
{
  copy_.add_unique_tag(tag_name);
}

inline void
ExtrudeMesh::add_unique_tag(iBase_TagHandle tag_handle)
{
  copy_.add_unique_tag(tag_handle);
}

inline std::set<iBase_EntitySetHandle> &
ExtrudeMesh::unique_sets()
{
  return copy_.unique_sets();
}
*/

inline void
ExtrudeMesh::connect_up_dots(iBase_EntityHandle *src, int size,
                             iBase_TagHandle local_tag, int *norms, int *inds,
                             int *offs, iBase_EntityHandle *pre,
                             iBase_EntityHandle *post)
{
  connect_up_dots(src, size, local_tag,
                  norms, inds, offs, pre,
                  norms, inds, offs, post);
}

#endif
