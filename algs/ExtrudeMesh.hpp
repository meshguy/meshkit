#ifndef EXTRUDEMESH_HPP
#define EXTRUDEMESH_HPP

#include "iMesh_extensions.h"
#include "CopyVerts.hpp"
#include "CopyMesh.hpp"
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

  CESets & extrude_sets();

  void update_sets();
  void reset_sets();

  void translate(iBase_EntityHandle *src, int size, int steps, const double *dx,
                 bool copy_faces = false);
  void translate(iBase_EntitySetHandle src, int steps, const double *dx,
                 bool copy_faces = false);
  void translate(iBase_EntityHandle *src, iBase_EntityHandle *dest, int size,
                 int steps);
  void translate(iBase_EntitySetHandle src, iBase_EntitySetHandle dest,
                 int steps);

  void rotate(iBase_EntityHandle *src, int size, int steps,
              const double *origin, const double *z, double angle,
              bool copy_faces = false);
  void rotate(iBase_EntitySetHandle src, int steps, const double *origin,
              const double *z, double angle, bool copy_faces = false);
  void rotate(iBase_EntityHandle *src, iBase_EntityHandle *dest, int size,
              int steps, const double *origin, const double *z, double angle);
  void rotate(iBase_EntitySetHandle src, iBase_EntitySetHandle dest, int steps,
              const double *origin, const double *z, double angle);


  void extrude(iBase_EntityHandle *src, int size, int steps,
               const CopyVerts &trans, bool copy_faces = false);
  void extrude(iBase_EntitySetHandle src, int steps, const CopyVerts &trans,
               bool copy_faces = false);
  void extrude(iBase_EntityHandle *src, iBase_EntityHandle *dest, int size,
               int steps, const CopyVerts &trans);
  void extrude(iBase_EntitySetHandle src, iBase_EntitySetHandle dest, int steps,
               const CopyVerts &trans);

  // Forwards from CopyMesh
  iBase_TagHandle copy_tag();

  void add_unique_tag(const std::string &tag_name);
  void add_unique_tag(iBase_TagHandle tag_handle);

  CESets & copy_sets();
  CESets & expand_sets();
  std::set<iBase_EntitySetHandle> & unique_sets();
private:
  void do_extrusion(iBase_EntitySetHandle src, iBase_EntitySetHandle dest,
                    bool use_dest, int inner_rows, const CopyVerts &trans);

  std::vector<int> get_normals(iBase_EntityHandle *verts, int *indices,
                               int *offsets, int size, double *dv);

  void connect_the_dots(
    iBase_EntityHandle *src, int size, iBase_TagHandle local_tag, 
    int *pre_norms,  int *pre_inds,  int *pre_offs,  iBase_EntityHandle *pre,
    int *post_norms, int *post_inds, int *post_offs, iBase_EntityHandle *post);

  iMesh_Instance impl_;
  CopyMesh copy_;

  bool updated_set_;
  LocalTag extrude_tag_;
  CESets extrude_sets_;
};

inline iBase_TagHandle
ExtrudeMesh::extrude_tag()
{
  return extrude_tag_;
}

inline iBase_TagHandle
ExtrudeMesh::copy_tag()
{
  return copy_.copy_tag();
}

inline void
ExtrudeMesh::add_unique_tag(const std::string &tag_name)
{
  copy_.add_unique_tag(tag_name);
}

inline void
ExtrudeMesh::add_unique_tag(iBase_TagHandle tag_handle)
{
  copy_.add_unique_tag(tag_handle);
}

inline CESets &
ExtrudeMesh::extrude_sets()
{
  return extrude_sets_;
}

inline CESets &
ExtrudeMesh::copy_sets()
{
  return copy_.copy_sets();
}

inline CESets &
ExtrudeMesh::expand_sets()
{
  return copy_.expand_sets();
}

inline std::set<iBase_EntitySetHandle> &
ExtrudeMesh::unique_sets()
{
  return copy_.unique_sets();
}

#endif
