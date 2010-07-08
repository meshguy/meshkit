#ifndef EXTRUDEMESH_HPP
#define EXTRUDEMESH_HPP

#include "iMesh_extensions.h"
#include "CopyVerts.hpp"
#include "CopyMesh.hpp"
#include "LocalTag.hpp"

class ExtrudeMesh
{
public:
  explicit ExtrudeMesh(iMesh_Instance impl);
  virtual ~ExtrudeMesh();

  iMesh_Instance impl() const
  {
    return impl_;
  }

  iBase_TagHandle copy_tag();
  iBase_TagHandle extrude_tag();

  int add_copy_tag   (const std::string &tag_name, const char *tag_val = NULL);
  int add_copy_tag   (iBase_TagHandle tag_handle,  const char *tag_val = NULL);
  int add_expand_tag (const std::string &tag_name, const char *tag_val = NULL);
  int add_expand_tag (iBase_TagHandle tag_handle,  const char *tag_val = NULL);
//  int add_unique_tag (const std::string &tag_name);
//  int add_unique_tag (iBase_TagHandle tag_handle);
  int add_extrude_tag(const std::string &tag_name, const char *tag_val = NULL);
  int add_extrude_tag(iBase_TagHandle tag_handle,  const char *tag_val = NULL);

  std::set<iBase_EntitySetHandle> & copy_sets();
  std::set<iBase_EntitySetHandle> & expand_sets();
//  std::set<iBase_EntitySetHandle> & unique_sets();
  std::set<iBase_EntitySetHandle> & extrude_sets();

  int update_sets();
  int reset_sets();

  int translate(iBase_EntityHandle *src, int size, int steps, const double *dx,
                bool copy_faces = false);
  int translate(iBase_EntitySetHandle src, int steps, const double *dx,
                bool copy_faces = false);

  int rotate(iBase_EntityHandle *src, int size, int steps, const double *origin,
             const double *z, double angle, bool copy_faces = false);
  int rotate(iBase_EntitySetHandle src, int steps, const double *origin,
             const double *z, double angle, bool copy_faces = false);

  int extrude(iBase_EntityHandle *src, int size, int steps,
              const CopyVerts &trans, bool copy_faces = false);
  int extrude(iBase_EntitySetHandle src, int steps, const CopyVerts &trans,
              bool copy_faces = false);

  // Forwards from CopyMesh
/*  int add_copy_expand_list(iBase_EntitySetHandle *ce_sets, int num_ce_sets,
                           int copy_or_expand);
                           int reset_ce_lists();*/

private:
  int do_extrusion(iBase_EntitySetHandle src, iBase_EntitySetHandle dest,
           bool use_dest, int inner_rows, const CopyVerts &trans);

  int * get_normals(iBase_EntityHandle *verts, int *indices, int *offsets,
                    int size, double *dv);

  void connect_higher_dots(iBase_EntityHandle *src, int size,
                           iBase_TagHandle local_extrude_tag, int *normals,
                           int *indices, int *offsets, iBase_EntityHandle *pre,
                           iBase_EntityHandle *post);

  iMesh_Instance impl_;

  struct tag_data
  {
    tag_data(iBase_TagHandle tag, char *value)
      : tag(tag), value(value)
    {}

    iBase_TagHandle tag;
    char *value;
  };

  bool updated_set_;
  LocalTag copy_tag_;
  std::vector<tag_data> copy_tags_;
  std::set<iBase_EntitySetHandle> copy_sets_;
  std::vector<tag_data> expand_tags_;
  std::set<iBase_EntitySetHandle> expand_sets_;


  LocalTag extrude_tag_;
  std::vector<tag_data> extrude_tags_;
  std::set<iBase_EntitySetHandle> extrude_sets_;
};

inline iBase_TagHandle
ExtrudeMesh::extrude_tag()
{
  return extrude_tag_;
}

inline iBase_TagHandle
ExtrudeMesh::copy_tag()
{
  return copy_tag_;
}

inline std::set<iBase_EntitySetHandle> &
ExtrudeMesh::copy_sets()
{
  return copy_sets_;
}

inline std::set<iBase_EntitySetHandle> &
ExtrudeMesh::expand_sets()
{
  return expand_sets_;
}

inline std::set<iBase_EntitySetHandle> &
ExtrudeMesh::extrude_sets()
{
  return extrude_sets_;
}

/*inline int
ExtrudeMesh::add_copy_expand_list(iBase_EntitySetHandle *ce_sets, 
                                  int num_ce_sets, int copy_or_expand)
{
  return copy_.add_copy_expand_list(ce_sets, num_ce_sets, copy_or_expand);
}

inline int
ExtrudeMesh::add_unique_tag(const std::string &tag_name)
{
  return copy_.add_unique_tag(tag_name);
}

inline int
ExtrudeMesh::add_unique_tag(iBase_TagHandle tag_handle)
{
  return copy_.add_unique_tag(tag_handle);
}

inline std::set<iBase_EntitySetHandle> &
ExtrudeMesh::unique_sets()
{
  return copy_.unique_sets();
  }*/

#endif
