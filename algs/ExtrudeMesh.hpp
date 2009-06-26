#ifndef EXTRUDEMESH_HPP
#define EXTRUDEMESH_HPP

#include "iMesh_extensions.h"
#include "CopyVerts.hpp"
#include "CopyMesh.hpp"

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
    int add_copy_expand_list(iBase_EntitySetHandle *ce_sets, int num_ce_sets,
                             int copy_or_expand);
    int reset_ce_lists();
    int add_copy_tag  (const std::string &tag_name, const char *tag_val = NULL);
    int add_copy_tag  (iBase_TagHandle tag_handle,  const char *tag_val = NULL);
    int add_expand_tag(const std::string &tag_name, const char *tag_val = NULL);
    int add_expand_tag(iBase_TagHandle tag_handle,  const char *tag_val = NULL);
    int add_unique_tag(const std::string &tag_name);
    int add_unique_tag(iBase_TagHandle tag_handle);
    std::set<iBase_EntitySetHandle> & copy_sets();
    std::set<iBase_EntitySetHandle> & expand_sets();
    std::set<iBase_EntitySetHandle> & unique_sets();

    int translate(iBase_EntityHandle *src,int size,int steps,const double *dx,
                  bool copy_faces = false);
    int translate(iBase_EntitySetHandle src,int steps,const double *dx,
                  bool copy_faces = false);
    int translate(iBase_EntityHandle *src,iBase_EntityHandle *dest,int size,
                  int steps);
    int translate(iBase_EntitySetHandle src,iBase_EntitySetHandle dest,
                  int steps);

    int rotate(iBase_EntityHandle *src,int size,int steps,const double *origin,
               const double *z,double angle,bool copy_faces = false);
    int rotate(iBase_EntitySetHandle src,int steps,const double *origin,
               const double *z,double angle,bool copy_faces = false);
    int rotate(iBase_EntityHandle *src,iBase_EntityHandle *dest,int size,
               int steps,const double *origin,const double *z,double angle);
    int rotate(iBase_EntitySetHandle src,iBase_EntitySetHandle dest,int steps,
               const double *origin,const double *z,double angle);


    int extrude(iBase_EntityHandle *src,int size,int steps,
                const CopyVerts &trans,bool copy_faces = false);
    int extrude(iBase_EntitySetHandle src,int steps,const CopyVerts &trans,
                bool copy_faces = false);
    int extrude(iBase_EntityHandle *src,iBase_EntityHandle *dest,int size,
                int steps,const CopyVerts &trans);
    int extrude(iBase_EntitySetHandle src,iBase_EntitySetHandle dest,
                int steps,const CopyVerts &trans);
private:
    int do_extrusion(iBase_EntitySetHandle src,iBase_EntitySetHandle dest,
                     bool use_dest,int inner_rows,const CopyVerts &trans);

    int * get_normals(iBase_EntityHandle *verts,int *indices,int *offsets,
                      int size,double *dv);

    void connect_the_dots(int *pre_normals, int *pre_indices, int *pre_offsets,
                          iBase_EntityHandle *pre,
                          int *post_normals,int *post_indices,int *post_offsets,
                          iBase_EntityHandle *post,
                          int size);

    iMesh_Instance impl_;
    CopyMesh copy_;
};

inline iBase_TagHandle
ExtrudeMesh::copy_tag()
{
    return copy_.copy_tag();
}

inline int
ExtrudeMesh::add_copy_expand_list(iBase_EntitySetHandle *ce_sets, 
                                  int num_ce_sets,int copy_or_expand)
{
    return copy_.add_copy_expand_list(ce_sets,num_ce_sets,copy_or_expand);
}

inline int
ExtrudeMesh::reset_ce_lists()
{
    return copy_.reset_ce_lists();
}

inline int
ExtrudeMesh::add_copy_tag(const std::string &tag_name,const char *tag_val)
{
    return copy_.add_copy_tag(tag_name,tag_val);
}

inline int
ExtrudeMesh::add_copy_tag(iBase_TagHandle tag_handle,const char *tag_val)
{
    return copy_.add_copy_tag(tag_handle,tag_val);
}

inline int
ExtrudeMesh::add_expand_tag(const std::string &tag_name,const char *tag_val)
{
    return copy_.add_expand_tag(tag_name,tag_val);
}

inline int
ExtrudeMesh::add_expand_tag(iBase_TagHandle tag_handle,const char *tag_val)
{
    return copy_.add_expand_tag(tag_handle,tag_val);
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
ExtrudeMesh::copy_sets()
{
    return copy_.copy_sets();
}

inline std::set<iBase_EntitySetHandle> &
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
