#ifndef EXTRUDEMESH_HPP
#define EXTRUDEMESH_HPP

#include "iMesh_extensions.h"

class ExtrudeMesh
{
public:
    explicit ExtrudeMesh(iMesh_Instance impl);
    virtual ~ExtrudeMesh();

    iMesh_Instance impl() const
    {
        return impl_;
    }

    // TODO: add some copy/extend tag stuff here because metadata is important
    // I guess!

    int translate(iBase_EntityHandle *src,int size,double *dv,int steps);
    int translate(iBase_EntityHandle *src,iBase_EntityHandle *dest,int size,
                  int steps);

    int translate(iBase_EntitySetHandle src,double *dv,int steps);
    int translate(iBase_EntitySetHandle src,iBase_EntitySetHandle dest,
                  int steps);
private:
    int * get_normals(iBase_EntityHandle *verts,int *indices,int *offsets,
                      int size,double *dv);
    void copy_rows(iBase_EntityHandle *src,int src_size,double *dv,int steps,
                   iBase_EntityHandle *dest,int dest_alloc);
    void connect_the_dots(int *pre_normals, int *pre_indices, int *pre_offsets,
                          iBase_EntityHandle *pre,
                          int *post_normals,int *post_indices,int *post_offsets,
                          iBase_EntityHandle *post,
                          int size);

    // technically not necessary since CopyMesh holds an iMesh instance
    iMesh_Instance impl_;
};

#endif
