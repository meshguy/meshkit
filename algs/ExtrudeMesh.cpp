#include "ExtrudeMesh.hpp"

#include <cassert>
#include <iostream>

ExtrudeMesh::ExtrudeMesh(iMesh_Instance mesh) : impl_(mesh), copy_(mesh)
{}

ExtrudeMesh::~ExtrudeMesh()
{}

// This is not very C++-y! But it's simple, so I'm ok with it!
double * cross(double *res,double *a,double *b)
{
    res[0] = a[1]*b[2] - a[2]*b[1];
    res[1] = a[2]*b[0] - a[0]*b[2];
    res[2] = a[0]*b[1] - a[1]*b[0];
    return res;
}

double dot(double *a,double *b)
{
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

int ExtrudeMesh::translate(iBase_EntityHandle *src,int size,double *dv,
                           int steps)
{
    int err;
    iBase_EntitySetHandle set;
    iMesh_createEntSet(impl_,false,&set,&err);
    assert(!err);
    
    iMesh_addEntArrToSet(impl_,src,size,set,&err);
    assert(!err);

    int ret = translate(set,dv,steps);

    iMesh_destroyEntSet(impl_,set,&err);
    assert(!err);

    return ret;
}

int ExtrudeMesh::translate(iBase_EntityHandle *src,iBase_EntityHandle *dest,
                           int size,int steps)
{
    return 0;
}

int ExtrudeMesh::translate(iBase_EntitySetHandle src,double *dv,int steps)
{
    int err;

    iBase_EntityHandle *ents=0; int ent_alloc=0, ent_size=0;
    iBase_EntityHandle *adj=0;  int adj_alloc=0, adj_size=0;
    int *indices=0;             int ind_alloc=0, ind_size=0;
    int *offsets=0;             int off_alloc=0, off_size=0;

    iMesh_getAdjEntIndices(impl_,src,iBase_FACE,iMesh_ALL_TOPOLOGIES,
                           iBase_VERTEX,
                           &ents,    &ent_alloc, &ent_size,
                           &adj,     &adj_alloc, &adj_size,
                           &indices, &ind_alloc, &ind_size,
                           &offsets, &off_alloc, &off_size,
                           &err);
    assert(!err);

    int *normals = get_normals(adj,indices,offsets,ent_size,dv);

    int total_count;
    iMesh_getNumOfType(impl_,0,iBase_VERTEX,&total_count,&err);
    assert(!err);

    // TODO: implement my own copy since CopyMesh is overkill
    for(int i=1; i<=steps; i++)
    {
        double v[] = {dv[0]*i, dv[1]*i, dv[2]*i};
        copy_.copy_move_entities(adj,adj_size,v);
    }

    iBase_EntityHandle *all=0;
    int all_alloc=0,all_size=0;
    iMesh_getEntities(impl_,0,iBase_VERTEX,iMesh_POINT,
                      &all,&all_alloc,&all_size,&err);
    assert(err==0);

    iBase_EntityHandle *kids = all+total_count;

    // Make the first row of volumes
    connect_the_dots(normals,indices,offsets,adj,
                     normals,indices,offsets,kids,ent_size);

    // Make the subsequent rows of volumes
    for(int i=1; i<steps; i++)
    {
        connect_the_dots(normals,indices,offsets,kids,
                         normals,indices,offsets,kids+adj_size,ent_size);
        kids += adj_size;
    }

    free(all);
    free(normals);
    free(ents);
    free(adj);
    free(indices);
    free(offsets);
}

int ExtrudeMesh::translate(iBase_EntitySetHandle src,iBase_EntitySetHandle dest,
                           int steps)
{
    int err;

    iBase_EntityHandle *ents=0; int ent_alloc=0, ent_size=0;
    iBase_EntityHandle *adj=0;  int adj_alloc=0, adj_size=0;
    int *indices=0;             int ind_alloc=0, ind_size=0;
    int *offsets=0;             int off_alloc=0, off_size=0;

    iMesh_getAdjEntIndices(impl_,src,iBase_FACE,iMesh_ALL_TOPOLOGIES,
                           iBase_VERTEX,
                           &ents,    &ent_alloc, &ent_size,
                           &adj,     &adj_alloc, &adj_size,
                           &indices, &ind_alloc, &ind_size,
                           &offsets, &off_alloc, &off_size,
                           &err);
    assert(!err);

    iBase_EntityHandle *ents2=0; int ent2_alloc=0, ent2_size=0;
    iBase_EntityHandle *adj2=0;  int adj2_alloc=0, adj2_size=0;
    int *indices2=0;             int ind2_alloc=0, ind2_size=0;
    int *offsets2=0;             int off2_alloc=0, off2_size=0;

    iMesh_getAdjEntIndices(impl_,dest,iBase_FACE,iMesh_ALL_TOPOLOGIES,
                           iBase_VERTEX,
                           &ents2,    &ent2_alloc, &ent2_size,
                           &adj2,     &adj2_alloc, &adj2_size,
                           &indices2, &ind2_alloc, &ind2_size,
                           &offsets2, &off2_alloc, &off2_size,
                           &err);
    assert(!err);

    double dv[] = {0,0,0};
    int *normals  = get_normals(adj, indices, offsets, ent_size, dv);
    int *normals2 = get_normals(adj2,indices2,offsets2,ent2_size,dv);

    if(steps == 1)
    {
        connect_the_dots(normals, indices, offsets, adj,
                         normals2,indices2,offsets2,adj2,ent_size);
    }
    else
    {
        int total_count;
        iMesh_getNumOfType(impl_,0,iBase_VERTEX,&total_count,&err);
        assert(!err);

        // TODO: implement my own copy since CopyMesh is overkill
        for(int i=1; i<steps; i++)
        {
            double v[] = {dv[0]*i, dv[1]*i, dv[2]*i};
            copy_.copy_move_entities(adj,adj_size,v);
        }

        iBase_EntityHandle *all=0;
        int all_alloc=0,all_size=0;
        iMesh_getEntities(impl_,0,iBase_VERTEX,iMesh_POINT,
                          &all,&all_alloc,&all_size,&err);
        assert(err==0);

        iBase_EntityHandle *kids = all+total_count;

        // Make the first row of volumes
        connect_the_dots(normals,indices,offsets,adj,
                         normals,indices,offsets,kids,ent_size);

        // Make the inner rows of volumes
        for(int i=1; i<steps-1; i++)
        {
            connect_the_dots(normals,indices,offsets,kids,
                             normals,indices,offsets,kids+adj_size,ent_size);
            kids += adj_size;
        }

        // Make the final row of volumes
        connect_the_dots(normals, indices, offsets, kids,
                         normals2,indices2,offsets2,adj2,ent_size);

        free(all);
    }

    free(normals); free(normals2);
    free(ents);    free(ents2);
    free(adj);     free(adj2);
    free(indices); free(indices2);
    free(offsets); free(offsets2);
}

// calculate the normals for each face (1 = towards v, -1 = away from v)
// TODO: this can fail with non-convex polygons
int * ExtrudeMesh::get_normals(iBase_EntityHandle *verts,int *indices,
                               int *offsets,int size,double *dv)
{
    int err;
    int *normals = (int*)malloc(size*sizeof(int));
    for(int i=0; i<size; i++)
    {
        double res[3],a[3],b[3];
        double *coords=0;
        int coord_alloc=0,coord_size=0;

        iBase_EntityHandle curr_verts[3];
        for(int j=0; j<3; j++)
            curr_verts[j] = verts[indices[ offsets[i]+j ]];

        iMesh_getVtxArrCoords(impl_,curr_verts,3,iBase_INTERLEAVED,
                              &coords,&coord_alloc,&coord_size,&err);
        for(int j=0; j<3; j++)
        {
            a[j] = coords[1*3 + j] - coords[0*3 + j];
            b[j] = coords[2*3 + j] - coords[1*3 + j];
        }

        normals[i] = (dot( cross(res,a,b),dv ) > 0) ? 1 : -1;
    }

    return normals;
}

void ExtrudeMesh::connect_the_dots(int *pre_norms, int *pre_inds,
                                   int *pre_offs, iBase_EntityHandle *pre,
                                   int *post_norms,int *post_inds,
                                   int *post_offs,iBase_EntityHandle *post,
                                   int size)
{
    int err;

    for(int i=0; i<size; i++)
    {
        int count = pre_offs[i+1] - pre_offs[i];

        // If the normal is facing in the wrong direction (away from the
        // translation) we add the vertices in reverse order. Otherwise, we
        // go in the usual order.
        int pre_first  = (pre_norms [i] == 1) ? pre_offs [i] : pre_offs [i+1];
        int post_first = (post_norms[i] == 1) ? post_offs[i] : post_offs[i+1];
        int pre_d  = pre_norms [i];
        int post_d = post_norms[i];

        int status;
        iBase_EntityHandle out;

        iBase_EntityHandle *nodes = new iBase_EntityHandle[count*2];
        for(int j=0; j<count; j++)
        {
            nodes[j]       = pre [ pre_inds [pre_first  + pre_d *j] ];
            nodes[j+count] = post[ post_inds[post_first + post_d*j] ];
        }

        if(count == 4)      // quad
            iMesh_createEnt(impl_,iMesh_HEXAHEDRON,nodes,8,&out,&status,&err);
        else if(count == 3) // tri
            iMesh_createEnt(impl_,iMesh_PRISM,nodes,6,&out,&status,&err);
        else
            assert(false);

        assert(status==0 && err==0);
        delete[] nodes;
    }
}

int main()
{
    int err;
    iMesh_Instance mesh;
    iMesh_newMesh("",&mesh,&err,0);
    iBase_EntitySetHandle root_set;
    iMesh_getRootSet(mesh,&root_set,&err);

    ExtrudeMesh *ext = new ExtrudeMesh(mesh);

    double verts[] = {
        0, 0, 0,
        1, 0, 0,
        1, 1, 0,
        0, 1, 0,
        2, 1, 0,
        9, 9, 9,
    };
    iBase_EntityHandle *ents=0;
    int size=0,alloc=0;

    iMesh_createVtxArr(mesh,6,iBase_BLOCKED,verts,3*6,&ents,&size,&alloc,&err);
    assert(err == 0);

    iBase_EntityHandle quad;
    int status;
    iMesh_createEnt(mesh,iMesh_QUADRILATERAL,ents,4,&quad,&status,&err);
    assert(err == 0);

    iBase_EntityHandle tri;
    iBase_EntityHandle tri_ents[] = { ents[1], ents[2], ents[4] };
    iMesh_createEnt(mesh,iMesh_TRIANGLE,tri_ents,3,&tri,&status,&err);
    assert(err == 0);

    iBase_EntityHandle faces[] = {quad, tri};
    double v[] = {0,0,1};
    int steps = 5;
    ext->translate(faces,2,v,steps);

    int count;
    iMesh_getNumOfType(mesh,0,iBase_VERTEX,&count,&err);
    assert(err == 0 && count == 5*(steps+1)+1);

    iMesh_getNumOfType(mesh,0,iBase_FACE,&count,&err);
    assert(err == 0 && count == 2);

    iMesh_getNumOfType(mesh,0,iBase_REGION,&count,&err);
    assert(err == 0 && count == steps*2);

    const char *file = "fun.h5m";
    iMesh_save(mesh,root_set,file,"",&err,strlen(file),0);
    assert(err == 0);

    delete ext;
    iMesh_dtor(mesh,&err);
}
