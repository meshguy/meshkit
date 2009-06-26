#include "ExtrudeMesh.hpp"

#include <cassert>
#include <iostream>

ExtrudeMesh::ExtrudeMesh(iMesh_Instance mesh) : impl_(mesh),copy_(mesh)
{}

ExtrudeMesh::~ExtrudeMesh()
{}

// This is not very C++-y! But it's simple, so I'm ok with it!
static double * cross(double *res,double *a,double *b)
{
    res[0] = a[1]*b[2] - a[2]*b[1];
    res[1] = a[2]*b[0] - a[0]*b[2];
    res[2] = a[0]*b[1] - a[1]*b[0];
    return res;
}

static double dot(double *a,double *b)
{
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

static double * vtx_diff(double *res,iMesh_Instance mesh,iBase_EntityHandle a,
                         iBase_EntityHandle b)
{
    int err;
    double xa[3],xb[3];

    iMesh_getVtxCoord(mesh,a,xa+0,xa+1,xa+2,&err);
    assert(!err);
    iMesh_getVtxCoord(mesh,b,xb+0,xb+1,xb+2,&err);
    assert(!err);

    for(int i=0; i<3; i++)
        res[i] = xa[i]-xb[i];
    return res;
}

int ExtrudeMesh::translate(iBase_EntityHandle *src,int size,int steps,
                           const double *dx,bool copy_faces)
{
    return extrude(src,size,steps,CopyMoveVerts(impl_,dx,steps),copy_faces);
}

int ExtrudeMesh::translate(iBase_EntitySetHandle src,int steps,const double *dx,
                           bool copy_faces)
{
    return extrude(src,steps,CopyMoveVerts(impl_,dx,steps),copy_faces);
}

int ExtrudeMesh::translate(iBase_EntityHandle *src,iBase_EntityHandle *dest,
                           int size,int steps)
{
    int err;
    iBase_EntitySetHandle src_set,dest_set;
    iMesh_createEntSet(impl_,false,&src_set,&err);
    assert(!err);
    iMesh_createEntSet(impl_,false,&dest_set,&err);
    assert(!err);
    
    iMesh_addEntArrToSet(impl_,src,size,src_set,&err);
    assert(!err);
    iMesh_addEntArrToSet(impl_,dest,size,dest_set,&err);
    assert(!err);

    int ret = translate(src_set,dest_set,steps);

    iMesh_destroyEntSet(impl_,src_set,&err);
    assert(!err);
    iMesh_destroyEntSet(impl_,dest_set,&err);
    assert(!err);

    return ret;
}

int ExtrudeMesh::translate(iBase_EntitySetHandle src,iBase_EntitySetHandle dest,
                           int steps)
{
    int err;

    // Deduce the per-step displacement vector "dx"
    // Note: we assume that src and dest are the same shape, etc.
    double dx[3];
    double coords[2][3];

    iBase_EntitySetHandle ends[] = {src,dest};
    for(int i=0; i<2; i++)
    {
        iMesh_EntityIterator iter;
        iMesh_initEntIter(impl_,ends[i],iBase_FACE,iMesh_ALL_TOPOLOGIES,&iter,
                          &err);
        assert(!err);

        iBase_EntityHandle face;
        int has_data;
        iMesh_getNextEntIter(impl_,iter,&face,&has_data,&err);
        assert(!err);
        assert(has_data);

        iMesh_endEntIter(impl_,iter,&err);
        assert(!err);

        iBase_EntityHandle *verts=0;
        int verts_alloc=0,verts_size=0;
        iMesh_getEntAdj(impl_,face,iBase_VERTEX,&verts,&verts_alloc,
                        &verts_size,&err);
        assert(!err);

        iMesh_getVtxCoord(impl_,verts[0],coords[i]+0,coords[i]+1,coords[i]+2,
                          &err);
        assert(!err);

        free(verts);
    }

    for(int i=0; i<3; i++)
        dx[i] = (coords[1][i]-coords[0][i])/steps;

    return extrude(src,dest,steps,CopyMoveVerts(impl_,dx));
}

int ExtrudeMesh::rotate(iBase_EntityHandle *src,int size,int steps,
                        const double *origin,const double *z,double angle,
                        bool copy_faces)
{
    return extrude(src,size,steps,CopyRotateVerts(impl_,origin,z,angle,steps),
                   copy_faces);
}

int ExtrudeMesh::rotate(iBase_EntitySetHandle src,int steps,
                        const double *origin,const double *z,double angle,
                        bool copy_faces)
{
    return extrude(src,steps,CopyRotateVerts(impl_,origin,z,angle),copy_faces);
}

int ExtrudeMesh::rotate(iBase_EntityHandle *src,iBase_EntityHandle *dest,
                        int size,int steps,const double *origin,
                        const double *z,double angle)
{
    return extrude(src,dest,size,steps,CopyRotateVerts(impl_,origin,z,angle));
}

int ExtrudeMesh::rotate(iBase_EntitySetHandle src,iBase_EntitySetHandle dest,
                        int steps,const double *origin,const double *z,
                        double angle)
{
    return extrude(src,dest,steps,CopyRotateVerts(impl_,origin,z,angle));
}

int ExtrudeMesh::extrude(iBase_EntityHandle *src,int size,int steps,
                         const CopyVerts &trans,bool copy_faces)
{
    int err;
    iBase_EntitySetHandle set;
    iMesh_createEntSet(impl_,false,&set,&err);
    assert(!err);
    
    iMesh_addEntArrToSet(impl_,src,size,set,&err);
    assert(!err);

    int ret = extrude(set,steps,trans,copy_faces);

    iMesh_destroyEntSet(impl_,set,&err);
    assert(!err);

    return ret;
}

int ExtrudeMesh::extrude(iBase_EntitySetHandle src,int steps,
                         const CopyVerts &trans,bool copy_faces)
{
    if(copy_faces)
    {
        int err;

        iBase_EntitySetHandle parent;
        iMesh_createEntSet(impl_,false,&parent,&err);
        assert(!err);
        iMesh_addEntSet(impl_,src,parent,&err);
        assert(!err);
        
        copy_.add_copy_expand_list(&src,   1,CopyMesh::COPY);
        copy_.add_copy_expand_list(&parent,1,CopyMesh::EXPAND);
        copy_.copy_transform_entities(src,trans,0,0,0);

        iMesh_rmvEntSet(impl_,src,parent,&err);
        assert(!err);

        iBase_EntitySetHandle dest;
        iBase_EntitySetHandle *tmp = &dest;
        int tmp_alloc=1,tmp_size=0;
        iMesh_getEntSets(impl_,parent,0,&tmp,&tmp_alloc,&tmp_size,&err);
        assert(!err);

        int ret = do_extrusion(src,dest,true,steps-1,trans);

        iMesh_destroyEntSet(impl_,parent,&err);
        assert(!err);
        iMesh_destroyEntSet(impl_,dest,&err);
        assert(!err);

        return ret;
    }
    else
        return do_extrusion(src,0,false,steps,trans);
}

int ExtrudeMesh::extrude(iBase_EntityHandle *src,iBase_EntityHandle *dest,
                         int size,int steps,const CopyVerts &trans)
{
    int err;
    iBase_EntitySetHandle src_set,dest_set;
    iMesh_createEntSet(impl_,false,&src_set,&err);
    assert(!err);
    iMesh_createEntSet(impl_,false,&dest_set,&err);
    assert(!err);
    
    iMesh_addEntArrToSet(impl_,src,size,src_set,&err);
    assert(!err);
    iMesh_addEntArrToSet(impl_,dest,size,dest_set,&err);
    assert(!err);

    int ret = extrude(src_set,dest_set,steps,trans);

    iMesh_destroyEntSet(impl_,src_set,&err);
    assert(!err);
    iMesh_destroyEntSet(impl_,dest_set,&err);
    assert(!err);

    return ret;
}

int ExtrudeMesh::extrude(iBase_EntitySetHandle src,iBase_EntitySetHandle dest,
                         int new_rows,const CopyVerts &trans)
{
    return do_extrusion(src,dest,true,new_rows-1,trans);
}

int ExtrudeMesh::do_extrusion(iBase_EntitySetHandle src,
                              iBase_EntitySetHandle dest,bool use_dest,
                              int new_rows,const CopyVerts &trans)
{
    assert(new_rows > 0 || use_dest);

    int err;

    iBase_EntityHandle *ents=0; int ent_alloc=0, ent_size=0;
    iBase_EntityHandle *adj=0;  int adj_alloc=0, adj_size=0;
    int *indices=0;             int ind_alloc=0, ind_size=0;
    int *offsets=0;             int off_alloc=0, off_size=0;

    iMesh_getAdjEntIndices(impl_,src,iBase_ALL_TYPES,iMesh_ALL_TOPOLOGIES,
                           iBase_VERTEX,
                           &ents,    &ent_alloc, &ent_size,
                           &adj,     &adj_alloc, &adj_size,
                           &indices, &ind_alloc, &ind_size,
                           &offsets, &off_alloc, &off_size,
                           &err);
    assert(!err);

    double dx[3];
    iBase_EntityHandle *curr;
    iBase_EntityHandle *next = adj;
    int *normals = 0;

    if(new_rows > 0)
    {
        int row_alloc = adj_size,row_size;
        curr = new iBase_EntityHandle[row_alloc];
        next = new iBase_EntityHandle[row_alloc];
        trans(1,adj,adj_size,&next,&row_alloc,&row_size);

        vtx_diff(dx,impl_,next[0],adj[0]);
        normals = get_normals(adj,indices,offsets,ent_size,dx);

        // Make the first set of volumes
        connect_the_dots(normals,indices,offsets,adj,
                         normals,indices,offsets,next,ent_size);

        // Make the inner volumes
        for(int i=2; i<=new_rows; i++)
        {
            std::swap(curr,next);
            trans(i,adj,adj_size,&next,&row_alloc,&row_size);
            connect_the_dots(normals,indices,offsets,curr,
                             normals,indices,offsets,next,ent_size);
        }
    }

    if(use_dest)
    {
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

        vtx_diff(dx,impl_,adj2[indices2[ offsets2[0] ]],
                          next[indices [ offsets [0] ]]);
        if(!normals)
            normals = get_normals(adj,indices,offsets,ent_size,dx);
        int *normals2 = get_normals(adj2,indices2,offsets2,ent2_size,dx);

        connect_the_dots(normals, indices, offsets, next,
                         normals2,indices2,offsets2,adj2,ent_size);

        free(normals2);
        free(ents2);
        free(adj2);
        free(indices2);
        free(offsets2);
    }

    if(new_rows > 0)
    {
        delete curr;
        delete next;
    }

    free(normals);
    free(ents);
    free(adj);
    free(indices);
    free(offsets);

    return 0;
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
        if(offsets[i+1] - offsets[i] > 2)
        {
            for(int j=0; j<3; j++)
                curr_verts[j] = verts[indices[ offsets[i]+j ]];

            iMesh_getVtxArrCoords(impl_,curr_verts,3,iBase_INTERLEAVED,
                                  &coords,&coord_alloc,&coord_size,&err);
            for(int j=0; j<3; j++)
            {
                a[j] = coords[1*3 + j] - coords[0*3 + j];
                b[j] = coords[2*3 + j] - coords[1*3 + j];
            }

            normals[i] = (dot( cross(res,a,b),dv ) > 0) ? 1:-1;
        }
        else if(offsets[i+1] - offsets[i] == 2)
        {
            normals[i] = 1; // TODO
        }
        else
            assert(false);
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
    using namespace std;

    for(int i=0; i<size; i++)
    {
        int count = pre_offs[i+1] - pre_offs[i];

        // If the normal is facing in the wrong direction (away from the
        // translation) we add the vertices in reverse order. Otherwise, we go
        // in the usual order. If count is 2, then we are creating quads and so
        // need to swap the order of the post set of verts.
        
        int dx = pre_norms [i];
        int dy = post_norms[i] * (count == 2 ? -1:1);
        int x  = (dx == 1) ? pre_offs [i] : pre_offs [i+1]-1;
        int y  = (dy == 1) ? post_offs[i] : post_offs[i+1]-1;

        iBase_EntityHandle *nodes = new iBase_EntityHandle[count*2];
        for(int j=0; j<count; j++)
        {
            nodes[j]       = pre [ pre_inds [x + dx*j] ];
            nodes[j+count] = post[ post_inds[y + dy*j] ];
        }

        int status;
        iBase_EntityHandle out;

        if(count == 4)      // quad
            iMesh_createEnt(impl_,iMesh_HEXAHEDRON,nodes,8,&out,&status,&err);
        else if(count == 3) // tri
            iMesh_createEnt(impl_,iMesh_PRISM,nodes,6,&out,&status,&err);
        else if(count == 2) // line
            iMesh_createEnt(impl_,iMesh_QUADRILATERAL,nodes,4,&out,&status,
                            &err);
        else
            assert(false);

        assert(status==0 && err==0);
        delete[] nodes;
    }
}

void test1()
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

    iMesh_createVtxArr(mesh,6,iBase_INTERLEAVED,verts,3*6,&ents,&size,&alloc,
                       &err);
    assert(err == 0);

    iBase_EntityHandle quad;
    int status;
    iMesh_createEnt(mesh,iMesh_QUADRILATERAL,ents,4,&quad,&status,&err);
    assert(err == 0);

    iBase_EntityHandle line;
    iMesh_createEnt(mesh,iMesh_LINE_SEGMENT,ents,2,&line,&status,&err);
    assert(err == 0);

    iBase_EntityHandle tri;
    iBase_EntityHandle tri_ents[] = { ents[1], ents[2], ents[4] };
    iMesh_createEnt(mesh,iMesh_TRIANGLE,tri_ents,3,&tri,&status,&err);
    assert(err == 0);

    iBase_EntityHandle faces[] = {quad, tri, line};
    double v[] = {0,0,1};
    int steps = 5;
    ext->translate(faces,3,steps,v);

    // VisIt doesn't like 1-d objects in pseudocolor volume graphs
    iMesh_deleteEnt(mesh,line,&err);
    assert(err == 0);

    int count;
    iMesh_getNumOfType(mesh,0,iBase_VERTEX,&count,&err);
    assert(err == 0 && count == 5*(steps+1)+1);

    iMesh_getNumOfType(mesh,0,iBase_FACE,&count,&err);
    assert(err == 0 && count == 2+1*steps);

    iMesh_getNumOfType(mesh,0,iBase_REGION,&count,&err);
    assert(err == 0 && count == steps*2);

    const char *file = "fun.vtk";
    iMesh_save(mesh,root_set,file,"",&err,strlen(file),0);
    assert(err == 0);

    delete ext;
    iMesh_dtor(mesh,&err);
}

void test2()
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

        0, 0, 1,
        1, 0, 1,
        1, 1, 1,
        0, 1, 1,
        2, 1, 1,
    };
    iBase_EntityHandle *ents=0;
    int size=0,alloc=0;

    iMesh_createVtxArr(mesh,11,iBase_INTERLEAVED,verts,3*11,&ents,&size,&alloc,
                       &err);
    assert(err == 0);

    iBase_EntityHandle quad[2];
    int status;
    iMesh_createEnt(mesh,iMesh_QUADRILATERAL,ents+0,4,quad+0,&status,&err);
    assert(err == 0);
    iMesh_createEnt(mesh,iMesh_QUADRILATERAL,ents+6,4,quad+1,&status,&err);
    assert(err == 0);

    iBase_EntityHandle tri[2];
    iBase_EntityHandle tri_ents[] = { ents[1], ents[2], ents[4],
                                      ents[7], ents[8], ents[10] };
    iMesh_createEnt(mesh,iMesh_TRIANGLE,tri_ents+0,3,tri+0,&status,&err);
    assert(err == 0);
    iMesh_createEnt(mesh,iMesh_TRIANGLE,tri_ents+3,3,tri+1,&status,&err);
    assert(err == 0);

    iBase_EntityHandle pre[]  = {quad[0], tri[0]};
    iBase_EntityHandle post[] = {quad[1], tri[1]};
    double v[] = {0,0,1};
    int steps = 5;
    ext->translate(pre,post,2,steps);

    int count;
    iMesh_getNumOfType(mesh,0,iBase_VERTEX,&count,&err);
    assert(err == 0 && count == 5*(steps+1)+1);

    iMesh_getNumOfType(mesh,0,iBase_FACE,&count,&err);
    assert(err == 0 && count == 4);

    iMesh_getNumOfType(mesh,0,iBase_REGION,&count,&err);
    assert(err == 0 && count == steps*2);

    const char *file = "fun2.vtk";
    iMesh_save(mesh,root_set,file,"",&err,strlen(file),0);
    assert(err == 0);

    delete ext;
    iMesh_dtor(mesh,&err);
}

void test3()
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
    };
    iBase_EntityHandle *ents=0;
    int size=0,alloc=0;

    iMesh_createVtxArr(mesh,4,iBase_INTERLEAVED,verts,3*4,&ents,&size,&alloc,
                       &err);
    assert(err == 0);

    iBase_EntityHandle quad;
    int status;
    iMesh_createEnt(mesh,iMesh_QUADRILATERAL,ents,4,&quad,&status,&err);
    assert(err == 0);

    iBase_EntityHandle faces[] = {quad};
    double v[] = {0,0,1};
    int steps = 5;
    ext->translate(faces,1,steps,v);

    int count;
    iMesh_getNumOfType(mesh,0,iBase_VERTEX,&count,&err);
    assert(err == 0 && count == 4*(steps+1));

    iMesh_getNumOfType(mesh,0,iBase_FACE,&count,&err);
    assert(err == 0 && count == 1);

    iMesh_getNumOfType(mesh,0,iBase_REGION,&count,&err);
    assert(err == 0 && count == steps);

    const char *file = "easy.vtk";
    iMesh_save(mesh,root_set,file,"",&err,strlen(file),0);
    assert(err == 0);

    delete ext;
    iMesh_dtor(mesh,&err);
}

#include <cmath>
void test4()
{
    int err;
    iMesh_Instance mesh;
    iMesh_newMesh("",&mesh,&err,0);
    iBase_EntitySetHandle root_set;
    iMesh_getRootSet(mesh,&root_set,&err);

    ExtrudeMesh *ext = new ExtrudeMesh(mesh);

    double verts[] = {
        0, 0, 0,
        0, 1, 0,
        0, 1, 1,
        0, 0, 1,
    };
    iBase_EntityHandle *ents=0;
    int size=0,alloc=0;

    iMesh_createVtxArr(mesh,4,iBase_INTERLEAVED,verts,3*4,&ents,&size,&alloc,
                       &err);
    assert(err == 0);

    iBase_EntityHandle quad;
    int status;
    iMesh_createEnt(mesh,iMesh_QUADRILATERAL,ents,4,&quad,&status,&err);
    assert(err == 0);

    iBase_EntityHandle faces[] = {quad};

    iBase_EntitySetHandle set;
    iMesh_createEntSet(mesh,false,&set,&err);
    assert(!err);
    
    iMesh_addEntArrToSet(mesh,faces,1,set,&err);
    assert(!err);

    int steps = 200;
    double origin[] = {0,-3,0};
    double z[] = {1,1,1};
    double angle = 2*3.14159/steps;
    ext->rotate(set,steps,origin,z,angle);

    int count;
    iMesh_getNumOfType(mesh,0,iBase_VERTEX,&count,&err);
    assert(err == 0 && count == 4*(steps+1));

    iMesh_getNumOfType(mesh,0,iBase_FACE,&count,&err);
    assert(err == 0 && count == 1);

    iMesh_getNumOfType(mesh,0,iBase_REGION,&count,&err);
    assert(err == 0 && count == steps);

    iMesh_destroyEntSet(mesh,set,&err);
    assert(!err);

    const char *file = "rotate.vtk";
    iMesh_save(mesh,root_set,file,"",&err,strlen(file),0);
    assert(err == 0);

    delete ext;
    iMesh_dtor(mesh,&err);
}

int main()
{
    using namespace std;
    cout << "Test 1...";
    test1();
    cout << " Passed!" << endl;

    cout << "Test 2...";
    test2();
    cout << " Passed!" << endl;

    cout << "Test 3...";
    test3();
    cout << " Passed!" << endl;

    cout << "Test 4...";
    test4();
    cout << " Passed!" << endl;
}
