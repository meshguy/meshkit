#ifndef COPYVERTS_HPP
#define COPYVERTS_HPP

#include "iBase.h"
#include "iMesh.h"

class CopyVerts
{
public:
    CopyVerts(iMesh_Instance impl);

    void operator ()(int n,iBase_EntityHandle *src,int src_size,
                     iBase_EntityHandle **dest,int *dest_alloc,int *dest_size)
        const;
protected:
    virtual void transform(int n,int i,double *coords) const = 0;
private:
    iMesh_Instance impl_;
};

class CopyMoveVerts : public CopyVerts
{
public:
    CopyMoveVerts(iMesh_Instance impl,double *dv);
protected:
    virtual void transform(int n,int i,double *coords) const;
private:
    double dv_[3];
};

class CopyRotateVerts : public CopyVerts
{
public:
    CopyRotateVerts(iMesh_Instance impl,double *origin,double *z,double theta);
protected:
    virtual void transform(int n,int i,double *coords) const;
private:
    double origin_[3],z_[3];
    double theta_;
};

#endif
