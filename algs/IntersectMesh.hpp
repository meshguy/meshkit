/*
 * IntersectMesh.hpp
 *
 *  Created on: Aug 13, 2010
 *      Author: iulian
 *
 *  When 2 meshes are discretizing the same domain, their intersection
 *  could be of interest; it would help in coupling the solutions for 2 different
 *  physics problem.
 *  This case solves 2D mesh triangular intersections;
 *  it should be extended to 3d someday
 *
 *  Input : 2 triangle sets in 2 different meshes, and one pair of triangles
 *  in different meshes that intersect
 *  Output: a new mesh instance, and a set of triangles. Each triangle is completely
 *  included in a triangle in each of original meshes
 *
 */

#ifndef INTERSECTMESH_HPP_
#define INTERSECTMESH_HPP_

#include "iBase.h"
#include "iMesh.h"

class IntersectMesh {
public:
	IntersectMesh( iMesh_Instance mesh1,  iBase_EntitySetHandle set1,
			iMesh_Instance mesh2,  iBase_EntitySetHandle set2);
	virtual ~IntersectMesh();

	// from a seed, compute the triangles from intersection
	int compute(iBase_EntityHandle t1, iBase_EntityHandle t2 ,
			iMesh_Instance outMesh, iBase_EntitySetHandle set);
private:
	iMesh_Instance _mesh1;
	iBase_EntitySetHandle _set1;
	iMesh_Instance _mesh2;
	iBase_EntitySetHandle _set2;
};

#endif /* INTERSECTMESH_HPP_ */
