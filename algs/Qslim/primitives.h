/*
 * primitives.h
 *
 *  Created on: Mar 22, 2010
 *      Author: iulian
 */

#ifndef PRIMITIVES_H_
#define PRIMITIVES_H_

#include "std.h"
#include "3D.h"
#include "Mat4.h"
#include "MBInterface.hpp"
#include "QslimOptions.h"
#include <vector>

extern int validFaceCount; // defined in QslimDecimation
extern MBTag validTag;   // defined in QslimDecimation
// the plane data will be stored for each triangle, and recomputed for each triangle
// it will be needed only for the -m option active (so do not recompute if you do not need it)
extern MBTag planeDataTag;  // defined in QslimDecimation
extern int ehIsValid(MBEntityHandle v);  // defined in QslimDecimation; maybe we should pass Interface too
extern QslimOptions opts;

extern void filterValid(MBInterface * mb, std::vector<MBEntityHandle> & io);
extern int classifyVertex(MBInterface * mb, MBEntityHandle v);

extern MBErrorCode contractionRegion(MBInterface * mb, MBEntityHandle v1, MBEntityHandle v2,
		std::vector<MBEntityHandle> & changed);

extern Vec3 getVec3FromMBVertex(MBInterface * mb, MBEntityHandle v);

// this will just retrieve it from the tag data
extern Plane trianglePlane(MBInterface * mb, MBEntityHandle tri);

extern void computeTrianglePlane (MBInterface * mb, MBEntityHandle tri);

extern MBErrorCode contract (MBInterface * mb, MBEntityHandle v0, MBEntityHandle v1, Vec3 & vnew, std::vector<MBEntityHandle> & changed );

#endif /* PRIMITIVES_H_ */
