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
#include <vector>

extern int validFaceCount; // defined in QslimDecimation

extern int classifyVertex(MBInterface * mb, MBEntityHandle v);

extern MBErrorCode contractionRegion(MBInterface * mb, MBEntityHandle v1, MBEntityHandle v2,
		std::vector<MBEntityHandle> & changed);

extern Vec3 getVec3FromMBVertex(MBInterface * mb, MBEntityHandle v);

extern Plane trianglePlane(MBInterface * mb, MBEntityHandle tri);

extern MBErrorCode contract (MBInterface * mb, MBEntityHandle v0, MBEntityHandle v1, Vec3 & vnew, std::vector<MBEntityHandle> & changed );

#endif /* PRIMITIVES_H_ */
