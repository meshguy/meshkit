/*
 * SmoothBase.hpp
 *
 *  Created on: Jul 19, 2010
 *      Author: iulian
 */

#ifndef MESHKIT_SMOOTHBASE_HPP
#define MESHKIT_SMOOTHBASE_HPP

#include "MBInterface.hpp"

class SmoothBase {
public:
	SmoothBase(MBInterface * mbIn, MBEntityHandle setIn, MBInterface * mbOut);
	SmoothBase() { };// this one should not be called, actually
	virtual ~SmoothBase();


	MBEntityHandle out_set()  { return _oSet; }
protected:
	MBInterface * _mb; // in
	MBEntityHandle _set; // in
	MBInterface * _mbOut; // output mesh instance
	MBEntityHandle _oSet; // output set, which will contain the mesh entities; it will be a mesh set for face and vertex
	// and an ordered set for curve
};

#endif /* SMOOTHBASE_HPP_ */
