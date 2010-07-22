/*
 * SmoothBase.cpp
 *
 *  Created on: Jul 19, 2010
 *      Author: iulian
 */

#include "SmoothBase.hpp"

SmoothBase::SmoothBase(MBInterface * mbIn, MBEntityHandle setIn, MBInterface * mbOut):
      _mb(mbIn), _set(setIn), _mbOut(mbOut) {

	// create new entity set _oSet, which will store the new quad mesh
	//_mbOut->create_meshset(MESHSET_SET, _oSet);
	// TODO Auto-generated constructor stub

}

SmoothBase::~SmoothBase() {
	// TODO Auto-generated destructor stub
}

