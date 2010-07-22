/*
 * SmoothVertex.hpp
 *
 *  Created on: Jul 19, 2010
 *      Author: iulian
 */

#include "MBInterface.hpp"
#include "SmoothBase.hpp"

#ifndef SMOOTHVERTEX_HPP_
#define SMOOTHVERTEX_HPP_

class SmoothVertex : public SmoothBase
{
public:
	SmoothVertex(MBInterface * mb, MBEntityHandle vset, MBInterface * mbo);
	virtual ~SmoothVertex();

	MBErrorCode create_mesh_vertex();

	MBEntityHandle get_new_node();
private :

};

#endif /* SMOOTHVERTEX_HPP_ */
