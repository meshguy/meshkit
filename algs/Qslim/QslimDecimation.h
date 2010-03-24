/*
 * QslimDecimation.h
 *
 *  Created on: Mar 10, 2010
 *      Author: iulian
 */

#ifndef QSLIMDECIMATION_H_
#define QSLIMDECIMATION_H_


#include "iMesh.h"
#include "QslimOptions.h"
#include "std.h"
//#include "AdjModel.h"

class QslimDecimation {
public:
	QslimDecimation (iMesh_Instance mesh, iBase_EntitySetHandle root_set);
	virtual ~QslimDecimation();

int decimate(QslimOptions & opts);

private:
 iMesh_Instance m_mesh;
 iBase_EntitySetHandle m_InitialSet;

 int Init ();
// Model * m_model;
};

#endif /* QSLIMDECIMATION_H_ */
