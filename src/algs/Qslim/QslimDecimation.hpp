/*
 * QslimDecimation.hpp
 *
 */

#ifndef QSLIMDECIMATION_H_
#define QSLIMDECIMATION_H_

#include "moab/Interface.hpp"
#include "moab/Range.hpp"
#include "meshkit/QslimOptions.hpp"
#include "std.h"
//#include "AdjModel.h"

namespace MeshKit {

class QslimDecimation {
public:
	QslimDecimation (moab::Interface * mb, moab::EntityHandle root_set);
	virtual ~QslimDecimation();

int decimate(QslimOptions & opts, moab::Range & oRange);

private:
 moab::Interface * _mb;
 moab::EntityHandle _InitialSet;

 int Init ();
// Model * m_model;
};

}  // namespace MeshKit

#endif /* QSLIMDECIMATION_H_ */
