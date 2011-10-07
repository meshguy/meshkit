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

 int Init ();
};

}  // namespace MeshKit

#endif /* QSLIMDECIMATION_H_ */
