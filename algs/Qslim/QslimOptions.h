/*
 * QslimOptions.h
 *
 *  Created on: Mar 10, 2010
 *      Author: iulian
 */

#ifndef QSLIMOPTIONS_H_
#define QSLIMOPTIONS_H_

#define PLACE_ENDPOINTS 0
#define PLACE_ENDORMID  1
#define PLACE_LINE      2
#define PLACE_OPTIMAL   3

#define OUTPUT_NONE         0x00
#define OUTPUT_CONTRACTIONS 0x01
#define OUTPUT_QUADRICS     0x02
#define OUTPUT_COST         0x04
#define OUTPUT_VERT_NOTES   0x08
#define OUTPUT_FACE_NOTES   0x10
#define OUTPUT_MODEL_DEFN   0x20
#define OUTPUT_ALL          0xFFFFFFFF

#include <fstream>

class QslimOptions {
public:
	QslimOptions();
	virtual ~QslimOptions();
	int face_target;
	double error_tolerance;


	int     will_use_plane_constraint;
	int     will_use_vertex_constraint;

	int     will_preserve_boundaries;
	int     will_preserve_mesh_quality;
	int     will_constrain_boundaries;
	double boundary_constraint_weight;

	int     will_weight_by_area;
	int     height_fields;
	int     timingIntervals;

	int useDelayedDeletion;

	int placement_policy;
	double pair_selection_tolerance;
	std::ostream * logfile;
	int selected_output;
};

#endif /* QSLIMOPTIONS_H_ */
