/*
 * QslimOptions.hpp
 *
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

//===========================================================================//
  /*!
   * \class QslimOptions
   * \brief Place holder for decimation options

   *  Options for decimation can be passed with the test driver.
   *
   * -o \<file\>      Output final model to given file, using moab output conventions
   *
   * -s \<count\>      Set the target number of faces.
   *
   * -i \<timings\>   Intervals for timing reports.
   *
   * -e \<thresh\> Set the maximum error tolerance.
   *
   * -t \<t\>        Set pair selection tolerance.
   *
   * -Q[pv]        Select what constraint quadrics to use [default=p].
   *
   * -On       Optimal placement policy.
                 0=endpoints, 1=endormid, 2=line, 3=optimal [default]

   *  -B \<weight\> Use boundary preservation planes with given weight.
   *
   *  -b       preserve boundary (do not use with -B option)
   *
   *  -m            Preserve mesh quality.
   *
   *  -a            Enable area weighting.
   *
   *  -p          Height fields positivity. Used for height fields, assume
                       triangles are originally positively oriented.

   *  -d          Use delayed deletion, as opposed to merging
   *
   *  -c           keep costs in a (sparse!!!!) tag
   *

 */
class QslimOptions {
public:
	QslimOptions();
	virtual ~QslimOptions();
	//! decimation stops when number of triangle is less than face_target
	int face_target;

	//! decimation can stop if the error is bigger than the given tolerance
	double error_tolerance;

	//
	int     will_use_plane_constraint;
	int     will_use_vertex_constraint;

	int     will_preserve_boundaries;
	int     will_preserve_mesh_quality;
	int     will_constrain_boundaries;
	double boundary_constraint_weight;

	int     will_weight_by_area;
	int     height_fields;
	int     timingIntervals;
	int     plotCost;

	int useDelayedDeletion;

	int placement_policy;
	double pair_selection_tolerance;
	std::ostream * logfile;
	int selected_output;
};

#endif /* QSLIMOPTIONS_H_ */
