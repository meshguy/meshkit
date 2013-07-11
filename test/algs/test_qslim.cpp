/*! \file test_qslim.cpp \test
 *
 * test_qslim.cpp
 *
 */

#include <iostream>

#include <time.h>
#include <stdlib.h>
#include <cstring>

#include "meshkit/MKCore.hpp"
#include "meshkit/MeshOp.hpp"
#include "meshkit/ModelEnt.hpp"
#include "meshkit/QslimMesher.hpp"
#include "meshkit/QslimOptions.hpp"

using namespace MeshKit;

#include "TestUtil.hpp"

MKCore *mk;

std::string usage_string =
"-o <file>	    Output final model to given file.\n"
"-s <count>	    Set the target number of faces.\n"
"-i <timings>   Intervals for timing reports.\n"
"-e <thresh>	Set the maximum error tolerance.\n"
"-t <t>		    Set pair selection tolerance.\n"
"-Q[pv]		    Select what constraint quadrics to use [default=p].\n"
"-On		    Optimal placement policy.\n"
"			    0=endpoints, 1=endormid, 2=line, 3=optimal [default]\n"
"-B <weight>	Use boundary preservation planes with given weight.\n"
"-b             preserve boundary (do not use with -B option) \n"
"-m		        Preserve mesh quality.\n"
"-a		        Enable area weighting.\n"
"-p 	        Height fields positivity. Used for height fields, assume \n"
"                   triangles are originally positively oriented. \n"
"-d 	        do not use delayed deletion \n"
"-c           keep costs in a (sparse!!!!) tag \n"
"-r           create the range with resulting triangles, and delete the original elements"
"-k           keep the topology unchanged, like holes, tunnels."
"\n";

std::string logging_usage_string =
"-l <file>	Log all simplification operations to given file.\n"
"-L[xqcvfdA]	Set information to be output.\n"
"			x=contractions, q=quadrics, c=cost\n"
"			v=vert. notes, f=face notes\n"
"			d=model defn, A=All.\n";

int main(int argc, char* argv[])
{
  // check command line arg
  const char *filename = 0;
  const char *outfile = 0;

  //process options
  QslimOptions options;
  // process all options, as given with regular qslim
  // the last argument must be the input file
  // one argument must be -o <the output file>, otherwise it is the standard output
  // in the process, a mesh file is read (instanced) and the
  // set of triangles is passed
  // it could be the root set for the mesh
  std::string fstr;
  if (argc<=1)
  {
	  std::cout<<usage_string;
	  std::cout << "\n\n";
	  fstr=TestDir + "/partBed.smf";
	  std::cout<< "default arguments: -s 4500 -B 1000 -p -o out.smf " << fstr <<" \n";
	  options.face_target = 4500;
	  options.will_constrain_boundaries = true;
	  options.boundary_constraint_weight = 1000;
	  options.height_fields = 1;

	  filename = fstr.c_str();
	  // ostr = "out.smf";
	  // outfile = ostr.c_str();
  }
  else
  {
	  // the file should be the last argument; if not, we have a problem, we will not be able to read it
	  filename = argv[argc-1];
	  int i=1;// will loop through arguments, and process them
	  for (i=1; i<argc-1 ; i++)
	  {
		  if (argv[i][0]=='-')
		  {
			  switch (argv[i][1])
			  {
				  case 'l':
				  {
					  // open a log file
					  options.logfile = new std::ofstream(argv[i+1]);
					  if (!options.logfile)
					  {
						  // cannot open log file, exit
						  std::cout << "can't open log file\n";
						  return 1;
					  }
					  i++;
					  break;
				  }
				  case 'L':
					  // parse the logging options
				  {
					  char *c;
					  int errflg = 0;
					  c = argv[i]+2;// skip -L
					  while( *c )
						{
							if( *c=='x' )
								options.selected_output |= OUTPUT_CONTRACTIONS;
							else if( *c=='q' )
								options.selected_output |= OUTPUT_QUADRICS;
							else if( *c=='c' )
								options.selected_output |= OUTPUT_COST;
							else if( *c=='v' )
								options.selected_output |= OUTPUT_VERT_NOTES;
							else if( *c=='f' )
								options.selected_output |= OUTPUT_FACE_NOTES;
							else if( *c=='d' )
								options.selected_output |= OUTPUT_MODEL_DEFN;
							else if( *c=='A' )
								options.selected_output |= OUTPUT_ALL;
							else
							errflg++;
							c++;
						}
					  if (errflg>0)
					  {
						  std::cout<<logging_usage_string;
						  std::cout<< "!! Ignore the error it\n";

					  }
					  break;
				  }
				  case 'a':
				  {
					  options.will_weight_by_area = 1;
					  break;
				  }
				  case 'B':
				  {
					  options.boundary_constraint_weight = atof(argv[i+1]);
					  options.will_constrain_boundaries = 1;
					  i++;
					  break;
				  }
				  case 's':
				  {
					  options.face_target = atoi(argv[i+1]);
					  i++;
					  break;
				  }
				  case 't':
				  {
					  options.pair_selection_tolerance = atof(argv[i+1]);
					  i++;
					  break;
				  }
				  case 'o':
				  {
					  outfile = argv[i+1];
					  i++;
					  break;
				  }

				  case 'O':
				  {
					  char * c ;
					  c = argv[i]+2;
					  int ival = atoi(c);
					  if( ival < 0 || ival > PLACE_OPTIMAL )
						  ival=PLACE_OPTIMAL;
					  options.placement_policy = ival;
					  break;
				  }

				  case 'Q':
				  {
					  options.will_use_plane_constraint = false;
					  options.will_use_vertex_constraint = false;
					  options.will_use_plane_constraint = false;

					  char * c ;
					  c = argv[i]+2;
					  while( *c )
					  {
						 if( *c=='p' )
							options.will_use_plane_constraint = true;
						 else if( *c=='v' )
							options.will_use_vertex_constraint = true;
						 else
						 {
							std::cout<< "wrong input " << argv[i] << std::endl;
						 }
						 c++;
					  }
					  break;
				  }

				  case 'e':
				  {
					 options.error_tolerance = atof(argv[i+1]);
					 i++;
					 break;
				  }

				  case 'b':
				  {
					 options.will_preserve_boundaries = true;
					 break;
				  }

				  case 'm':
				  {
					 options.will_preserve_mesh_quality = true;
					 break;
				  }

				  case 'p':
				  {
					 options.height_fields = true;
					 break;
				  }
				  case 'i':
				  {
					  options.timingIntervals=atof(argv[i+1]);
					  i++;
					  break;
				  }
				  case 'd':
				  {
					  options.useDelayedDeletion = false;
					  break;
				  }
				  case 'c':
				  {
				     options.plotCost = 1;
				     break;
				  }
				  case 'r':
				  {
				    options.create_range = true;
				    break;
				  }
				  case 'k':
				  {
				    options.topology = 1;
				    break;
				  }
				  default :
				  {
					 std::cout << "unsupported wrong input argument " << argv[i] <<std::endl;
					 return 0;
				  }
			  }
		  }
	  }
  }
  // initialize everything

  mk = new MKCore();
  mk->load_mesh(filename);
  MEntVector selection, dum;
  mk->get_entities_by_dimension(2, dum);
  selection.push_back(*dum.rbegin());// push just the last one retrieved from core

  QslimMesher *qm = (QslimMesher*) mk->construct_meshop("QslimMesher", selection);

  qm->set_options(options);

  mk->setup_and_execute();

  if(outfile)
  {
    std::cout << "writing output to " << outfile << std::endl;
    mk->moab_instance()->write_mesh(outfile);// write everything left
  }

  return 0;
}
