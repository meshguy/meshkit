/*********************************************
June,10
Reactor Mesh Assembler
Argonne National Laboratory

Driver program to generate core mesh by 
copy/move/merge'ing assembly meshes 
as specified in the input file
*********************************************/
#include "crgen.hpp"
#include "clock.hpp"

int main (int argc, char *argv[])
{
  int err = 0;

  // the one and only Core!
  CCrgen TheCore; 
 
  err = TheCore.banner ();
  ERRORR("Failed in creating banner", 1);

  // start the timer 
  CClock Timer;
  std::string szDateTime;
  Timer.GetDateTime (szDateTime);
  std::cout << "\nStarting out at : " << szDateTime << "\n";
  
  // read inputs, create and write makefile 
  err = TheCore.prepareIO (argc, argv);
  ERRORR("Failed in preparing i/o.", 1);

  // read inputs and open makefile for writing 
  err = TheCore.load_meshes ();
  ERRORR("Failed in loading meshes.", 1);

  err = TheCore.copy_move ();
  ERRORR("Failed in copy moving meshes.", 1);

  err = TheCore.merge_nodes ();
  ERRORR("Failed in merging nodes.", 1);

  err = TheCore.assign_gids ();
  ERRORR("Failed in assigning global ids.", 1);

  err = TheCore.save ();
  ERRORR("Failed in saving file.", 1);

  // get the current date and time
  Timer.GetDateTime (szDateTime);
  std::cout << "Ending at : " << szDateTime;
 
  // compute the elapsed time
  std::cout << "Elapsed wall clock time: " << Timer.DiffTime ()
	    << " seconds or " << (Timer.DiffTime ())/60.0 << " mins\n";

  return 0;
}
