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
  int run_flag = 1; 

  // the one and only Core!
  CCrgen TheCore; 
 
  err = TheCore.banner ();
  ERRORR("Failed in creating banner", 1);

  // start the timer 
  CClock Timer;
  std::string szDateTime;
  Timer.GetDateTime (szDateTime);
  std::cout << "\nStarting out at : " << szDateTime << "\n";
  
  // read inputs and create makefile
  err = TheCore.prepareIO (argc, argv);
  ERRORR("Failed in preparing i/o.", 1);

  if (argc > 1){
    if(argv[1][0] == '-' && argv[1][1] == 'm'){
      run_flag = 0;
    }
  }

  // copy, move, merge, extrude, assign gids, save and close
  if (run_flag == 1){
    err = TheCore.load_meshes ();
    ERRORR("Failed to load meshes.", 1);

    err = TheCore.copy_move ();
    ERRORR("Failed in copy move routine.", 1);

    err = TheCore.merge_nodes ();
    ERRORR("Failed to merge nodes.", 1);

    err = TheCore.extrude();
    ERRORR("Failed to extrude.", 1);

    err = TheCore.assign_gids ();
    ERRORR("Failed to assign global ids.", 1);

    err = TheCore.save ();
    ERRORR("Failed to save o/p file.", 1);
    
    err = TheCore.close ();
    ERRORR("Failed to dellocate.", 1);
  }

  // get the current date and time
  Timer.GetDateTime (szDateTime);
  std::cout << "Ending at : " << szDateTime;
 
  // compute the elapsed time
  std::cout << "Elapsed wall clock time: " << Timer.DiffTime ()
	    << " seconds or " << (Timer.DiffTime ())/60.0 << " mins\n";

  return 0;
}
