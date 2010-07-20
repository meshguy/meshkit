/*********************************************
Mar,10
Reactor Geometry Generator
Argonne National Laboratory

A program to generate cubit journal files for 
generating nuclear reactor geometries specified 
in the input file
*********************************************/
#include "nrgen.hpp"
#include "clock.hpp"

int main (int argc, char *argv[])
{
  int err = 0;
  // the one and only NR!
  CNrgen TheNR; 

  // start the timer 
  CClock Timer;
  std::string szDateTime;
  Timer.GetDateTime (szDateTime);
  std::cout << "\nStarting out at : " << szDateTime << "\n";
  
  // show program banner
  TheNR.Banner (std::cout); 
  
  // prepare for I/O
  err =  TheNR.PrepareIO (argc, argv);
  ERRORR("Error in function PrepareIO", err);
 
  //count pin cylinders and cell material, needed for setting array size before actual read
  err = TheNR.ReadInputPhase1 ();
  ERRORR("Error in function ReadInputPhase1", err);

  // read the problem size and create pincell
  TheNR.ReadAndCreate ();
  ERRORR("Error in function ReadAndCreate", err);
  
  // create the .jou file
  TheNR.CreateCubitJournal();
  ERRORR("Error in function CreateCubitJournal", err);

  // get the current date and time
  Timer.GetDateTime (szDateTime);
  std::cout << "Ending at : " << szDateTime;
 
  // compute the elapsed time
  std::cout << "Elapsed wall clock time: " << Timer.DiffTime ()
	    << " seconds or " << (Timer.DiffTime ())/60.0 << " mins\n";

  // close input and output files
  TheNR.TerminateProgram ();
  ERRORR("Error in function TerminateProgram", err);

  std::cout << "\nxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"<<std::endl;
  
  return 0;
}
