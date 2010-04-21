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
  int err;
  CNrgen TheNR; // the one and only NR!

  // show program banner
  TheNR.Banner (std::cout); 

  // prepare for I/O
  TheNR.PrepareIO (argc, argv);
  
  // start the timer 
  CClock Timer;
  std::string szDateTime;
  Timer.GetDateTime (szDateTime);
  std::cout << "\nStarting out at : " << szDateTime << "\n";
  
  //count pin cylinders and cell material, needed for setting array size before actual read
  TheNR.CountPinCylinders ();
  
  // read the problem size and create pincell
  err = TheNR.ReadAndCreate ();
  if (err!=1)
    std::cout << "Error in function ReadAndCreate\n";
  
  // create the .jou file
  err = TheNR.CreateCubitJournal();
  if (err!=1)
    std::cout << "Error in function CreateCubitJournal\n";

  // get the current date and time
  Timer.GetDateTime (szDateTime);
  std::cout << "Ending at : " << szDateTime;
 
  // compute the elapsed time
  std::cout << "Elapsed wall clock time: " << Timer.DiffTime ()
	    << " seconds or " << (Timer.DiffTime ())/60.0 << " mins\n";

  // close input and output files
  err = TheNR.TerminateProgram ();
  if (err!=1)
    std::cout << "Error while terminating \n";
  std::cout << "\nxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"<<std::endl;
  
  return 0;
}
