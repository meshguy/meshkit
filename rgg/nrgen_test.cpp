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
  
  CNrgen TheNR; // the one and only NR!

   // switch to set as 1 for surface only
  TheNR.m_nPlanar = 0;
  // show program banner
  TheNR.Banner (std::cout); 

  // Prepare for I/O
  TheNR.PrepareIO (argc, argv);
  
  // start the timer 
  CClock Timer;
  std::string szDateTime;
  Timer.GetDateTime (szDateTime);
  std::cout << "\nStarting out at : " << szDateTime << "\n";
  
  //count pin cylinders and cell material, needed for setting array size before actual read
  TheNR.CountPinCylinders ();
  
  // read the problem size and create pincell
  int i = TheNR.ReadAndCreate ();
  if (i!=1)
    std::cout << "Error in function ReadAndCreate";
  
  // get the current date and time
  Timer.GetDateTime (szDateTime);
  std::cout << "Ending at : " << szDateTime;
 
  // compute the elapsed time -----------------------------------------------
  std::cout << "Elapsed wall clock time: " << Timer.DiffTime ()
	    << " seconds or " << (Timer.DiffTime ())/60.0 << " mins\n";

  // Close input and output files
  int j = TheNR.TerminateProgram ();
  if (j!=1)
    std::cout << "Error ";
  std::cout << "\nxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"<<std::endl;
  
  return 0;
}
