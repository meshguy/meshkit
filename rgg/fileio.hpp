/*********************************************
Dec,09
Reactor Geometry Generator
Argonne National Laboratory

File I/O
*********************************************/
#ifndef __RGG_FILEIO_H__
#define __RGG_FILEIO_H__

#include <iostream>
#include <fstream>	
#include <string>
#include <vector>

void OpenInputFileByName (const std::string& szPrompt, 
                          std::ifstream& IFile, 
                          const std::ios::openmode&);
void OpenOutputFileByName (const std::string& szPrompt,
                           std::ofstream& OFile,
                           const std::ios::openmode&);
void Rewind (std::ifstream& IOFile);

#endif
