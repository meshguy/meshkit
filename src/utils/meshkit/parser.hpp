/*********************************************
Dec,09
Reactor Geometry Generator
Argonne National Laboratory

CParser class definition.
*********************************************/
#ifndef __RGG_PARSER_H__
#define __RGG_PARSER_H__

#include <fstream>
#include <string>
#include <vector>
#include <iostream>

class CParser
{
public:
  CParser ();
  ~CParser ();

  bool ReadNextLine (std::ifstream& FileInput, int& nL,
		     std::string& szInputString,
		     const int MAXCHARS, 
		     const std::string& szComment,
		     bool bLowerCase = true);
  void GetTokens (const std::string& input, const std::string& delims, 
		  std::vector<std::string>& tokens);
  void FilterComment(std::string& input, const std::string& szComment);
  void RemoveToken (std::string& input);

  bool EatWhiteSpace (std::string& input);
private:

};

#endif
