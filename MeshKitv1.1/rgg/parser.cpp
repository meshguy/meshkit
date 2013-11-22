/*********************************************
Reactor Geometry Generator
Argonne National Laboratory

Contains CParser class implementation.
*********************************************/
#include <cctype>
#include "parser.hpp"
//#include <math.h>
#include <stdlib.h>

CParser::CParser ()
// ---------------------------------------------------------------------------
// Function: default constructor
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
}

CParser::~CParser ()
// ---------------------------------------------------------------------------
// Function: destructor
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
}

bool CParser::ReadNextLine (std::ifstream& FileInput, int& nLineNum, 
			    std::string& szInputString, const int MAXCHARS,
			    const std::string& szComment, bool bLowerCase)
// ---------------------------------------------------------------------------
// Function: reads the next line skipping over the comment lines
//           and converts all alphabets to lower case if requested
// Input:    file istream, line #, string to hold the input line,
//           max. # of characters expected in each input line,
//           comment character(s) at the beginning of a comment line,
//           lowercase conversion option
// Output:   updated values of line # and the string
//           return value is true if successful
//                           false if an error state is encountered
// Restriction: Cannot read a line over 256 characters
// ---------------------------------------------------------------------------
{
  int flag = 0;
  int flag1 =0;
  bool bWhSpc = false;
  int tokenfound = 1;
  const int MAXCH = 1500;
  char szInp[MAXCH];
  char szTemp [MAXCH];
  std::vector<std::string> tokens;

  // enough capacity to read and store?
  if (MAXCHARS > MAXCH)
    return false;

  // comment character(s)
  int nCLen = static_cast<int>(szComment.length());
  // read the line (skip over comment lines)
  for(;;)
    {
      ++nLineNum;
      FileInput.getline (szInp, MAXCHARS);
//       // end-of-file?
//       if (FileInput.eof())
// 	return false;
      if (FileInput.fail())
	FileInput.clear (FileInput.rdstate() & ~std::ios::failbit);
      // unrecoverable error?
      if (FileInput.bad())
	return false;

      // successful read
      szInputString = szInp;
      GetTokens(szInputString, " ", tokens);
      bWhSpc = EatWhiteSpace(szInputString);
      if ((szInputString.substr(0,nCLen) != szComment)&& (bWhSpc ==false)){	   
	szInputString = szInp;
	GetTokens(szInputString, " ", tokens);
	for(int i=0; i< abs(tokens.size()); i++){
	  std::string temptoken = tokens[i];
	  if (temptoken == "&")
	    flag1 = 1;
	}

	//Filter the comment tokens
	//  FilterComment(szInputString, szComment);

	//if "&" is found continue to read the next line      
	std::string szTempString = szInputString;

	// check if line is continued &
	while(flag1 ==1 && tokenfound == 1){	
	  GetTokens(szTempString, " ", tokens);
	  for(int i=1; i<=abs(tokens.size()); i++){
	    std::string temptoken = tokens[i-1];
	    if (temptoken == "&"){
	      tokenfound = 1;
	      flag = 1;
	    }
	    else{
	      if(flag==1)
		flag = 1;//do nothing token already found
	      else
		tokenfound = 0;
	    }
	  }
	  if(tokenfound ==1){
	    ++nLineNum;
	    RemoveToken(szInputString);
	    //- getting more tokens and add to the existing 
	    FileInput.getline (szTemp, MAXCHARS);
	    // end-of-file?
	    if (FileInput.eof())
	      return false;
	    if (FileInput.fail())
	      FileInput.clear (FileInput.rdstate() & ~std::ios::failbit);
	    // unrecoverable error?
	    if (FileInput.bad())
	      return false;
	    // successful read 
	    szTempString = szTemp;
	    FilterComment(szTempString, szComment);
	    szInputString+=" ";
	    szInputString+=szTemp;
	  }
	  else{
	    break;//while loop ents
	  }
	  flag = 0;
	}
	// while loop ends
	// convert to lower case?
	if (bLowerCase){
	  for (int i=0; i < static_cast<int>(szInputString.length()); i++)
	    szInputString[i] = tolower(szInputString[i]);
	}
	break;
      }
    }
  return true;
}

void CParser::GetTokens (const std::string& input, 
			 const std::string& delims, 
			 std::vector<std::string>& tokens)
// ----------------------------------------------------------------------------
// Function: Parses the input line and gets the tokens
// Input:    string, delimiters
// Output:   vector containing the tokens
// ----------------------------------------------------------------------------
{
  std::string::size_type beg_index, end_index;

  // clear the vector that will store the tokens
  tokens.clear();

  // get location of the first character that is not a delimiter
  beg_index = input.find_first_not_of(delims);

  // loop while the beginning index is not the end of string
  while (beg_index != std::string::npos)
    {
      // get location of the next delimiter
      end_index = input.find_first_of (delims, beg_index);

      // if this location is the end of string then adjust the value
      // as string length
      if (end_index == std::string::npos) end_index = input.length();

      // save the string between beg_index and end_index
      tokens.push_back (input.substr(beg_index,end_index-beg_index));

      // get location of the next character that is not a delimiter
      beg_index = input.find_first_not_of (delims, end_index);
    }
}


void CParser:: FilterComment (std::string& input,  const std::string& szComment)

// ----------------------------------------------------------------------------
// Function: Parses the input line and gets the tokens
// Input:    string, delimiters
// Output:   vector containing the tokens
// ----------------------------------------------------------------------------
{
  // remove comment from the line obtained
  int i;
  std::vector<std::string> tokens;
  std::string tempInput;
  GetTokens(input, " ", tokens);
  for(i=0; i<abs(tokens.size()); i++){
    std::string temptoken = tokens[i];
    if(temptoken == szComment){
      break;
    }
    else{
      tempInput+=temptoken;
      if(i!=(abs(tokens.size())-1)){ //indent{
	if (tokens[i+1]!=szComment)
	  tempInput+=" ";
      }
      else{
	tempInput+=""; // no indent
      }
    }
  }
  input = tempInput.c_str();
}

void CParser:: RemoveToken (std::string& input)

// ----------------------------------------------------------------------------
// Function: Parses the input line and gets the tokens
// Input:    string, delimiters
// Output:   vector containing the tokens
// ----------------------------------------------------------------------------
{
  // remove comment from the line obtained
  int i;
  std::vector<std::string> tokens;
  std::string tempInput;
  GetTokens(input, " ", tokens);
  for(i=0; i<abs(tokens.size()); i++){
    std::string temptoken = tokens[i];
    if(temptoken == "&"){
      break;
    }
    else{
      tempInput+=temptoken;
      if(i!=(abs(tokens.size())-1)){ //indent{
	if (tokens[i+1]!="&")
	  tempInput+=" ";
      }
      else{
	tempInput+=""; // no indent
      }
    }
  }
  input = tempInput.c_str();
}

bool CParser:: EatWhiteSpace (std::string& input)

// ----------------------------------------------------------------------------
// Function: Parses the input line and gets the tokens
// Input:    string, delimiters
// Output:   vector containing the tokens
// ----------------------------------------------------------------------------
{
  // remove comment from the line obtained
  int i;
  std::vector<std::string> tokens;
  std::string tempInput;
  GetTokens(input, " ", tokens);
  for(i=0; i<abs(tokens.size()); i++){
    std::string temptoken = tokens[i];
    if(temptoken != " " && temptoken !="!"){
      return false;
      break;  
    }
    else if(temptoken =="!"){
      return true;
      break;
    }
  }
  return false;
}
