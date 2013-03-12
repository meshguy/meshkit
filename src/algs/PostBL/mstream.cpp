#include "mstream.hpp"
mstream::mstream ()
// ---------------------------------------------------------------------------
// Function: default constructor
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
}

mstream::~mstream ()
// ---------------------------------------------------------------------------
// Function: destructor
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
}

mstream& mstream::operator<< (std::ostream& (*pfun)(std::ostream&))
// ---------------------------------------------------------------------------
// Function: helps enable endline - std::endl
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
  pfun(coss);
  pfun(std::cout);
  return *this;
}
