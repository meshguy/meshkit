#include <iostream>
#include <fstream>
/*!
 * \class mstream
 * \brief Print something on console and file using redirection symbol <<
 *
 **/
class mstream
{
  public:
  std::ofstream coss;
  mstream();
  ~mstream();
  mstream& operator<< (std::ostream& (*pfun)(std::ostream&));
};

template <class T>
mstream& operator<< (mstream& st, T val)
{
  st.coss << val;
  std::cout << val;
  return st;
};
