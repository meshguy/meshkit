/*!
 * \class mstream
 * \brief Print something on console and file using redirection symbol <<,
 *        open a file in your program:  <mstream object>.coss.open (m_LogName.c_str(), std::ios::out);
 *
 **/
#ifndef __MESHKIT_MSTREAM_H__
#define __MESHKIT_MSTREAM_H__

#include <iostream>
#include <fstream>

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

#endif
