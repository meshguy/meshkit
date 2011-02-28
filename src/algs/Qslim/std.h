#ifndef GFX_STD_INCLUDED // -*- C++ -*-
#define GFX_STD_INCLUDED

/************************************************************************

  Standard base include file for all gfx-based programs.  This defines
  various common stuff that is used elsewhere.

  $Id: std.h,v 1.8 1997/06/25 14:12:31 garland Exp $

 ************************************************************************/


#include <stdlib.h>
#include <string>
#include <math.h>
#include <iostream>


#ifndef FEQ_EPS
#define FEQ_EPS 1e-6
#define FEQ_EPS2 1e-12
#endif
inline bool FEQ(double a,double b,double eps=FEQ_EPS) { return fabs(a-b)<eps; }
inline bool FEQ(float a,float b,float eps=FEQ_EPS) { return fabsf(a-b)<eps; }

#ifndef GFX_NO_AXIS_NAMES
enum Axis {X=0, Y=1, Z=2, W=3};
#endif


#define fatal_error(s)  report_error(s,__FILE__,__LINE__)

#ifdef assert
#  undef assert
#endif
#define  assert(i)  (i)?((void)NULL):assert_failed(# i,__FILE__,__LINE__)

#ifdef SAFETY
#  define AssertBound(t) assert(t)
#else
#  define AssertBound(t)
#endif

//
// Define the report_error and assert_failed functions.
//
inline void report_error(char *msg,char *file,int line)
{
    std::cerr << msg << " ("<<file<<":"<<line<<")"<<std::endl;
    exit(1);
}

inline void assert_failed(char *text,char *file,int line)
{
	std::cerr << "Assertion failed: {" << text <<"} at ";
	std::cerr << file << ":" << line << std::endl;
    abort();
}

inline void assertf(int test, char *msg,
                    char *file=__FILE__, int line=__LINE__)
{
    if( !test )
    {
    	std::cerr << "Assertion failed: " << msg << " at ";
    	std::cerr << file << ":" << line << std::endl;
        abort();
    }
}

// GFX_STD_INCLUDED
#endif
