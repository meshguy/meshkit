// Copyright (C) 2005, 2009 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id: hs071_main.cpp 1864 2010-12-22 19:21:02Z andreasw $
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2005-08-10


#include "IASolver.hpp"

#include <stdio.h>

int main(int argv, char* argc[])
{
  IASolver ia_solver;
  // setup
  bool status = ia_solver.solve();

  return (int) status;
}
