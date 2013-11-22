/**\file test_testutil.cpp
 *\brief Make sure test suite utilities work (a meta-test?)
 */
 
#include "TestUtil.hpp"
#include <iostream>

void test_fail()      { CHECK(false); }
void test_abort()     { abort();      }
void test_exception() { throw 5;      }
void test_success()   {               }

int main()
{
  int fail = RUN_TEST(test_fail);
  int abrt = RUN_TEST(test_abort);
  int expt = RUN_TEST(test_exception);
  int succ = RUN_TEST(test_success);
  
  if (!fail) 
    std::cerr << "TestSuite did not register failed test!" << std::endl;
  if (!abrt) 
    std::cerr << "TestSuite did not regsiter aborted test!" << std::endl;
  if (!abrt) 
    std::cerr << "TestSuite did not regsiter uncaught exception!" << std::endl;
  if (succ) 
    std::cerr << "TestSuite did incorrectly registered failure for successful test!" << std::endl;
  
  int success = !fail + !abrt + !expt + succ;
  return success;
}
