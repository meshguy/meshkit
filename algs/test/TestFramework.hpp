#ifndef TESTFRAMEWORK_HPP
#define TESTFRAMEWORK_HPP

#include <iostream>

#define ERROR(str) do {                           \
    std::cerr << (str) << std::endl;              \
    return 1;                                     \
  } while(false)

#define CHECK_TRUE(expr, str) do {                \
    if (!(expr)) {                                \
      ERROR(str);                                 \
    }                                             \
  } while(false)

#define CHECK_ERR(str) do {                       \
    if (err != iBase_SUCCESS) {                   \
      ERROR(str);                                 \
    }                                             \
  } while(false)

#define CHECK_THROW(expr) do {                    \
    try {                                         \
      expr;                                       \
    }                                             \
    catch(const std::exception &e) {              \
      ERROR(e.what());                            \
    }                                             \
  } while(false)

#define RUN_TEST(function) do {                   \
    std::cerr << "  " << #function << "... ";     \
    if(function() == 0) {                         \
      std::cerr << "passed" << std::endl;         \
    }                                             \
    else {                                        \
      std::cerr << "failed" << std::endl;         \
      result = 1;                                 \
    }                                             \
  } while(false)

#endif
