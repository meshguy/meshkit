#include "meshkit/Error.hpp"
#include <stdio.h>
#include <stdarg.h>

namespace MeshKit {

Error::Error( int err, const char* format, ... )
{
  char buffer[512];
  va_list args;
  va_start(args, format);
  vsnprintf(buffer, sizeof(buffer)-1, format, args);
  errDescription = buffer;
  va_end(args);
}

}
