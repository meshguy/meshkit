#ifndef MK_DEFINES_H
#define MK_DEFINES_H

#include "iBase.h"
#include "MBTypes.h"

#define ERROR(a) {if (iBase_SUCCESS != err) std::cerr << a << std::endl;}
#define ERRORR(a,b) {if (iBase_SUCCESS != err) {std::cerr << a << std::endl; return b;}}
#define ERRORRF(a) {if (iBase_SUCCESS != err) {std::cerr << a << std::endl; return false;}}
#define MBERRORR(a, b) {if (MB_SUCCESS != rval) {std::cerr << a << std::endl; return b;}}

const iBase_ErrorType iBase_ERROR_MAP[MB_FAILURE+1] = 
{
  iBase_SUCCESS, // MB_SUCCESS = 0,
  iBase_INVALID_ENTITY_HANDLE, // MB_INDEX_OUT_OF_RANGE,
  iBase_INVALID_ENTITY_TYPE, // MB_TYPE_OUT_OF_RANGE,
  iBase_MEMORY_ALLOCATION_FAILED, // MB_MEMORY_ALLOCATION_FAILED,
  iBase_INVALID_ENTITY_HANDLE, // MB_ENTITY_NOT_FOUND,
  iBase_NOT_SUPPORTED, // MB_MULTIPLE_ENTITIES_FOUND,
  iBase_TAG_NOT_FOUND, // MB_TAG_NOT_FOUND,
  iBase_FILE_NOT_FOUND, // MB_FILE_DOES_NOT_EXIST,
  iBase_FILE_WRITE_ERROR, // MB_FILE_WRITE_ERROR,
  iBase_NOT_SUPPORTED, // MB_NOT_IMPLEMENTED,
  iBase_TAG_ALREADY_EXISTS, // MB_ALREADY_ALLOCATED,
  iBase_FAILURE, // MB_VARIABLE_DATA_LENGTH,
  iBase_FAILURE, // MB_INVALID_SIZE,
  iBase_NOT_SUPPORTED, // MB_UNSUPPORTED_OPERATION,
  iBase_INVALID_ARGUMENT, // MB_UNHANDLED_OPTION
  iBase_FAILURE // MB_FAILURE};
};
#endif
