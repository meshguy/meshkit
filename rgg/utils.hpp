#ifndef UTILS_HPP
#define UTILS_HPP

#include "CopyMesh.hpp"
#include "CopyGeom.hpp"
#include <iostream>

#define ERROR(a) {if (iBase_SUCCESS != err) std::cerr << a << std::endl;}
#define ERRORR(a,b) {if (iBase_SUCCESS != err) {std::cerr << a << std::endl; return b;}}

int get_copy_sets(CopyMesh *cm, iBase_EntitySetHandle orig_set,
                  const char **tag_names, const char **tag_vals, int num_tags);

int get_expand_sets(CopyMesh *cm, iBase_EntitySetHandle orig_set,
                    const char **tag_names, const char **tag_vals,
                    int num_tags);

int extend_expand_sets(CopyMesh *cm);

// Geometric versions

int get_copy_sets(CopyGeom *cg, iBase_EntitySetHandle orig_set,
                  const char **tag_names, const char **tag_vals, int num_tags);

int get_expand_sets(CopyGeom *cg, iBase_EntitySetHandle orig_set,
                    const char **tag_names, const char **tag_vals,
                    int num_tags);

int extend_expand_sets(CopyGeom *cg);

#endif
