/**
 * Copyright 2006 Sandia Corporation.  Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Coroporation, the U.S. Government
 * retains certain rights in this software.
 * 
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * 
 */
#ifndef CAMAL_INTERFACE_HPP
#define CAMAL_INTERFACE_HPP

#include "camel.hpp"

bool CAMAL_mesh_entity(CMEL *cmel,
                       iBase_EntityHandle gentity, double mesh_size,
                       int mesh_intervals, const bool force_intervals,
                       std::vector<iBase_EntityHandle> &new_entities,
                       const bool quadMesh = false);

bool CAMAL_mesh_trimmed_surface(CMEL * cmel, iBase_EntityHandle surface,
      double mesh_size, std::vector<iBase_EntityHandle> &new_entities,
      std::vector<double> trimmingBoundary, const bool quadMesh = false);

#endif
