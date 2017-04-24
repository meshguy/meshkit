# A from-scratch CMake build of lemon for meshkit. Lemon is a stand-alone
# project and does lots of things not necessary inside of meshkit.

set(lemon_srcs
  arg_parser.cc
  base.cc
  color.cc
  lp_base.cc
  lp_skeleton.cc
  random.cc
  lemon/bits/windows.cc)

set(lemon_headers
  lemon/adaptors.h
  lemon/arg_parser.h
  lemon/assert.h
  lemon/bellman_ford.h
  lemon/bfs.h
  lemon/bin_heap.h
  lemon/binomial_heap.h
  lemon/bucket_heap.h
  lemon/capacity_scaling.h
  lemon/cbc.h
  lemon/circulation.h
  lemon/clp.h
  lemon/color.h
  lemon/concept_check.h
  lemon/connectivity.h
  lemon/config.h
  lemon/core.h
  lemon/cost_scaling.h
  lemon/counter.h
  lemon/cplex.h
  lemon/cycle_canceling.h
  lemon/dfs.h
  lemon/dheap.h
  lemon/dijkstra.h
  lemon/dim2.h
  lemon/dimacs.h
  lemon/edge_set.h
  lemon/elevator.h
  lemon/error.h
  lemon/euler.h
  lemon/fib_heap.h
  lemon/fractional_matching.h
  lemon/full_graph.h
  lemon/glpk.h
  lemon/gomory_hu.h
  lemon/graph_to_eps.h
  lemon/grid_graph.h
  lemon/hartmann_orlin_mmc.h
  lemon/howard_mmc.h
  lemon/hypercube_graph.h
  lemon/karp_mmc.h
  lemon/kruskal.h
  lemon/hao_orlin.h
  lemon/lgf_reader.h
  lemon/lgf_writer.h
  lemon/list_graph.h
  lemon/lp.h
  lemon/lp_base.h
  lemon/lp_skeleton.h
  lemon/maps.h
  lemon/matching.h
  lemon/math.h
  lemon/min_cost_arborescence.h
  lemon/nauty_reader.h
  lemon/network_simplex.h
  lemon/pairing_heap.h
  lemon/path.h
  lemon/planarity.h
  lemon/preflow.h
  lemon/quad_heap.h
  lemon/radix_heap.h
  lemon/radix_sort.h
  lemon/random.h
  lemon/smart_graph.h
  lemon/soplex.h
  lemon/static_graph.h
  lemon/suurballe.h
  lemon/time_measure.h
  lemon/tolerance.h
  lemon/unionfind.h
  lemon/bits/windows.h
  lemon/bits/alteration_notifier.h
  lemon/bits/array_map.h
  lemon/bits/bezier.h
  lemon/bits/default_map.h
  lemon/bits/edge_set_extender.h
  lemon/bits/enable_if.h
  lemon/bits/graph_adaptor_extender.h
  lemon/bits/graph_extender.h
  lemon/bits/map_extender.h
  lemon/bits/path_dump.h
  lemon/bits/solver_bits.h
  lemon/bits/traits.h
  lemon/bits/variant.h
  lemon/bits/vector_map.h
  lemon/concepts/digraph.h
  lemon/concepts/graph.h
  lemon/concepts/graph_components.h
  lemon/concepts/heap.h
  lemon/concepts/maps.h
  lemon/concepts/path.h)

add_library(lemon
  ${lemon_srcs}
  ${lemon_headers})
include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR})

install(
  TARGETS   "lemon"
  EXPORT    meshkit
  RUNTIME   DESTINATION bin
  LIBRARY   DESTINATION lib
  ARCHIVE   DESTINATION lib
  COMPONENT Runtime)
