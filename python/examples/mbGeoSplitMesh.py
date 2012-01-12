from MeshKit import *

file_name = "../data/partBed.smf"
outfile = "out.h5m"
print file_name 

mk = MKCore()

mk.load_mesh(file_name)

surfs = MEntVector()
mk.get_entities_by_dimension(2, surfs)

# no options for this op
mk.construct_meshop("MBGeomOp", surfs)

mk.setup_and_execute()

# round 2; after having the geometry model created, invoke camal mesher on
# mesh based geometry

mk.clear_graph()

# maybe this is not needed?
mk.delete_model_entities()

# clear model ents
surfs.clear()

# moab model is still loaded; populate now with model entities without geometry
# start the split
mk.populate_model_ents(-1, 0, -1)

mk.get_entities_by_dimension(2, surfs)

qm = mk.construct_meshop("MBSplitOp", surfs)
print "start the split for surf with id " , surfs[0].id()
direction = [0.,  0.,  1. ]

# we know we have one surface () with gid 1
gid = 1
closed = 1
min_dot = 0.8
qm.set_options(gid, direction[0], direction[1], direction[2], closed, min_dot)
qm.add_points(537200,  7680200,  -2000.)
qm.add_points(537800,  7680100,  -2000.)
qm.add_points(537700,  7680980,  -2000.)
qm.add_points(537400 , 7680930 , -2000.)
qm.add_points(537100 ,7680900 , -2000.)

mk.setup_and_execute()

# the moab model will stay in memory
mk.clear_graph()
surfs.clear()
mk.delete_model_entities()

# here, we should be done with splitting
indx = mk.initialize_mesh_based_geometry(0)

mk.get_entities_by_dimension(2, surfs)
# we should have 2 surfaces, one internal
print "surfs size:" , surfs.size()
surfs2 = MEntVector()

# use only the second surface
surfs2.push_back(surfs[1])
print " mesh surface with id: " , surfs2[0].id()
mk.construct_meshop("CAMALPaver", surfs2)

esize = SizingFunction(mk, -1, 100)

# we know only one surface
surfs2[0].sizing_function_index(esize.core_index())

mk.setup_and_execute()

mk.remove_mesh_based_geometry(indx)

mk.save_mesh_from_model_ents("out.h5m", surfs2)





