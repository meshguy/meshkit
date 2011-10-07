from MeshKit import *


file_name = "../data/pts.h5m"
outfile = "out.h5m"

mk = MKCore()

mk.load_mesh(file_name)

surfs = MEntVector()
mk.get_entities_by_dimension(2, surfs)

tm = mk.construct_meshop("TriangleMesher", surfs)
opts = 'pc'
direction = 3
fretting = 500.
tm.set_options(opts, direction, fretting)

#prepare qslim run
qopts = QslimOptions()
qopts.will_constrain_boundaries = 1
qopts.face_target = 4500
qopts.boundary_constraint_weight = 1000
qopts.height_fields = 1

qm = mk.construct_meshop("QslimMesher", surfs)
qm.set_options(qopts)

mk.add_arc(tm, qm )

mk.setup_and_execute()

# now, delete any model ents we had so far, and start over
mk.clear_graph()
mk.delete_model_entities()
# clear model ents
surfs.clear()

mk.save_mesh("AfterQslim.h5m")
# this will create a model ent of dimension 2, that will be geometrized
# it can be then cropped !!! it is the second part of load_mesh ()
mk.populate_model_ents(-1, 0, -1)

mk.get_entities_by_dimension(2, surfs)
print "before geometrization, surfs.size()=", surfs.size()
qm = mk.construct_meshop("MBGeomOp", surfs)

mk.setup_and_execute()

# so now, after another clear, we have what we want, a geometrized surface
# now, delete any model ents we had so far, and start over
mk.clear_graph()
mk.delete_model_entities()
# clear model ents
surfs.clear()

# start the split
mk.populate_model_ents(-1, 0, -1)
mk.get_entities_by_dimension(2, surfs)
qm = mk.construct_meshop("MBSplitOp", surfs)
print "start the split for surf with id " , surfs[0].id()
direction = [0.,  0.,  1. ]

# we know we have one surface () with gid 1
gid = surfs[0].id()
closed = 1
qm.set_options(gid, direction[0], direction[1], direction[2], closed)
qm.add_points(537200,  7680200,  -2000.)
qm.add_points(537800,  7680100,  -2000.)

qm.add_points(537700,  7680980,  -2000.)
qm.add_points(537400 , 7680930 , -2000.)
qm.add_points(537100 ,7680900 , -2000.)

mk.setup_and_execute()

# the moab model will stay in memory (splitted sets)
mk.clear_graph()
surfs.clear()
mk.delete_model_entities()

mk.initialize_mesh_based_geometry()

mk.get_entities_by_dimension(2, surfs)
# we should have 2 surfaces, one internal
print "surfs size:" , surfs.size()
surfs2 = MEntVector()

# use only the second surface
surfs2.push_back(surfs[1])
print " mesh surface with id: " , surfs2[0].id()
mk.construct_meshop("CAMALPaver", surfs2)

esize = SizingFunction(mk, -1, 80)

# we know only one surface
surfs2[0].sizing_function_index(esize.core_index())

mk.setup_and_execute()

mk.save_mesh_from_model_ents(outfile, surfs2)
