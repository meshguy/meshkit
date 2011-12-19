from MeshKit import *

file_name = "../data/partBed.smf"
outfile = "out.h5m"
print file_name 

mk = MKCore()

mk.load_mesh(file_name)

surfs = MEntVector()
mk.get_entities_by_dimension(2, surfs)

qm = mk.construct_meshop("MBGeomOp", surfs)

mk.setup_and_execute()

# round 2; after having the geometry model created, invoke camal mesher on
# mesh based geometry

mk.clear_graph()
# maybe this is not needed?
mk.delete_model_entities()

# moab model is still loaded, interpret it as geometry

# here it is a little oxymoron; pass a set we know it is empty
indx = mk.convert_db_to_mesh_based_geometry()

surfs.clear()
mk.get_entities_by_dimension(2, surfs)

print "surfs size:" , surfs.size()
mk.construct_meshop("CAMALPaver", surfs)

esize = SizingFunction(mk, -1, 100)

# we know only one surface
surfs[0].sizing_function_index(esize.core_index())

mk.setup_and_execute()

mk.remove_mesh_based_geometry(indx)

mk.save_mesh_from_model_ents("out.h5m", surfs)





