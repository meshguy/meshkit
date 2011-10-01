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

mk.save_mesh(outfile)




