from MeshKit import *


file_name = "../data/TriangleInput.h5m"
outfile = "tout.h5m"

mk = MKCore()

mk.load_mesh(file_name)

surfs = MEntVector()
mk.get_entities_by_dimension(2, surfs)

qm = mk.construct_meshop("TriangleMesher", surfs)
opts = 'pc'
direction = 3
fretting = 3.
qm.set_options(opts, direction, fretting)

mk.setup_and_execute()

mk.save_mesh(outfile)
