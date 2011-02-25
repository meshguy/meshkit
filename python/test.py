from MeshKit import *

mk = MKCore()
for i in range(mk.num_meshops()):
    op = mk.meshop_proxy(i)
    print op.name()
