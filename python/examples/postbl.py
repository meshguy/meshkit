from MeshKit import *
import sys
mk = MKCore()
vols = MEntVector()
pb = mk.construct_meshop("PostBL", vols)
TestDir = '../data'
pb.PrepareIO(sys.argv, TestDir)
mk.setup_and_execute()

