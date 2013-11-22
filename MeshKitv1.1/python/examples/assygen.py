from MeshKit import *
import sys
mk = MKCore()
vols = MEntVector()
ag = mk.construct_meshop("AssyGen", vols)
TestDir = '../data'
ag.PrepareIO(sys.argv, TestDir)
mk.setup_and_execute()

