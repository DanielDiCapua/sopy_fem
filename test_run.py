from os.path import join
import sys
from sopy_fem.sopy_fem_run import sopy_fem_run

def test_run():
    args = sys.argv[1:]
    
    if len(args) == 0:
        exampleType = input("Introduce an example: ")
    else:
        exampleType = args[0]
        
    dataFileName = join("./sopy_fem/Examples/", exampleType, "data.json")
    
    sopy_fem_run(dataFileName)

test_run()

'''
Choose a value for "exampleType" from among the following options to run a example:
  -  structural_FRAME02
  -  structural_TRUSS02
  -  mechanics_BR02
  -  dynamics_FRAME02
  -  thermal_BR02
  -  mechanics_QU04
  -  electrical_BR02
  -  mechanics_TR03
  -  dynamics_TRUSS02
'''
