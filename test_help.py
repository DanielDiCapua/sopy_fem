import sys
from sopy_fem.sopy_fem_help import sopy_fem_help

def test_help():
    args = sys.argv[1:]
    exampleType = ""
    basePath = "sopy_fem/Examples/"
    if (len(args)  != 0):
        exampleType = args[0]
    sopy_fem_help(exampleType, basePath)

test_help()

'''
Choose a value for "exampleType" from among the following options to print the json file into the terminal:
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
