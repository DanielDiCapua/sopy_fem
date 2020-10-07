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