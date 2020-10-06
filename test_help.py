import sys
from sopy_fem.sopy_fem_help import sopy_fem_help

def test_help():
    args = sys.argv[1:]
    exampleType = args[0]
    sopy_fem_help(exampleType)

test_help()