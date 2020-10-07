import sys
from sopy_fem.sopy_fem_run import sopy_fem_run

def test_run():
    args = sys.argv[1:]
    dataFileName = ""
    if (len(args)  != 0):
        dataFileName = args[0]
    sopy_fem_run(dataFileName)

test_run()
