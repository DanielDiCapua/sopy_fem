import sopy_fem.globalvars as globalvars
from sopy_fem.read_data import read_data
from sopy_fem.initialization import initialization
from sopy_fem.assembly import assembly
from sopy_fem.solver import solve
from sopy_fem.postprocess import postprocess


def main():
    read_data()
    initialization()
    assembly()
    solve()
    postprocess()


if __name__ == '__main__':
    main()
    