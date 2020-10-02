import globalvars
from read_data import read_data
from initialization import initialization
from assembly import assembly
from solver import solve
from postprocess import postprocess


def main():
    read_data()
    initialization()
    assembly()
    solve()
    postprocess()


if __name__ == '__main__':
    main()
    