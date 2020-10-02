import globalvars
import numpy as np
from sopy_fem.initialization import initialization

def solve():
    n = globalvars.num_unknows
    neq = globalvars.neq

    #Incluiding the effect of the known displacements u_known in load vector bu
    bu = globalvars.asload[:n]
    Auk = globalvars.astiff[:n, n:]
    bu -= np.matmul(Auk, globalvars.u_known)
    num_known = neq - n

    #Solve system of equations => Auu·u_unknows = bu
    Auu = globalvars.astiff[:n,:n]
    u_unknows = np.linalg.solve(Auu, bu)

    globalvars.u_vec[:n] = u_unknows[:n]
    globalvars.u_vec[n:neq] = globalvars.u_known[:num_known]

    #Reactions calculation
    bk = globalvars.asload[n:]
    Aku = globalvars.astiff[n:, :n]
    Akk = globalvars.astiff[n:, n:]
    globalvars.react_vec = np.matmul(Aku, u_unknows) + np.matmul(Akk, globalvars.u_known) - bk
