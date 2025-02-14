import numpy as np
import scipy.linalg
import math
import sopy_fem.globalvars as globalvars
from sopy_fem.initialization import initialization
from sopy_fem.assembly import dynamics_loads_assembly

def solve():
    n = globalvars.num_unknows
    neq = globalvars.neq

    if(globalvars.data["AnalysisType"] == "StaticAnalysis"):
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
    else:
        Auu = globalvars.astiff[:n,:n]
        Muu = globalvars.amassmat[:n,:n]
        eigvals, eigvecs = scipy.linalg.eig(Auu, Muu)
        eigvals_ord, eigvecs_ord = orderingEigvalsAndEigevecs(eigvals, eigvecs)
        numModes = globalvars.data["Dynamic_Analysis_Description"]["Num_Modes"]
        for imode in range(numModes):
            globalvars.natFreqVec[imode] = math.sqrt(eigvals_ord[imode])/(2.0*math.pi)
            globalvars.vibrationModes[imode,:n] =  eigvecs_ord[:neq, imode]
            A = globalvars.vibrationModes[imode,:]
            MxA = np.matmul(globalvars.amassmat, globalvars.vibrationModes[imode,:])
            AtxMxA = np.dot(A, MxA)
            globalvars.vibrationModes[imode,:neq] = A/math.sqrt(AtxMxA)
        DynamicsAnalyis()


def DynamicsAnalyis():
    numModes = globalvars.data["Dynamic_Analysis_Description"]["Num_Modes"]
    numInc = globalvars.data["Dynamic_Analysis_Description"]["Num_increments"]
    deltaT = globalvars.data["Dynamic_Analysis_Description"]["DeltaT"]
    nu = globalvars.data["Dynamic_Analysis_Description"]["Damping_ratio"]
    neq = globalvars.neq
    gamma = 0.5
    beta = 0.25
    for imode in range(numModes):
        omega = 2.0 * math.pi * globalvars.natFreqVec[imode]
        A = globalvars.vibrationModes[imode,:neq]
        m = 1.0
        c = 2.0 * nu * omega
        k = omega ** 2.0
        a_old = 0.0
        v_old = 0.0
        u_old = 0.0
        kmod = k + m / (beta * deltaT ** 2.0) + c * gamma / (beta * deltaT)
        for istep in range(numInc):
            t = istep * deltaT
            dynamics_loads_assembly(t)
            pt = np.dot(A, globalvars.asload)
            pt_mod = pt + m * (u_old / (beta * deltaT ** 2.0) + v_old / (beta * deltaT) + a_old * (1.0 / (2.0 * beta) - 1.0))
            pt_mod += c * (u_old * gamma / (beta * deltaT) + v_old * ((gamma / beta) - 1.0) + a_old * deltaT * ((gamma / (2.0 * beta)) - 1.0))
            u_new = pt_mod / kmod
            v_new = (gamma / (beta * deltaT)) * (u_new - u_old) + (1.0-(gamma / beta)) * v_old + (1.0-(gamma / (2.0 * beta))) * deltaT * a_old
            a_new = (1.0 / (beta * deltaT ** 2.0)) * (u_new - u_old - v_old * deltaT) - ((1.0 / (2.0 * beta)) - 1.0) * a_old
            u_old = u_new
            v_old = v_new
            a_old = a_new
            globalvars.modal_disp[imode, istep] = u_new
            
    for istep in range(numInc):
        for imode in range(numModes):
            A = globalvars.vibrationModes[imode,:neq]
            globalvars.dynamics_uvec[istep,:] += globalvars.modal_disp[imode, istep] * A


def orderingEigvalsAndEigevecs(eigvals, eigvecs):
    num_rows, num_columns = eigvecs.shape
    eigvals_ord = np.zeros((num_rows), dtype='float')
    eigvecs_ord = np.zeros((num_rows, num_columns), dtype='float')
    num_columns += 1
    big_mat = np.zeros((num_rows, num_columns), dtype='float')
    for irow in range(num_rows):
        big_mat[irow, 0] = eigvals[irow].real
        for icol in range(1,num_columns):
            big_mat[irow, icol] = eigvecs[(icol - 1), irow]
    
    big_mat_sorted = big_mat[big_mat[:, 0].argsort()]

    for irow in range(num_rows):
        eigvals_ord[irow] = big_mat_sorted[irow, 0]
        for icol in range(1,num_columns):
            eigvecs_ord[(icol - 1), irow] = big_mat_sorted[irow, icol]
    
    return eigvals_ord, eigvecs_ord
    

    
    

            
            
            

        

