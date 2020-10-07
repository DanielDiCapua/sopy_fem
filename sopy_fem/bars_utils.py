import numpy as np
import math
import sopy_fem.globalvars as globalvars
from sopy_fem.initialization import initialization
from sopy_fem.dmat import giveLocalStiffness

def ElemTrusses_IntForces():
    globalvars.results["ElemIntForces"] = []
    mesh = globalvars.data["Mesh"]
    nodeList = mesh["Nodes"]
    elemList = mesh["Elements"]
    ElemIntForces = np.zeros(len(elemList), dtype=float)
    node = np.array([0, 0])
    idire = np.zeros((2), dtype=int)
    global_disp = np.zeros((2, 2), dtype=float)
    local_disp = np.zeros((2), dtype=float)
    for ielem, elem in enumerate(elemList):
        mat_id = elem["MaterialId"] - 1
        material = globalvars.data["Materials"][mat_id]
        node1 = elem["Connectivities"][0] - 1
        x1 = nodeList[node1]["x"]
        y1 = nodeList[node1]["y"]
        node2 = elem["Connectivities"][1] - 1
        x2 = nodeList[node2]["x"]
        y2 = nodeList[node2]["y"]
        length = math.sqrt((x2 - x1)** 2 + (y2 - y1)** 2)
        cos_alpha = (x2 - x1) / length
        sin_alpha = (y2 - y1) / length
        k = giveLocalStiffness(material, length)

        for inode in range(2):
            node[inode] = elem["Connectivities"][inode] - 1
            for igl in range(globalvars.ndof):
                idire = globalvars.madgln[node[inode], igl]
                global_disp[inode, igl] = globalvars.u_vec[idire]
            local_disp[inode]=global_disp[inode,0]*cos_alpha + global_disp[inode,1]*sin_alpha
        
        axialForce = k * (local_disp[1] - local_disp[0])
        ElemIntForces[ielem] = axialForce
        axialForcesFormatted = "{:10.4e}".format(axialForce)
        globalvars.results["ElemIntForces"].append(axialForcesFormatted)

    return ElemIntForces

def ElemBars_IntFluxes():
    globalvars.results["ElemIntFluxes"] = []
    mesh = globalvars.data["Mesh"]
    nodeList = mesh["Nodes"]
    elemList = mesh["Elements"]
    ElemIntFluxes = np.zeros(len(elemList), dtype=float)
    node = np.array([0, 0])
    idire = np.zeros((2), dtype=int)
    var_vec = np.zeros((2), dtype=float)
    for ielem, elem in enumerate(elemList):
        mat_id = elem["MaterialId"] - 1
        material = globalvars.data["Materials"][mat_id]
        node1 = elem["Connectivities"][0] - 1
        x1 = nodeList[node1]["x"]
        y1 = nodeList[node1]["y"]
        node2 = elem["Connectivities"][1] - 1
        x2 = nodeList[node2]["x"]
        y2 = nodeList[node2]["y"]
        length = math.sqrt((x2 - x1)** 2 + (y2 - y1)** 2)
        k = giveLocalStiffness(material, length)

        for inode in range(2):
            node[inode] = elem["Connectivities"][inode] - 1
            idire = globalvars.madgln[node[inode], 0]
            var_vec[inode] = globalvars.u_vec[idire]
        
        if (globalvars.data["ProblemType"] == "Thermal"):
            k *= (-1.0)
        flux = k * (var_vec[1] - var_vec[0])
        ElemIntFluxes[ielem] = flux
        fluxFormatted = "{:10.4e}".format(flux)
        globalvars.results["ElemIntFluxes"].append(fluxFormatted)

    return ElemIntFluxes



