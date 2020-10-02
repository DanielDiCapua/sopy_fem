import numpy as np
import math
import sopy_fem.globalvars as globalvars
from sopy_fem.initialization import initialization
from sopy_fem.solids_utils import Bmat_TR03, Calc_Bmat
from sopy_fem.dmat import dmat_Solids2D, giveLocalStiffness
from sopy_fem.gauss_quadrature import Set_Ngauss, GaussQuadrature
from sopy_fem.solids_utils import derivCartesian, GiveNdof, GiveNnodes

def assembly():
    stiff_assembly()
    loads_assembly()

def stiff_assembly():
    mesh_type = globalvars.data["Mesh"]["ElemType"]
    for elem in globalvars.data["Mesh"]["Elements"]:
        if (mesh_type == "BAR02"):
            rigimat = stiffness_BAR02(elem)
        elif (mesh_type == "TRUSS02"):
            rigimat = stiffness_TRUSS02(elem)
        elif (mesh_type == "TR03"):
            rigimat = stiffness_TR03(elem)
        elif (mesh_type == "QU04"):
            rigimat = stiffness_QU04(elem)
        
        conec_list = elem["Connectivities"]
        doflist = np.zeros((len(conec_list)*globalvars.ndof), dtype=int)
        idof = 0
        for inode in range(len(conec_list)):
            id_node = conec_list[inode]-1
            for igl in range(globalvars.ndof):
                doflist[idof] =  globalvars.madgln[id_node, igl]
                idof += 1
        
        assamk(doflist, rigimat)


def loads_assembly():
    Loads = globalvars.data["Loads"]

    if ("Point_Loads" in Loads):
        doflist = np.zeros((globalvars.ndof), dtype=int)
        fvect = np.zeros((globalvars.ndof), dtype=float)
        for pointLoad in Loads["Point_Loads"]:
            id_node = pointLoad["Node"] - 1
            idof = 0
            for igl in range(globalvars.ndof):
                doflist[idof] = globalvars.madgln[id_node, igl]
                fvect[idof] = pointLoad["Values"][igl]
                idof += 1
            assamf(doflist, fvect)   

    if ("Surface_Loads" in Loads):
        print("There is not surface loads yet...")

        
def stiffness_BAR02(elem):
    mat_id = elem["MaterialId"] - 1
    material = globalvars.data["Materials"][mat_id]
    node1 = elem["Connectivities"][0] - 1
    x1 = globalvars.data["Mesh"]["Nodes"][node1]["x"]
    y1 = globalvars.data["Mesh"]["Nodes"][node1]["y"]
    node2 = elem["Connectivities"][1] - 1
    x2 = globalvars.data["Mesh"]["Nodes"][node2]["x"]
    y2 = globalvars.data["Mesh"]["Nodes"][node2]["y"]
    l = math.sqrt((x2 - x1)** 2 + (y2 - y1)** 2)
    k = giveLocalStiffness(material, l)
    rigimat = np.zeros((2, 2), dtype=float)
    rigimat[0, 0] = k
    rigimat[0, 1] = -k
    rigimat[1, 0] = -k
    rigimat[1, 1] = k
    return rigimat


def stiffness_TRUSS02(elem):
    mat_id = elem["MaterialId"] - 1
    material = globalvars.data["Materials"][mat_id]
    Young = material["Young"]
    Area = material["Area"]
    node1 = elem["Connectivities"][0] - 1
    x1 = globalvars.data["Mesh"]["Nodes"][node1]["x"]
    y1 = globalvars.data["Mesh"]["Nodes"][node1]["y"]
    node2 = elem["Connectivities"][1] - 1
    x2 = globalvars.data["Mesh"]["Nodes"][node2]["x"]
    y2 = globalvars.data["Mesh"]["Nodes"][node2]["y"]
    length = math.sqrt((x2 - x1)** 2 + (y2 - y1)** 2)
    cos_alpha = (x2 - x1) / length
    sin_alpha = (y2 - y1) / length
    
    k = Young * Area / length
    rigimat = np.zeros((4,4), dtype=float)
    rigimat[0, 0] = k * cos_alpha ** 2
    rigimat[0, 1] = k * cos_alpha * sin_alpha
    rigimat[0, 2] = -k * cos_alpha ** 2
    rigimat[0, 3] = -k * cos_alpha * sin_alpha
    rigimat[1, 0] = k*cos_alpha*sin_alpha
    rigimat[1, 1] = k * sin_alpha ** 2
    rigimat[1, 2] = -k * cos_alpha * sin_alpha
    rigimat[1, 3] = -k * sin_alpha ** 2
    rigimat[2, 0] = -k * cos_alpha ** 2
    rigimat[2, 1] = -k * cos_alpha * sin_alpha
    rigimat[2, 2] = k*cos_alpha** 2
    rigimat[2, 3] = k * cos_alpha * sin_alpha
    rigimat[3, 0] = -k * cos_alpha * sin_alpha
    rigimat[3, 1] = -k * sin_alpha ** 2
    rigimat[3, 2] = k * cos_alpha * sin_alpha
    rigimat[3, 3] = k * sin_alpha ** 2
    return rigimat

def stiffness_TR03(elem):
    mat_id = elem["MaterialId"] - 1
    material = globalvars.data["Materials"][mat_id]
    Dmat = dmat_Solids2D(material)

    Bmat, Area = Bmat_TR03(elem)
    
    DmatxB = np.matmul(Dmat, Bmat)
    thickness = 1.0
    if(material["Plane_Type"] == "Plane_Stress"):
        thickness =material["Thickness"]
    rigimat = np.matmul(Bmat.transpose(), DmatxB) * Area * thickness
    return rigimat


def stiffness_QU04(elem):
    ProblemType = globalvars.data["ProblemType"]
    mat_id = elem["MaterialId"] - 1
    material = globalvars.data["Materials"][mat_id]
    Dmat = dmat_Solids2D(material)
    thickness = 1.0
    if(material["Plane_Type"] == "Plane_Stress"):
        thickness = material["Thickness"]
    Dmat *= thickness
    return stiffnessMatCalc(elem, "QU04", ProblemType, Dmat)

def stiffnessMatCalc(elem, ElemType, ProblemType, Dmat):
    ngauss = Set_Ngauss(ElemType)
    Nodes = globalvars.data["Mesh"]["Nodes"]
    ndof = GiveNdof(ElemType, ProblemType)
    nnodes = GiveNnodes(ElemType)
    size = nnodes * ndof
    rigimat = np.zeros((size,size), dtype=float)
    pos_pg, W = GaussQuadrature(ElemType)

    for igauss in range(ngauss):
        det_Jac, derivCart = derivCartesian(elem, ElemType, Nodes, pos_pg[igauss])
        Bmat = Calc_Bmat(ElemType, derivCart, ProblemType)
        DxBmat = np.matmul(Dmat, Bmat)
        rigimat_pg = np.matmul(Bmat.transpose(), DxBmat) * det_Jac * W[igauss]
        rigimat += rigimat_pg
    return rigimat


def assamk(doflist, rigimat):
  ndof = len(doflist)
  for idofn in range(ndof):
      ipos = doflist[idofn]
      for jdofn in range(ndof):
        jpos = doflist[jdofn]
        globalvars.astiff[ipos, jpos] += rigimat[idofn, jdofn]


def assamf(doflist, fvect):
  ndof = len(doflist)
  for idofn in range(ndof):
      ipos = doflist[idofn]
      globalvars.asload[ipos] += fvect[idofn]
