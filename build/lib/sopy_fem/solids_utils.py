import numpy as np
import math
import sopy_fem.globalvars as globalvars
from sopy_fem.initialization import initialization
from sopy_fem.gauss_quadrature import GaussQuadrature, Set_Ngauss
from sopy_fem.dmat import dmat_Solids2D

def GiveNdime(ElemType):
    if(ElemType == "BAR02" or ElemType == "BAR03" or ElemType == "TRUSS02"):
        ndime=1
    elif(ElemType == "TR03" or ElemType == "TR06" or ElemType == "QU04" or ElemType == "QU08" or ElemType == "QU09"):
      ndime=2
    return ndime

def GiveNComp(ElemType, ProblemType):
  if(ProblemType == "Structural_Mechanics"):
    if(ElemType == "BAR02" or ElemType == "BAR03" or ElemType == "TRUSS02"):
        ncomp=1
    elif(ElemType == "TR03" or ElemType == "TR06" or ElemType == "QU04" or ElemType == "QU08" or ElemType == "QU09"):
        ncomp = 3
  elif (ProblemType == "Thermal"):
    if(ElemType == "BAR02" or ElemType == "BAR03" or ElemType == "TRUSS02"):
      ncomp=1
    elif(ElemType == "TR03" or ElemType == "TR06" or ElemType == "QU04" or ElemType == "QU08" or ElemType == "QU09"):
      ncomp = 2
  return ncomp

def GiveNdof(ElemType, ProblemType):
  if(ProblemType == "Structural_Mechanics"):
    if(ElemType == "BAR02" or ElemType == "BAR03"):
      ndof=1
    elif(ElemType == "TRUSS02" or ElemType == "TR03" or ElemType == "TR06" or ElemType == "QU04" or ElemType == "QU08" or ElemType == "QU09"):
      ndof=2
  elif (ProblemType == "Thermal"):
    if(ElemType == "BAR02" or ElemType == "BAR03"):
      ndof=1
    elif(ElemType == "TR03" or ElemType == "TR06" or ElemType == "QU04" or ElemType == "QU08" or ElemType == "QU09"):
      ndof=2
  return ndof

def GiveNnodes(ElemType):
  if(ElemType == "BAR02" or ElemType == "TRUSS02"):
    nnodes = 2
  elif(ElemType == "BAR03"):
    nnodes = 3
  elif(ElemType == "TR03"):
    nnodes = 3
  elif(ElemType == "TR06"):
    nnodes = 6
  elif(ElemType == "QU04"):
    nnodes = 4
  elif(ElemType == "QU08"):
    nnodes = 8
  elif(ElemType == "QU09"):
    nnodes = 9
  return nnodes

def FuncForm(ElemType, pos_pg, igauss):
  ndime = GiveNdime(ElemType)
  if(ndime == 1):
    s = pos_pg[igauss, 0]
    t = 0
  elif(ndime == 2):
    s = pos_pg[igauss, 0]
    t = pos_pg[igauss, 1]

  nnodes = GiveNnodes(ElemType)
  fform = np.zeros((nnodes), dtype=float)
  if(ElemType == "BAR02" or ElemType == "TRUSS02"):
    fform[0] = 0.5 * (1.0 - s)
    fform[1] = 0.5 * (1.0 + s)
  elif(ElemType == "BAR03"):
    fform[0] = 0.5 * s * (s - 1)
    fform[1] = (1 + s) * (1 - s)
    fform[2] = 0.5 * s * (1 + s)
  elif(ElemType == "TR03"):
    fform[0] = 1 - s - t
    fform[1] = s
    fform[2] = t
  elif(ElemType == "TR06"):
    st = s * t
    tt = t * t
    ss = s * s
    fform[0] = 1 - (s + t) * 3 + (ss + tt) * 2 + st * 4
    fform[1] = (s - ss - st) * 4
    fform[2] = ss * 2 - s
    fform[3] = 4 * st
    fform[4] = tt * 2 - t
    fform[5] = (t - tt - st) * 4
  elif(ElemType == "QU04"):
    st = s * t
    fform[0] = (1 - s - t + st) / 4
    fform[1] = (1 + s - t - st) / 4
    fform[2] = (1 + s + t + st) / 4
    fform[3] = (1 - s + t - st) / 4
  elif(ElemType == "QU08"):
    st = s * t
    tt = t * t
    ss = s * s
    sst = ss * t
    stt = tt * s
    fform[0] = (-1 + st + ss + tt - sst - stt) * 0.25
    fform[1] = (1 - t - ss + sst) * 0.5
    fform[2] = (-1 - st + ss + tt - sst + stt) * 0.25
    fform[3] = (1 + s - tt - stt) * 0.5
    fform[4] = (-1 + st + ss + tt + sst + stt) * 0.25
    fform[5] = (1 + t - ss - sst) * 0.5
    fform[6] = (-1 - st + ss + tt + sst - stt) * 0.25
    fform[7] = (1 - s - tt + stt) * 0.5
  elif(ElemType == "QU09"):
    st = s * t
    tt = t * t
    ss = s * s
    s9 = s - 1
    t9 = t - 1
    s1 = s + 1
    t1 = t + 1
    fform[0] = 0.25 * s9 * st * t9
    fform[1] = 0.5 * (1 - ss) * t * t9
    fform[2] = 0.25 * s1 * st * t9
    fform[3] = 0.5 * s * s1 * (1 - tt)
    fform[4] = 0.25 * s1 * st * t1
    fform[5] = 0.5 * (1 - ss) * t * t1
    fform[6] = 0.25 * s9 * st * t1
    fform[7] = 0.5 * s * s9 * (1 - tt)
    fform[8] = (1 - ss) * (1 - tt)
  return fform

def DerivNatCoord(ElemType, pos_pg):
  ndime = GiveNdime(ElemType)
  if(ndime == 1):
    s = pos_pg[0]
    t = 0
  elif(ndime == 2):
    s = pos_pg[0]
    t = pos_pg[1]

  nnodes = GiveNnodes(ElemType)
  deriv = np.zeros((ndime,nnodes), dtype=float)

  #deriv[ider,inode]
  if(ElemType == "BAR02" or ElemType == "TRUSS02"):
    deriv[0,0] =-0.5
    deriv[0, 1] = 0.5
  elif(ElemType == "BAR03"):
    deriv[0,0] = s - 0.5
    deriv[0,1] = -2 * s
    deriv[0, 2] = s - 0.5
  elif(ElemType == "TR03"):
    deriv[0,0] = -1
    deriv[0,1] =  1
    deriv[0,2] =  0
    deriv[1,0] = -1
    deriv[1,1] =  0
    deriv[1, 2] = 1
  elif(ElemType == "TR06"):
    deriv[0,0]=-3+(s+t)*4
    deriv[0,1]=(1-s*2-t)*4
    deriv[0,2]=s*4-1
    deriv[0,3]=t*4
    deriv[0,4]=0
    deriv[0,5]=-t*4
    deriv[1,0]=-3+(t+s)*4
    deriv[1,1]=-s*4
    deriv[1,2]=0
    deriv[1,3]=4*s
    deriv[1,4]=t*4-1
    deriv[1,5]=(1-t*2-s)*4
  elif(ElemType == "QU04"):
    deriv[0,0] = (-1 + t) / 4
    deriv[0,1] = (1 - t) / 4
    deriv[0,2] = (1 + t) / 4
    deriv[0,3] = (-1 - t) / 4
    deriv[1,0] = (-1 + s) / 4
    deriv[1,1] = (-1 - s) / 4
    deriv[1,2] = (1 + s) / 4
    deriv[1,3] = (1 - s) / 4
  elif(ElemType == "QU08"):
    s2 = 2 * s
    t2 = 2 * t
    st = s * t
    st2 = s * t * 2
    tt = t * t
    ss = s * s
    deriv[0,0] = (t + s2 - st2 - tt) * 0.25
    deriv[0,1] = -s + st
    deriv[0,2] = (-t + s2 - st2 + tt) * 0.25
    deriv[0,3] = (1.0 - tt) * 0.5
    deriv[0,4] = (t + s2 + st2 + tt) * 0.25
    deriv[0,5] = -s - st
    deriv[0,6] = (-t + s2 + st2 - tt) * 0.25
    deriv[0,7] = (-1.0 + tt) * 0.5
    deriv[1,0] = (s + t2 - ss - st2) * 0.25
    deriv[1,1] = (-1.0 + ss) * 0.5
    deriv[1,2] = (-s + t2 - ss + st2) * 0.25
    deriv[1,3] = -t - st
    deriv[1,4] = (s + t2 + ss + st2) * 0.25
    deriv[1,5] = (1.0 - ss) * 0.5
    deriv[1,6] = (-s + t2 + ss - st2) * 0.25
    deriv[1,7] = -t + st
  elif(ElemType == "QU09"):
    s2 = 2 * s
    t2 = 2 * t
    st = s * t
    st2 = s * t * 2
    tt = t * t
    ss = s * s
    s9 = s - 1
    t9 = t - 1
    s1 = s + 1
    t1 = t + 1
    deriv[0,0] = 0.25 * t * t9 * (-1 + s2)
    deriv[0,1] = -st * t9
    deriv[0,2] = 0.25 * (1 + s2) * t * t9
    deriv[0,3] = 0.5 * (1 + s2) * (1 - tt)
    deriv[0,4] = 0.25 * (1 + s2) * t * t1
    deriv[0,5] = -st * t1
    deriv[0,6] = 0.25 * (-1 + s2) * t * t1
    deriv[0,7] = 0.5 * (-1 + s2) * (1 - tt)
    deriv[0,8] = -s2 * (1 - tt)
    deriv[1,0] = 0.25 * (-1 + t2) * s * s9
    deriv[1,1] = 0.5 * (1 - ss) * (-1 + t2)
    deriv[1,2] = 0.25 * s * s1 * (-1 + t2)
    deriv[1,3] = -st * s1
    deriv[1,4] = 0.25 * s * s1 * (1 + t2)
    deriv[1,5] = 0.5 * (1 - ss) * (1 + t2)
    deriv[1,6] = 0.25 * s * s9 * (1 + t2)
    deriv[1,7] = -st * s9
    deriv[1,8] = -t2 * (1 - ss)
  return deriv

def det_Jacob(elem, ElemType, Nodes, pos_pg): 
  ndime = GiveNdime(ElemType)
  nnodes = elem.give_Nnodes()
  corel = np.zeros((nnodes,ndime), dtype=float)
  for inode in range(nnodes):
    id_node = elem["Connectivities"][inode] - 1
    for idime in range(ndime):
      if(idime == 0):
        corel[inode, idime] = Nodes[id_node]["x"]
      elif (idime == 1):
        corel[inode, idime] = Nodes[id_node]["y"]

  deriv_nat = DerivNatCoord(ElemType, pos_pg)
  jacobMat = np.matmul(deriv_nat,corel)
  det_jacobMat = np.linalg.det(jacobMat)
  if (det_jacobMat==0):
    raise Exception("Error: Jacobian negative")
  return det_jacobMat

def derivCartesian(elem, ElemType, Nodes, pos_pg):
  ndime = GiveNdime(ElemType)
  nnodes = GiveNnodes(ElemType)
  corel = np.zeros((nnodes, ndime), dtype=float)
  if (ElemType == "TRUSS02"):
    node1 = elem["Connectivities"][0] - 1
    x1 = globalvars.data["Mesh"]["Nodes"][node1]["x"]
    y1 = globalvars.data["Mesh"]["Nodes"][node1]["y"]
    node2 = elem["Connectivities"][1] - 1
    x2 = globalvars.data["Mesh"]["Nodes"][node2]["x"]
    y2 = globalvars.data["Mesh"]["Nodes"][node2]["y"]
    length = math.sqrt((x2 - x1)** 2 + (y2 - y1)** 2)

  for inode in range(nnodes):
    id_node = elem["Connectivities"][inode] - 1
    for idime in range(ndime):
      if (idime == 0):
        if (ElemType == "TRUSS02"):
          if (inode == 0):
            corel[inode, idime] = 0.0
          else:
            corel[inode, idime] = length
        else:
          corel[inode, idime] = Nodes[id_node]["x"]
      elif (idime == 1):
        corel[inode, idime] = Nodes[id_node]["y"]

  deriv_nat = DerivNatCoord(ElemType, pos_pg)
  jacobMat = np.matmul(deriv_nat, corel)
  det_jacobMat = np.linalg.det(jacobMat)
  if (det_jacobMat<=0):
    raise Exception("Error: Jacobian negative")
  JacobMat_inv = np.linalg.inv(jacobMat)
  derivCart = np.matmul(JacobMat_inv, deriv_nat)
  return det_jacobMat, derivCart
 

def Smooth_Strains():
    num_nodes = len(globalvars.data["Mesh"]["Nodes"])
    nodal_Eps = np.zeros((num_nodes, 3), dtype=float)
    node_counter = np.zeros((num_nodes), dtype=int)
    ProblemType = globalvars.data["ProblemType"]
    ElemType = globalvars.data["Mesh"]["ElemType"]
    nnode = GiveNnodes(ElemType)
    ngauss = Set_Ngauss(ElemType)
    ncomp = GiveNComp(ElemType, ProblemType)
    Eps_vec = np.zeros((ncomp, ngauss), dtype=float)
    num_elems = len(globalvars.data["Mesh"]["Elements"])
    for icomp in range(ncomp):
      globalvars.Eps_elem_vec[icomp] = np.zeros((num_elems, ngauss), dtype=float)
    ielem = 0
    for elem in globalvars.data["Mesh"]["Elements"]:
        if (ElemType == "TR03"):
          Esp_vec_TR03 = Strain_TR03(elem)
          for igauss in range(ngauss):
            Eps_vec[:ncomp,igauss] = Esp_vec_TR03[:ncomp]
        else:
          for igauss in range(ngauss):
            Eps_vec[:, igauss] = Strain_Solid(elem, ElemType, ProblemType, igauss)

        for igauss in range(ngauss):
          for icomp in range(ncomp):
            globalvars.Eps_elem_vec[icomp][ielem,igauss] = Eps_vec[icomp, igauss]
        ielem +=1
        
        for inode in range(nnode):
            id_node = elem["Connectivities"][inode] - 1
            node_counter[id_node] += 1
            nodal_Eps[id_node,:ncomp] += Eps_vec[:ncomp, inode]
    
    for inode in range(num_nodes):
        nodal_Eps[inode,:ncomp] /= node_counter[inode]

    return nodal_Eps

def Smooth_Stresses():
    num_nodes = len(globalvars.data["Mesh"]["Nodes"])
    nodal_Sigma = np.zeros((num_nodes, 3), dtype=float)
    node_counter = np.zeros((num_nodes), dtype=int)
    ProblemType = globalvars.data["ProblemType"]
    ElemType = globalvars.data["Mesh"]["ElemType"]
    nnode = GiveNnodes(ElemType)
    ngauss = Set_Ngauss(ElemType)
    ncomp = GiveNComp(ElemType, ProblemType)
    Eps_vec = np.zeros((ncomp, ngauss), dtype=float)
    Sigma_vec = np.zeros((ncomp, ngauss), dtype=float)
    num_elems = len(globalvars.data["Mesh"]["Elements"])
    for icomp in range(ncomp):
      globalvars.Sigma_elem_vec[icomp] = np.zeros((num_elems, ngauss), dtype=float)
    ielem = 0
    for elem in globalvars.data["Mesh"]["Elements"]:
        if (ElemType == "TR03"):
          Sigma_vec_TR03 = Stress_TR03(elem)
          for igauss in range(ngauss):
            Sigma_vec[:ncomp,igauss] = Sigma_vec_TR03[:ncomp]
        else:
          for igauss in range(ngauss):
            if(len(globalvars.Eps_elem_vec[0]) == 0):
              Eps_vec[:, igauss] = Strain_Solid(elem, ElemType, ProblemType, igauss)
            else:
              for icomp in range(ncomp):
                Eps_vec[icomp, igauss] = globalvars.Eps_elem_vec[icomp][ielem,igauss]
            Sigma_vec[:, igauss] = Stress_Solid(elem, ElemType, Eps_vec[:, igauss])

        for igauss in range(ngauss):
          for icomp in range(ncomp):
            globalvars.Sigma_elem_vec[icomp][ielem, igauss] = Sigma_vec[icomp, igauss]
        ielem += 1
        
        for inode in range(nnode):
            id_node = elem["Connectivities"][inode] - 1
            node_counter[id_node] += 1
            nodal_Sigma[id_node,:ncomp] += Sigma_vec[:ncomp, inode]
    
    for inode in range(num_nodes):
        nodal_Sigma[inode,:ncomp] /= node_counter[inode]

    return nodal_Sigma


def Smooth_Forces():
    num_nodes = len(globalvars.data["Mesh"]["Nodes"])
    nodal_Forces = np.zeros((num_nodes, 3), dtype=float)
    node_counter = np.zeros((num_nodes), dtype=int)
    ProblemType = globalvars.data["ProblemType"]
    ElemType = globalvars.data["Mesh"]["ElemType"]
    nnode = GiveNnodes(ElemType)
    ngauss = Set_Ngauss(ElemType)
    ncomp = GiveNComp(ElemType, ProblemType)
    Eps_vec = np.zeros((ncomp, ngauss), dtype=float)
    Sigma_vec = np.zeros((ncomp, ngauss), dtype=float)
    Forces_vec = np.zeros((ncomp, ngauss), dtype=float)
    ielem = 0
    for elem in globalvars.data["Mesh"]["Elements"]:
        mat_id = elem["MaterialId"] - 1
        material = globalvars.data["Materials"][mat_id]
        factor = 1.0
        if(ElemType == "TR03" or ElemType == "TR06" or ElemType == "QU04" or ElemType == "QU08" or ElemType == "QU09"):
          if (material["Plane_Type"] == "Plane_Stress"):
            factor = material["Thickness"]
        elif (ElemType == "TRUSS02" or ElemType == "BAR02" or ElemType == "BAR03"):
          factor = material["Area"]
        if (ElemType == "TR03"):
          Sigma_vec_TR03 = Stress_TR03(elem)
          Forces_vec[:ncomp,:] = Sigma_vec_TR03[:ncomp]*factor
        else:
          for igauss in range(ngauss):
            if (len(globalvars.Sigma_elem_vec[0]) == 0):
              Eps_vec[:, igauss] = Strain_Solid(elem, ElemType, ProblemType, igauss)
              Sigma_vec[:, igauss] = Stress_Solid(elem, ElemType, Eps_vec[:, igauss])
            else:
              for icomp in range(ncomp):
                Sigma_vec[icomp, igauss] = globalvars.Sigma_elem_vec[icomp][ielem, igauss]
            Forces_vec[:, igauss] = Sigma_vec[:, igauss] * factor
        ielem += 1
        
        for inode in range(nnode):
            id_node = elem["Connectivities"][inode] - 1
            node_counter[id_node] += 1
            nodal_Forces[id_node,:ncomp] += Forces_vec[:ncomp, inode]
    
    for inode in range(num_nodes):
        nodal_Forces[inode,:ncomp] /= node_counter[inode]

    return nodal_Forces   

def Bmat_TR03(elem):
    x = np.zeros((3), dtype=float)
    y = np.zeros((3), dtype=float)
    b = np.zeros((3), dtype=float)
    c = np.zeros((3), dtype=float)
    Bmat = np.zeros((3,6), dtype=float)
    for inode in range(3):
        id_node = elem["Connectivities"][inode] - 1
        x[inode] = globalvars.data["Mesh"]["Nodes"][id_node]["x"]
        y[inode] = globalvars.data["Mesh"]["Nodes"][id_node]["y"]

    b[0] = y[1] - y[2]
    b[1] = y[2] - y[0]
    b[2] = y[0] - y[1]
    c[0] = x[2] - x[1]
    c[1] = x[0] - x[2]
    c[2] = x[1] - x[0]
    
    coords_Mat = np.array([[1,x[0],y[0]], [1,x[1],y[1]],[1,x[2],y[2]]]) 
    Area = 0.5 * np.linalg.det(coords_Mat)
    for inode in range(3):
        icol1 = inode * 2
        icol2 = inode * 2 + 1
        Bmat[0, icol1] = b[inode] / (2.0 * Area)
        Bmat[1, icol2] = c[inode] / (2.0 * Area)
        Bmat[2, icol1] = c[inode] / (2.0 * Area)
        Bmat[2, icol2] = b[inode] / (2.0 * Area)    

    return Bmat, Area

def Strain_TR03(elem):
    u_vec = np.zeros((6), dtype=float)
    for inode in range(3):
        id_node = elem["Connectivities"][inode] - 1
        for igl in range(2):
            idire = globalvars.madgln[id_node, igl]
            ipos = inode*2 + igl
            u_vec[ipos] = globalvars.u_vec[idire]
    
    Bmat, _ = Bmat_TR03(elem)
    Eps_vec = np.matmul(Bmat, u_vec)
    return Eps_vec

def Strain_Solid(elem, ElemType, ProblemType, igauss):
    nnodes = GiveNnodes(ElemType)
    ndof = globalvars.ndof
    uvec_size = nnodes * ndof
    u_vec = np.zeros((uvec_size), dtype=float)
    Nodes = globalvars.data["Mesh"]["Nodes"]
    for inode in range(nnodes):
        id_node = elem["Connectivities"][inode] - 1
        for igl in range(ndof):
            idire = globalvars.madgln[id_node, igl]
            ipos = inode*ndof + igl
            u_vec[ipos] = globalvars.u_vec[idire]
    
    if(ElemType == "TR03"):
      Bmat, _ = Bmat_TR03(elem)
    else:
      pos_pg, _ = GaussQuadrature(ElemType)
      _, derivCart = derivCartesian(elem, ElemType, Nodes, pos_pg[igauss])
      Bmat = Calc_Bmat(ElemType, derivCart, ProblemType)
  
    Eps_vec = np.matmul(Bmat, u_vec)
    return Eps_vec

def Stress_TR03(elem):
    mat_id = elem["MaterialId"] - 1
    material = globalvars.data["Materials"][mat_id]
    Dmat = dmat_Solids2D(material)
    Eps_vec = Strain_TR03(elem)
    Sigma_vec = np.matmul(Dmat, Eps_vec)
    return Sigma_vec

def Stress_Solid(elem, ElemType, Eps_vec):
    mat_id = elem["MaterialId"] - 1
    material = globalvars.data["Materials"][mat_id]
    if(ElemType == "BAR02" or ElemType == "BAR03" or ElemType == "TRUSS02"):
      Young = material["Young"]
      Sigma_vec = Young * Eps_vec
    else:
      Dmat = dmat_Solids2D(material)
      Sigma_vec = np.matmul(Dmat, Eps_vec)
    return Sigma_vec

def Calc_Bmat(ElemType, derivCart, ProblemType):
  ndime = GiveNdime(ElemType)
  nnodes = GiveNnodes(ElemType)
  if(ProblemType == "Thermal"):
    nrows = ndime
    ncols = nnodes
  elif(ProblemType == "Structural_Mechanics"):
    if (ndime == 1):
      nrows = 1
      ncols = nnodes
    elif(ndime == 2):
      nrows = 3 #This case not include revolution solid case
      ncols = nnodes * 2
      
  B_matrix = np.zeros((nrows, ncols), dtype=float)

  if(ProblemType == "Thermal"):
    for idime in range(ndime):
      for inode in range(nnodes):
        B_matrix[idime,inode] = derivCart[idime,inode]
  elif(ProblemType == "Structural_Mechanics"):
    if (ndime == 1):
      for idime in range(ndime):
        for inode in range(nnodes):
          B_matrix[idime,inode] = derivCart[idime,inode]
    elif(ndime == 2):
      for inode in range(nnodes):
        col_ini = inode * ndime
        B_matrix[0, col_ini] = derivCart[0, inode]
        B_matrix[1, col_ini + 1] = derivCart[1, inode]
        B_matrix[2, col_ini] = derivCart[1, inode]
        B_matrix[2, col_ini + 1] = derivCart[0, inode]
  return B_matrix

def giveElemVolume(elem, elemType):
  mat_id = elem["MaterialId"] - 1
  material = globalvars.data["Materials"][mat_id]
  if(elemType == "BR02" or elemType == "TRUSS02"):
    area = material["Area"]
    node1 = elem["Connectivities"][0] - 1
    x1 = globalvars.data["Mesh"]["Nodes"][node1]["x"]
    y1 = globalvars.data["Mesh"]["Nodes"][node1]["y"]
    node2 = elem["Connectivities"][1] - 1
    x2 = globalvars.data["Mesh"]["Nodes"][node2]["x"]
    y2 = globalvars.data["Mesh"]["Nodes"][node2]["y"]
    length = math.sqrt((x2 - x1)** 2 + (y2 - y1)** 2)
    volume = area*length
  elif(elemType == "TR03" or elemType == "QU04"):
    thickness = 1.0
    if(material["Plane_Type"] == "Plane_Stress"):
        thickness = material["Thickness"]
    ngauss = Set_Ngauss(elemType)
    Nodes = globalvars.data["Mesh"]["Nodes"]
    pos_pg, W = GaussQuadrature(elemType)
    area = 0.0
    for igauss in range(ngauss):
      det_Jac, _ = derivCartesian(elem, elemType, Nodes, pos_pg[igauss])
      area += det_Jac * W[igauss]
    volume = area * thickness
  return volume



