import globalvars
import numpy as np

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
      ndof = 2
  elif (ProblemType == "Electrical"):
    if(ElemType == "BAR02" or ElemType == "BAR03"):
      ndof=1
  return ndof

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
  elif (ProblemType == "Electrical"):
    if(ElemType == "BAR02" or ElemType == "BAR03"):
      ncomp=1
  return ncomp

def initialization():
    mesh_type = globalvars.data["Mesh"]["ElemType"]
    ProblemType = globalvars.data["ProblemType"]
    globalvars.ndof = GiveNdof(mesh_type, ProblemType)
    globalvars.ncomp = GiveNComp(mesh_type, ProblemType)
    num_nodes = len(globalvars.data["Mesh"]["Nodes"])

    globalvars.madgln = np.full((num_nodes, globalvars.ndof), -2, dtype=int)

    nglpr = 0
    for ipres in range(len(globalvars.data["Constraints"])):
        id_node = globalvars.data["Constraints"][ipres]["Node"] - 1
        for igl in range(globalvars.ndof):
            if(globalvars.data["Constraints"][ipres]["Activation"][igl]):
                nglpr=nglpr+1
                globalvars.madgln[id_node, igl] = -1  #Flag value
    
    #Calculate madgln
    globalvars.neq = 0
    for elem in globalvars.data["Mesh"]["Elements"]:
      conec_list = elem["Connectivities"]
      for inode in range(len(conec_list)):
        id_node = conec_list[inode]-1
        for igl in range(globalvars.ndof):
          if (globalvars.madgln[id_node,igl] == -2):
            globalvars.madgln[id_node, igl] = globalvars.neq
            globalvars.neq += 1
            
    for ipres in range(len(globalvars.data["Constraints"])):
        id_node = globalvars.data["Constraints"][ipres]["Node"] - 1
        for igl in range(globalvars.ndof):
            if(globalvars.data["Constraints"][ipres]["Activation"][igl]):
                globalvars.madgln[id_node, igl] = globalvars.neq
                globalvars.neq += 1

    #Number of unknows
    globalvars.num_unknows = globalvars.neq - nglpr
    globalvars.astiff = np.zeros((globalvars.neq, globalvars.neq), dtype=float)
    globalvars.asload = np.zeros((globalvars.neq), dtype=float)
    globalvars.u_vec = np.zeros((globalvars.neq), dtype=float)

    num_known = globalvars.neq - globalvars.num_unknows
    globalvars.u_known = np.zeros((num_known), dtype=float)
    for ipres in range(len(globalvars.data["Constraints"])):
        id_node = globalvars.data["Constraints"][ipres]["Node"] - 1
        for igl in range(globalvars.ndof):
            if (globalvars.data["Constraints"][ipres]["Activation"][igl]):
                idire = globalvars.madgln[id_node, igl] - globalvars.num_unknows
                globalvars.u_known[idire] = globalvars.data["Constraints"][ipres]["Values"][igl]
    

    

