import math
import numpy as np

def GiveNdime(ElemType):
    if(ElemType == "BAR02" or ElemType == "BAR03"):
        ndime=1
    elif(ElemType == "TRUSS02" or ElemType == "TR03" or ElemType == "TR06" or ElemType == "QU04" or ElemType == "QU08" or ElemType == "QU09"):
      ndime=2
    return ndime

def Set_Ngauss(ElemType):
    if(ElemType == "BAR02" or ElemType == "TRUSS02"):
      ngauss=2
    elif(ElemType == "BAR03" or ElemType == "TR03"):
      ngauss=3
    elif(ElemType == "TR06"):
      ngauss=6
    elif(ElemType == "QU04"):
      ngauss=4
    elif(ElemType == "QU08" or ElemType == "QU09"):
      ngauss=9
    return ngauss

def GaussQuadrature(ElemType):
  ngauss = Set_Ngauss(ElemType)
  ndime = GiveNdime(ElemType)
  W = np.zeros((ngauss), dtype=float)
  pos = np.zeros((ngauss,ndime), dtype=float)
  if(ElemType == "BAR02" or ElemType == "TRUSS02"):
    pos[0,0] = -math.sqrt(3) / 3
    W[0] = 1.0
    pos[1,0] = math.sqrt(3) / 3
    W[1] = 1.0
  elif(ElemType == "BAR03"):
    pos[0,0] = -math.sqrt(3 / 5)
    W[0] = 5 / 9
    pos[1,0] = 0.0
    W[1] = 8 / 9
    pos[2,0] = math.sqrt(3 / 5)
    W[2] = 5 / 9
  elif(ElemType == "TR03"):
    pos[0,0] = 0.5
    pos[0,1] = 0
    W[0] = 1 / 6
    pos[1,0] = 0.5
    pos[1,1] = 0.5
    W[1] = 1 / 6
    pos[2,0] = 0
    pos[2,1] = 0.5
    W[2] = 1 / 6
  elif(ElemType == "TR06"):
    pos[0,0] = 0.0915762135
    pos[0,1] = 0.0915762135
    W[0] = 0.0549758719
    pos[1,0] = 0.8168475730
    pos[1,1] = 0.0915762135
    W[1] = 0.0549758719
    pos[2,0] = 0.0915762135
    pos[2,1] = 0.8168475730
    W[2] = 0.0549758719
    pos[3,0] = 0.4459484909
    pos[3,1] = 0.4459484909
    W[3] = 0.1116907949
    pos[4,0] = 0.1081030182
    pos[4,1] = 0.4459484909
    W[4] = 0.1116907949
    pos[5,0] = 0.4459484909
    pos[5,1] = 0.1081030182
    W[5] = 0.1116907949
  elif(ElemType == "QU04"):
    pos[0,0] = -math.sqrt(3) / 3
    pos[0,1] = -math.sqrt(3) / 3
    W[0] = 1.0
    pos[1,0] = math.sqrt(3) / 3
    pos[1,1] = -math.sqrt(3) / 3
    W[1] = 1.0
    pos[2,0] = math.sqrt(3) / 3
    pos[2,1] = math.sqrt(3) / 3
    W[2] = 1.0
    pos[3,0] = -math.sqrt(3) / 3
    pos[3,1] = math.sqrt(3) / 3
    W[3] = 1.0
  elif(ElemType == "QU08" or ElemType == "QU09"):
    pos[0,0] = -math.sqrt(3 / 5)
    pos[0,1] = -math.sqrt(3 / 5)
    W[0] = 25 / 81
    pos[1,0] = math.sqrt(3 / 5)
    pos[1,1] = -math.sqrt(3 / 5)
    W[1] = 25 / 81
    pos[2,0] = math.sqrt(3 / 5)
    pos[2,1] = math.sqrt(3 / 5)
    W[2] = 25 / 81
    pos[3,0] = -math.sqrt(3 / 5)
    pos[3,1] = math.sqrt(3 / 5)
    W[3] = 25 / 81
    pos[4,0] = 0
    pos[4,1] = -math.sqrt(3 / 5)
    W[4] = 40 / 81
    pos[5,0] = math.sqrt(3 / 5)
    pos[5,1] = 0
    W[5] = 40 / 81
    pos[6,0] = 0
    pos[6,1] = math.sqrt(3 / 5)
    W[6] = 40 / 81
    pos[7,0] = math.sqrt(3 / 5)
    pos[7,1] = 0
    W[7] = 40 / 81
    pos[8,0] = 0
    pos[8,1] = 0
    W[8] = 64 / 81
  return pos, W 
  