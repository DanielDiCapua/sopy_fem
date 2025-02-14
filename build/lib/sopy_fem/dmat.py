import numpy as np
import sopy_fem.globalvars as globalvars

def dmat_Solids2D(material):
    dmat = np.zeros((3, 3), dtype=float)
    Young = material["Young"]
    Poisson = material["Poisson"]
    G = Young/(2.0*(1.0+Poisson))

    if (material["Plane_Type"] == "Plane_Stress"):
        fact1 = Young/(1.0-Poisson**2)
        dmat[0, 0] = fact1
        dmat[0, 1] = fact1 * Poisson
        dmat[1, 0] = fact1 * Poisson
        dmat[1, 1] = fact1
        dmat[2,2] = G
    elif (material["Plane_Type"] == "Plane_Strain"):
        fact1 = Young * (1.0 - Poisson) / ((1.0 + Poisson) * (1 - 2.0 * Poisson))
        fact2= Young *  Poisson / ((1.0 + Poisson) * (1 - 2.0 * Poisson))
        dmat[0, 0] = fact1
        dmat[0, 1] = fact2
        dmat[1, 0] = fact2
        dmat[1, 1] = fact1
        dmat[2, 2] = G
        
    return dmat

def giveLocalStiffness(material, l):
    ProblemType = globalvars.data["ProblemType"]
    if (ProblemType == "Structural_Mechanics"):
        if("Young" in material and "Area" in material):
            Young = material["Young"]
            Area = material["Area"]
        else:
            raise Exception("The Area and Young definitions are mandatory for mechanics bars")
        k = Young * Area / l
    elif (ProblemType == "Thermal"):
        if("Thermal_Conductivity" in material and "Area" in material):
            Thermal_Cond = material["Thermal_Conductivity"]
            Area = material["Area"]
        else:
            raise Exception("The Area and Thermal_Conductivity definitions are mandatory for thermal bars")
        k = Thermal_Cond * Area / l
    elif (ProblemType == "Electrical"):
        if("Electrical_Conductivity" in material and "Area" in material):
            Electrical_Cond = material["Electrical_Conductivity"]
            Area = material["Area"]
        else:
            raise Exception("The Area and Electrical_Conductivity definitions are mandatory for electrical bars")
        k = Electrical_Cond * Area / l
    
    return k
        
        

        



    