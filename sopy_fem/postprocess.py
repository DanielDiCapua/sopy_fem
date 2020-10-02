import numpy as np
import json
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import sopy_fem.globalvars as globalvars
from sopy_fem.initialization import initialization
from sopy_fem.solids_utils import Smooth_Strains, Smooth_Stresses, Smooth_Forces, GiveNComp
from sopy_fem.plot_solids import plot_contour
from sopy_fem.plot_trusses import plotMesh, plotDeformed, plotElemIntFluxes, plotNodalBarResult

#Set plotting defaults
gray = '#757575'
plt.rcParams['image.cmap'] = "rainbow"
plt.rcParams["mathtext.fontset"] = "cm"
plt.rcParams["text.color"] = gray
plt.rcParams["font.size"] = 10
plt.rcParams["xtick.color"] = gray
plt.rcParams["ytick.color"] = gray
plt.rcParams["axes.labelcolor"] = gray
plt.rcParams["axes.edgecolor"] = gray
plt.rcParams["axes.spines.right"] = False
plt.rcParams["axes.spines.top"] = False

def postprocess():
    ElemType = globalvars.data["Mesh"]["ElemType"]
    if("Show_displacements" in globalvars.data["Postprocess"] and globalvars.data["Postprocess"]["Show_displacements"]):
        plotDisplacements()

    if(ElemType == "TR03" or ElemType == "TR06" or ElemType == "QU04" or ElemType == "QU08" or ElemType == "QU09"):
        if("Show_strains" in globalvars.data["Postprocess"] and globalvars.data["Postprocess"]["Show_strains"]):
            plotStrains()
        
        if("Show_stresses" in globalvars.data["Postprocess"] and globalvars.data["Postprocess"]["Show_stresses"]):
            plotStresses()

        if("Show_forces" in globalvars.data["Postprocess"] and globalvars.data["Postprocess"]["Show_forces"]):
            plotForces()
    elif (ElemType == "BAR02" or ElemType == "BAR03" or ElemType == "TRUSS02"):
        plt.rcParams["axes.spines.right"] = True
        plt.rcParams["axes.spines.top"] = True
        plotMesh()
        if("Show_forces" in globalvars.data["Postprocess"] and globalvars.data["Postprocess"]["Show_forces"]):
            plotElemIntFluxes("Axial Forces", "N")

        if ("Show_temperatures" in globalvars.data["Postprocess"] and globalvars.data["Postprocess"]["Show_temperatures"]):
            plotTemperatures()

        if ("Show_voltage" in globalvars.data["Postprocess"] and globalvars.data["Postprocess"]["Show_voltage"]):
            plotVoltage()

        if("Show_thermal_fluxes" in globalvars.data["Postprocess"] and globalvars.data["Postprocess"]["Show_thermal_fluxes"]):
            plotElemIntFluxes("Thermal Fluxes", "q (W)")

        if("Show_current_intensity" in globalvars.data["Postprocess"] and globalvars.data["Postprocess"]["Show_current_intensity"]):
            plotElemIntFluxes("Current Intensity", "I (A)")                


    if("Show_reactions" in globalvars.data["Postprocess"] and globalvars.data["Postprocess"]["Show_reactions"]):
        writeReactions()

    resultsFile = globalvars.dataFileName.replace(".json", ".res.json")
    with open(resultsFile, 'w') as json_file:
        json.dump(globalvars.results, json_file, indent = 4)
    
    plt.show()

def plotDisplacements():
    mesh = globalvars.data["Mesh"]
    ElemType = mesh["ElemType"]
    ndof = globalvars.ndof
    nodal_disp_x = []
    nodal_disp_y = []
    globalvars.results["Displacements"] = []
    for inode in range(len(globalvars.data["Mesh"]["Nodes"])):
        idire_x = globalvars.madgln[inode, 0]
        disp_x = globalvars.u_vec[idire_x]
        nodal_disp_x.append(disp_x)

        if(ndof == 1):
            nodal_disp_res = {
                "Node": inode+1,
                "Disp_x": "{:10.4e}".format(disp_x),
            }
        elif(ndof == 2):
            idire_y = globalvars.madgln[inode, 1]
            disp_y = globalvars.u_vec[idire_y]
            nodal_disp_y.append(disp_y)
            nodal_disp_res = {
                "Node": inode+1,
                "Disp_x": "{:10.4e}".format(disp_x),
                "Disp_y": "{:10.4e}".format(disp_y)
            }

        globalvars.results["Displacements"].append(nodal_disp_res)
    
    if(ElemType == "TR03" or ElemType == "TR06" or ElemType == "QU04" or ElemType == "QU08" or ElemType == "QU09"):
        plot_contour(mesh, nodal_disp_x, "Displacements x", r"$u_x$ (m)")
        plot_contour(mesh, nodal_disp_y, "Displacements y", r"$u_y$ (m)")
    elif (ElemType == "BAR02" or ElemType == "BAR03"):
        plt.rcParams["axes.spines.right"] = True
        plt.rcParams["axes.spines.top"] = True
        plotNodalBarResult("Displacements", r"$u_x$ (m)", nodal_disp_x)
    elif (ElemType == "TRUSS02"):
        plt.rcParams["axes.spines.right"] = True
        plt.rcParams["axes.spines.top"] = True
        plotNodalBarResult("Horizontal displacements", r"$u_x$ (m)", nodal_disp_x)
        plotNodalBarResult("Vertical displacements", r"$u_y$ (m)", nodal_disp_y)
        if ("Show_deformed" in globalvars.data["Postprocess"] and globalvars.data["Postprocess"]["Show_deformed"]):
            plotDeformed()


def plotStrains():
    mesh = globalvars.data["Mesh"]
    ncomp = globalvars.ncomp
    nodal_Epsilon = Smooth_Strains()

    globalvars.results["Strains"] = []
    for inode in range(len(globalvars.data["Mesh"]["Nodes"])):
        if (ncomp == 1):
            nodal_strains_res = {
                "Node": inode+1,
                "Epsilon_x": "{:10.4e}".format(nodal_Epsilon[inode,0]),
            }
        elif (ncomp == 3):
            nodal_strains_res = {
                "Node": inode+1,
                "Epsilon_x": "{:10.4e}".format(nodal_Epsilon[inode,0]),
                "Epsilon_y": "{:10.4e}".format(nodal_Epsilon[inode,1]),
                "Gama_xy":   "{:10.4e}".format(nodal_Epsilon[inode,2])
            }
        globalvars.results["Strains"].append(nodal_strains_res) 

    if(ncomp == 3):
        plot_contour(mesh, nodal_Epsilon[:,0], "Strains Epsilon_x", r"$\epsilon_{x}$")
        plot_contour(mesh, nodal_Epsilon[:,1], "Strains Epsilon_y", r"$\epsilon_{y}$")
        plot_contour(mesh, nodal_Epsilon[:,2], "Strains Gama_xy", r"$\gamma_{xy}$")


def plotStresses():
    mesh = globalvars.data["Mesh"]
    ncomp = globalvars.ncomp
    nodal_Sigma = Smooth_Stresses()

    globalvars.results["Stresses"] = []
    for inode in range(len(globalvars.data["Mesh"]["Nodes"])):
        if (ncomp == 1):
            nodal_stresses_res = {
                "Node": inode+1,
                "Sigma_x": "{:10.4e}".format(nodal_Sigma[inode,0]),
            }
        elif (ncomp == 3):
            nodal_stresses_res = {
                "Node": inode+1,
                "Sigma_x": "{:10.4e}".format(nodal_Sigma[inode,0]),
                "Sigma_y": "{:10.4e}".format(nodal_Sigma[inode,1]),
                "Tau_xy":  "{:10.4e}".format(nodal_Sigma[inode,2])
            }
        globalvars.results["Stresses"].append(nodal_stresses_res)        

    if(ncomp == 3):
        plot_contour(mesh, nodal_Sigma[:,0], "Stresses Sigma_x", r"$\sigma_{x}$")
        plot_contour(mesh, nodal_Sigma[:,1], "Stresses Sigma_y", r"$\sigma_{y}$")
        plot_contour(mesh, nodal_Sigma[:, 2], "Stresses Tau_xy", r"$\tau_{xy}$")
        

def plotForces():
    mesh = globalvars.data["Mesh"]
    ncomp = globalvars.ncomp
    nodal_Forces = Smooth_Forces()

    globalvars.results["Forces"] = []
    for inode in range(len(globalvars.data["Mesh"]["Nodes"])):
        if (ncomp == 1):
            nodal_forces_res = {
                "Node": inode+1,
                "N_x": "{:10.4e}".format(nodal_Forces[inode,0]),
            }
        elif (ncomp == 3):
            nodal_forces_res = {
                "Node": inode+1,
                "N_x": "{:10.4e}".format(nodal_Forces[inode,0]),
                "N_y": "{:10.4e}".format(nodal_Forces[inode,1]),
                "N_xy":  "{:10.4e}".format(nodal_Forces[inode,2])
            }
        globalvars.results["Forces"].append(nodal_forces_res)        

    if(ncomp == 3):
        plot_contour(mesh, nodal_Forces[:,0], "Internal Forces N_x", r"$N_{x}$")
        plot_contour(mesh, nodal_Forces[:,1], "Internal Forces N_y", r"$N_{y}$")
        plot_contour(mesh, nodal_Forces[:, 2], "Internal Forces N_xy", r"$N_{xy}$")

def plotTemperatures():
    mesh = globalvars.data["Mesh"]
    ElemType = mesh["ElemType"]
    nodal_temp = []
    globalvars.results["Temperatures"] = []
    for inode in range(len(globalvars.data["Mesh"]["Nodes"])):
        idire = globalvars.madgln[inode, 0]
        temp = globalvars.u_vec[idire]
        nodal_temp.append(temp)

        nodal_temp_res = {
            "Node": inode+1,
            "Temp": "{:10.4e}".format(temp),
        }

        globalvars.results["Temperatures"].append(nodal_temp_res)
    
    if(ElemType == "TR03" or ElemType == "TR06" or ElemType == "QU04" or ElemType == "QU08" or ElemType == "QU09"):
        plot_contour(mesh, nodal_temp, "Temperatures", r"$Temp$ ($^\circ$C)")
    elif (ElemType == "BAR02" or ElemType == "BAR03"):
        plotNodalBarResult("Temperatures", r"$Temp$ ($^\circ$C)", nodal_temp)
        
def plotVoltage():
    mesh = globalvars.data["Mesh"]
    ElemType = mesh["ElemType"]
    nodal_temp = []
    globalvars.results["Voltage"] = []
    for inode in range(len(globalvars.data["Mesh"]["Nodes"])):
        idire = globalvars.madgln[inode, 0]
        temp = globalvars.u_vec[idire]
        nodal_temp.append(temp)

        nodal_temp_res = {
            "Node": inode+1,
            "Voltage": "{:10.4e}".format(temp),
        }

        globalvars.results["Voltage"].append(nodal_temp_res)
    
    if (ElemType == "BAR02" or ElemType == "BAR03"):
        plotNodalBarResult("Voltage", r"$Voltage$ (V)", nodal_temp)

def writeReactions():
    globalvars.results["Reactions"] = []
    for ipres in range(len(globalvars.data["Constraints"])):
        node = globalvars.data["Constraints"][ipres]["Node"]
        id_node = node - 1
        if (globalvars.data["ProblemType"] == "Structural_Mechanics"):
            react_results_node = {
                "Node": node,
                "Rx": "{:10.4e}".format(0.0),
                "Ry": "{:10.4e}".format(0.0)
            }
            for igl in range(globalvars.ndof):
                if(globalvars.data["Constraints"][ipres]["Activation"][igl]):
                    idire = globalvars.madgln[id_node, igl] - globalvars.num_unknows
                    if(igl == 0):
                        react_results_node["Rx"] = "{:10.4e}".format(globalvars.react_vec[idire])
                    elif(igl == 1):
                        react_results_node["Ry"] ="{:10.4e}".format(globalvars.react_vec[idire])
        elif (globalvars.data["ProblemType"] == "Thermal"):
            idire = globalvars.madgln[id_node, 0] - globalvars.num_unknows
            react_results_node = {
                "Node": node,
                "Flux": "{:10.4e}".format(globalvars.react_vec[idire])
            }
        elif (globalvars.data["ProblemType"] == "Electrical"):
            idire = globalvars.madgln[id_node, 0] - globalvars.num_unknows
            react_results_node = {
                "Node": node,
                "Current_Intensity": "{:10.4e}".format(globalvars.react_vec[idire])
            }
        
        globalvars.results["Reactions"].append(react_results_node)
        

