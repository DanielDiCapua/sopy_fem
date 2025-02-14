import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib as mpl
import sopy_fem.globalvars as globalvars
from sopy_fem.initialization import initialization
from sopy_fem.bars_utils import ElemTrusses_IntForces, ElemBars_IntFluxes
from sopy_fem.plot_solids import plot_contour

def plotMesh(figTitle="Mesh"):
    window = plt.figure(figTitle)
    nodes = globalvars.data["Mesh"]["Nodes"]
    coords = np.zeros((len(nodes), 2), dtype=float)
    for i, node in enumerate(nodes):
        coords[i, 0] = node["x"]
        coords[i, 1] = node["y"]

    elements = globalvars.data["Mesh"]["Elements"]
    elemList= np.zeros((len(elements), 2), dtype=int)
    for i, elem in enumerate(elements):
        elemList[i, :] = elem["Connectivities"]

    fig, ax = plot_nodes(coords, color='k')
    fig, ax = plot_elems(coords, elemList, fig, ax, color='C7')
    label2d(coords, elemList, ax)
    plt.title("Mesh")
    window.tight_layout()

def plotDeformed(u_vec, figTitle):
    window = plt.figure(figTitle)
    nodes = globalvars.data["Mesh"]["Nodes"]
    coords = np.zeros((len(nodes), 2), dtype=float)
    scale_disp = globalvars.scale_disp
    if ("Deformed_scale" in globalvars.data["Postprocess"] and globalvars.data["Postprocess"]["Deformed_scale"]):
        scale_disp = globalvars.data["Postprocess"]["Deformed_scale"]
    for i, node in enumerate(nodes):
        idire_x = globalvars.madgln[i, 0]
        idire_y = globalvars.madgln[i, 1]
        disp_x = u_vec[idire_x]
        disp_y = u_vec[idire_y]
        coords[i, 0] = node["x"] + scale_disp*disp_x
        coords[i, 1] = node["y"] + scale_disp*disp_y

    elements = globalvars.data["Mesh"]["Elements"]
    E= np.zeros((len(elements), 2), dtype=int)
    for i, elem in enumerate(elements):
        E[i, :] = elem["Connectivities"]

    fig, ax = plot_nodes(coords, color='k')
    fig, ax = plot_elems(coords, E, fig, ax, color='C7')
    label2d(coords, E, ax)
    plt.title(figTitle.replace("_"," "))
    window.tight_layout()
    
def plotElemIntFluxes(figTitle, plotTitle):
    if(figTitle =="Axial Forces"):
        ElemIntFluxes = ElemTrusses_IntForces()
    elif (figTitle == "Thermal Fluxes" or figTitle == "Current Intensity"):
        ElemIntFluxes = ElemBars_IntFluxes()
    
    window = plt.figure(figTitle)
    nodes = globalvars.data["Mesh"]["Nodes"]
    coords = np.zeros((len(nodes), 2), dtype=float)
    for i, node in enumerate(nodes):
        coords[i, 0] = node["x"]
        coords[i, 1] = node["y"]

    elements = globalvars.data["Mesh"]["Elements"]
    elemList= np.zeros((len(elements), 2), dtype=int)
    for i, elem in enumerate(elements):
        elemList[i, :] = elem["Connectivities"]

    fig, ax = plot_nodes(coords, color='k')

    plt.title(plotTitle)
    contour_limit = np.zeros((2), dtype=float)
    min_flux = ElemIntFluxes.min()
    max_flux = ElemIntFluxes.max()
    if (abs(min_flux - max_flux) < 1.0e-10):
        if (min_flux < 0):
            max_flux = 0.0
        if (max_flux > 0):
            min_flux = 0.0
    contour_limit[0] = min_flux
    contour_limit[1] = max_flux

    con = plotTitle, ElemIntFluxes, contour_limit
    fig, ax = plot_elems(coords, elemList, fig, ax, color='C7', contour=con)
    labelBarValue(coords, elemList, ElemIntFluxes, ax)
    window.tight_layout()

def plot_nodes(coords, color='k', size=10):
    plt.scatter(coords[:,0],coords[:,1],marker='o',s=15*size,color=color,zorder=30)
    fig = plt.gcf()
    ax = plt.gca()
    return fig, ax

def plot_elems(coords, E, fig=None, ax=None, color='C0', contour=None, lim_scale=1.2):
 
    if contour is not None:
        clabel = contour[0]
        contour_lim = contour[2]
        contour = contour[1]
        
    for k,e in enumerate(E):
        NA = int(e[-2])
        NE = int(e[-1])
        
        if contour is not None:
            # build up colormap and normalizer
            colormap = plt.get_cmap('bwr', 10)
            colormap.set_over('0.95')
            colormap.set_under('0.6')
            
            if contour_lim is None:
                norm = mpl.colors.Normalize(vmin=min(contour), vmax=max(contour))
            else:
                norm = mpl.colors.Normalize(vmin=contour_lim[0], vmax=contour_lim[1])
            color = colormap(norm(contour[k]))
            
            # create a ScalarMappable and initialize a data structure
            s_m = mpl.cm.ScalarMappable(cmap=colormap, norm=norm)
            s_m.set_array([])
            
        if contour is not None:
                _, = plt.plot([coords[NA-1][0],coords[NE-1][0]],
                                    [coords[NA-1][1],coords[NE-1][1]],color='C7',zorder=1,linewidth=9.0)
        _, = plt.plot([coords[NA-1][0],coords[NE-1][0]],
                            [coords[NA-1][1],coords[NE-1][1]],color=color,zorder=2,linewidth=5.0)
        
    if contour is not None:
        plt.colorbar(s_m, ax=ax, ticks=np.linspace(contour_lim[0],contour_lim[1],11))
        plt.ylabel('CONTOUR = '+clabel)
        
    m = lim_scale*max(abs(plt.xlim()[0]),
                    abs(plt.xlim()[1]),
                    abs(plt.ylim()[0]),
                    abs(plt.ylim()[1]))
        
    if lim_scale < 0:
        m = -lim_scale
    else:
        alimx = (plt.xlim()[1] - plt.xlim()[0])/2
        alimy = (plt.ylim()[1] - plt.ylim()[0])/2
        mlimx = (plt.xlim()[1] + plt.xlim()[0])/2
        mlimy = (plt.ylim()[1] + plt.ylim()[0])/2
        alim = lim_scale*max(alimx,alimy)
        minx, maxx = (mlimx-alim, mlimx+alim)
        miny, maxy = (mlimy-alim, mlimy+alim)
        if abs(lim_scale) < 100:
            plt.xlim(minx, maxx)
            plt.ylim(miny, maxy)

    plt.gca().set_aspect('equal')
    if lim_scale < 0 and abs(lim_scale) < 100:
        plt.xlim(-m,m)
        plt.ylim(-m,m)

    plt.gcf().set_size_inches(8,6)
    
    return fig, ax

def label2d(coords, elemList, axes):
    """Plot the nodes and element label
    """
    x_min = coords[:,0].min()
    x_max = coords[:,0].max()
    y_min = coords[:,1].min()
    y_max = coords[:, 1].max()
    max_size = max((x_max - x_min), (y_max - y_min))
    offset_fact = 0.02
    offset_x = offset_fact * max_size
    offset_y = offset_fact * max_size
    if (x_min == x_max):
        offset_y = 0.0
    if (y_min == y_max):
        offset_x = 0.0

    for inode, node in enumerate(coords):
        x_pos = node[0] + offset_x
        y_pos = node[1] + offset_y
        axes.text( x_pos , y_pos, str(inode+1), color='b', size=10)

    for ele, con in enumerate(elemList):
        xm = (coords[con[0]-1][0] + coords[con[1]-1][0]) / 2 + offset_x
        ym = (coords[con[0]-1][1] + coords[con[1]-1][1]) / 2 + offset_y
        axes.text(xm, ym, str(ele + 1), color='g', size=10)
        
def labelBarValue(coords, elemList, valueList, axes):
    """Plot text with bar value
    """
    x_min = coords[:,0].min()
    x_max = coords[:,0].max()
    y_min = coords[:,1].min()
    y_max = coords[:, 1].max()
    max_size = max((x_max - x_min), (y_max - y_min))
    offset_fact = 0.02
    offset_x = offset_fact * max_size
    offset_y = offset_fact * max_size
    if (x_min == x_max):
        offset_y = 0.0
    if (y_min == y_max):
        offset_x = 0.0

    for ielem, con in enumerate(elemList):
        xm = (coords[con[0]-1][0] + coords[con[1]-1][0]) / 2 + offset_x
        ym = (coords[con[0] - 1][1] + coords[con[1] - 1][1]) / 2 + offset_y
        #valueFormatted = round(valueList[ielem],4)
        valueFormatted =  '{0:.5g}'.format(valueList[ielem])
        axes.text(xm, ym, valueFormatted, color='g', size=12)

def plotNodalBarResult(figTitle, plotTitle, nodalValues):
    window = plt.figure(figTitle)
    quadsMesh = createFictitiousQuadsMeshFromBars(globalvars.data["Mesh"])
    quadNodalValues = createQuadMeshNodalValues(nodalValues, globalvars.data["Mesh"]["Elements"])
    nodes = globalvars.data["Mesh"]["Nodes"]
    coords = np.zeros((len(nodes), 2), dtype=float)
    for i, node in enumerate(nodes):
        coords[i, 0] = node["x"]
        coords[i, 1] = node["y"]
    plot_nodes(coords, color='k')
    plot_contour(quadsMesh, quadNodalValues, figTitle, plotTitle)
    window.tight_layout()

def createFictitiousQuadsMeshFromBars(barMesh):
    barNodesList = barMesh["Nodes"]
    barElemsList = barMesh["Elements"]

    x_vec = np.zeros((len(barNodesList)), dtype=float)
    y_vec = np.zeros((len(barNodesList)), dtype=float)
    for inode, node in enumerate(barNodesList):
        x_vec[inode] = node["x"]
        y_vec[inode] = node["y"]
    x_min = x_vec.min()
    x_max = x_vec.max()
    y_min = y_vec.min()
    y_max = y_vec.max()
    max_size = max((x_max - x_min), (y_max - y_min))
    elementHeight = 0.04 * max_size

    quadsMesh = {
        "ElemType": "QU04",
        "Nodes": [],
        "Elements": []
    }

    for ielem, barElem in enumerate(barElemsList):
        node1 = barElem["Connectivities"][0] - 1
        x1 = barNodesList[node1]["x"]
        y1 = barNodesList[node1]["y"]
        node2 = barElem["Connectivities"][1] - 1
        x2 = barNodesList[node2]["x"]
        y2 = barNodesList[node2]["y"]
        l = math.sqrt((x2 - x1)** 2 + (y2 - y1)** 2)
        cos_alpha = (x2 - x1) / l
        sin_alpha = (y2 - y1) / l
        nx = -sin_alpha
        ny = cos_alpha

        materialId = barElem["MaterialId"]
        quadElem = {
            "Connectivities": [],
            "MaterialId": materialId
        }

        node_quad_1 = {
            "x": x1 - nx*elementHeight / 2,
            "y": y1 - ny*elementHeight / 2
        }
        quadsMesh["Nodes"].append(node_quad_1)

        node_quad_2 = {
            "x": x2 - nx*elementHeight / 2,
            "y": y2 - ny*elementHeight / 2
        }
        quadsMesh["Nodes"].append(node_quad_2)

        
        node_quad_3 = {
            "x": x2 + nx*elementHeight / 2,
            "y": y2 + ny*elementHeight / 2
        }
        quadsMesh["Nodes"].append(node_quad_3)

        node_quad_4 = {
            "x": x1 + nx*elementHeight / 2,
            "y": y1 + ny*elementHeight / 2
        }
        quadsMesh["Nodes"].append(node_quad_4)

        inode_quad_1 = ielem * 4 + 1
        quadElem["Connectivities"].append(inode_quad_1)
        inode_quad_2 = ielem * 4 + 2
        quadElem["Connectivities"].append(inode_quad_2)
        inode_quad_3 = ielem * 4 + 3
        quadElem["Connectivities"].append(inode_quad_3)
        inode_quad_4 = ielem * 4 + 4
        quadElem["Connectivities"].append(inode_quad_4)

        quadsMesh["Elements"].append(quadElem)
    
    return quadsMesh

def createQuadMeshNodalValues(barNodalValues, barElemsList):
    quadNodalValues = np.zeros((4*len(barElemsList)),dtype=float)
    for ielem, barElem in enumerate(barElemsList):
        nodeBar_1 = barElem["Connectivities"][0] - 1
        valueBarNode1 = barNodalValues[nodeBar_1]
        nodeBar_2 = barElem["Connectivities"][1] - 1
        valueBarNode2 = barNodalValues[nodeBar_2]

        nodeQuad_1 = ielem * 4
        quadNodalValues[nodeQuad_1] = valueBarNode1

        nodeQuad_2 = ielem * 4 + 1
        quadNodalValues[nodeQuad_2] = valueBarNode2

        nodeQuad_3 = ielem * 4 + 2
        quadNodalValues[nodeQuad_3] =valueBarNode2

        nodeQuad_4 = ielem * 4 + 3
        quadNodalValues[nodeQuad_4] = valueBarNode1

    return quadNodalValues
