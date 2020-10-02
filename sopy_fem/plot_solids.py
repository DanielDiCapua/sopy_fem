import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri

# converts quad elements into triangle elements
def quads_to_triangles(quadsMesh):
    trianglesMesh = {
        "ElemType": "TR03",
        "Nodes": quadsMesh["Nodes"],
        "Elements": []
    }

    for quadElem in quadsMesh["Elements"]:
        node = np.empty((4),dtype=int)
        node[:4] =  quadElem["Connectivities"][:4]
        triangleElem1 = {
            "Connectivities": [node[0], node[1], node[2]],
            "MaterialId": quadElem["MaterialId"]
        }
        trianglesMesh["Elements"].append(triangleElem1)

        triangleElem2 = {
            "Connectivities": [node[2], node[3], node[0]],
            "MaterialId": quadElem["MaterialId"]
        }
        trianglesMesh["Elements"].append(triangleElem2)
    return trianglesMesh


# plots a finite element mesh
def plot_fem_mesh(mesh):
    for element in mesh["Elements"]:
        x = []
        y = []
        for node in element["Connectivities"]:
            id_node = node - 1
            x.append(mesh["Nodes"][id_node]["x"])
            y.append(mesh["Nodes"][id_node]["y"])
        plt.fill(x, y, edgecolor='black', fill=False)

def plot_contour(mesh, nodal_values, figTitle, plotTitle):
    mesh_triangulation = mesh
    if (mesh["ElemType"] == "QU04"):
        mesh_triangulation = quads_to_triangles(mesh)
    
    x_coords = []
    y_coords = []
    for node in mesh_triangulation["Nodes"]:
        x_coords.append(node["x"])
        y_coords.append(node["y"])

    connectivities_list = []
    for element in mesh_triangulation["Elements"]:
        connectivities_elem = []
        for node in element["Connectivities"]:
            id_node = node - 1
            connectivities_elem.append(id_node)
        connectivities_list.append(connectivities_elem)

    # create an unstructured triangular grid instance
    triangulation = tri.Triangulation(x_coords, y_coords, connectivities_list)

    # plot the contours
    plt.figure(figTitle)
    plot_fem_mesh(mesh)
    plt.tricontourf(triangulation, nodal_values, 12)

    # show
    plt.title(plotTitle)
    plt.colorbar(orientation='vertical')
    plt.axis('equal')
