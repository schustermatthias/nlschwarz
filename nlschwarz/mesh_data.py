import nlfem
import numpy as np

import meshio


def read_gmsh(filename, method='CG', problem='dirichlet'):
    if problem not in {"dirichlet", "neumann"}:
        raise ValueError("Wrong argument for the variable 'problem'. Use 'dirichlet' or 'neumann' instead.")
    meshio_mesh = meshio.read(filename)
    elements = np.array(meshio_mesh.cells_dict["triangle"], dtype=np.long)
    vertices = np.array(meshio_mesh.points[:, :2])
    elementLabels = np.array(meshio_mesh.cell_data_dict["gmsh:physical"]["triangle"], dtype=np.long)

    try:
        lines = np.array(meshio_mesh.cells_dict["line"], dtype=int)
        # The following line may be needed in the future, if we have multiple boundaries (maybe including the boundary
        # of Omega)
        lineLabels = np.array(meshio_mesh.cell_data_dict["gmsh:physical"]["line"], dtype=np.long)
    except KeyError:
        lines = None
        lineLabels = None
    number_elements = len(elements)
    if method == "CG":
        vertexLabels = np.zeros(len(vertices), dtype=np.long)
        for i in range(number_elements):
            label = elementLabels[i]
            for k in elements[i]:
                vertexLabels[k] = np.max((vertexLabels[k], label))
    elif method == "DG":
        vertexLabels = np.zeros(3 * len(elements), dtype=np.int)
        for i in range(number_elements):
            vertexLabels[3 * i:3 * i + 2] = elementLabels[i]
    else:
        raise ValueError("Wrong argument for the variable 'method'. Use 'DG' or 'CG' instead.")

    if problem == "dirichlet":
        boundary_label = np.max(elementLabels)
        elementLabels[elementLabels == boundary_label] = -1.0
        vertexLabels[vertexLabels == boundary_label] = -1.0

    return elements, vertices, elementLabels, vertexLabels, lines, lineLabels


def get_mesh_dict(mesh_file, method='CG', problem='dirichlet'):
    elements, vertices, elementLabels, vertexLabels, lines, lineLabels = read_gmsh(mesh_file + ".msh", method, problem)
    mesh = {"elements": elements,
            "vertices": vertices,
            "elementLabels": elementLabels,
            "vertexLabels": vertexLabels,
            "lines": lines,
            "lineLabels": lineLabels,
            "outdim": 1}
    return mesh


def get_artificial_node_meshes(mesh_object, method='CG', problem="dirichlet"):
    """In this method the mesh-dictionary is set up. Furthermore, the interface-nodes are doubled and an artificial node
    is added. Elements are changed such that every element, that contains an interface node, now addresses the version
    of this node that corresponds to the label of this element. The two interfaces are now connected through cells with
    label 0 that additionally contain the artificial node."""
    if type(mesh_object) == str:
        elements, vertices, elementLabels, vertexLabels, lines, lineLabels = read_gmsh(mesh_object + ".msh", method,
                                                                                       problem=problem)
    elif type(mesh_object) == dict:
        elements, elementLabels = mesh_object["elements"], mesh_object["elementLabels"]
        vertices, vertexLabels = mesh_object["vertices"], mesh_object["vertexLabels"]
        lines, lineLabels = mesh_object["lines"], mesh_object["lineLabels"]
    else:
        raise ValueError("The variable mesh_object either has to be the name of an .msh file(without msh) or a "
                         "dictionary with the required information.")
    artificial_node = [0.0, 0.0]

    boundary_indices = np.unique(lines)

    # Substract boundary nodes that are on the Dirichlet boundary, extract Dirichlet boundary from lines
    inner_boundary_vertices = np.where([vertexLabels[boundary_indices[i]] > 0.0 for i in range(len(boundary_indices))])
    inner_boundary_vertices = np.array(inner_boundary_vertices[0])
    if len(inner_boundary_vertices) < len(boundary_indices):
        # TODO > or >=?
        boundary_line_indices = np.where([vertexLabels[lines[i][0]] > 0.0 or vertexLabels[lines[i][1]] > 0.0
                                          for i in range(len(lines))])
        lines = lines[boundary_line_indices]
        lineLabels = lineLabels[boundary_line_indices]

        boundary_indices = boundary_indices[inner_boundary_vertices]

    number_vertices = len(vertices)
    number_boundary_nodes = len(boundary_indices)

    # get indices of the new boundary nodes
    new_boundary_indices = np.arange(number_vertices, number_vertices + number_boundary_nodes)
    # create map with the information, which boundary nodes correspond to each other
    boundary_map = np.concatenate((boundary_indices.reshape((number_boundary_nodes, 1)),
                                   new_boundary_indices.reshape((number_boundary_nodes, 1))), axis=1, dtype=np.int)

    # add new vertices
    boundary_vertices = vertices[boundary_indices]
    boundary_vertices = np.append(boundary_vertices, [artificial_node], axis=0)
    vertices = np.append(vertices, boundary_vertices, axis=0)
    artificial_node_index = len(vertices) - 1

    # change interface nodes of elements with label 2 to the new interface nodes
    unique_lineLabels = np.unique(lineLabels)
    for label in unique_lineLabels:
        current_line_vertices = np.unique(lines[lineLabels == label])
        first_line = lines[lineLabels == label][0]

        # returns list of elements, which include the vertices of the first line
        elements_first_line = np.where([first_line[0] in elements[i] and first_line[1] in elements[i]
                                        for i in range(len(elementLabels))])
        max_elementLabel = np.max(elementLabels[elements_first_line])

        element_list = np.where(elementLabels == max_elementLabel)[0]
        for index in element_list:
            element = elements[index]
            for i in range(3):
                # both checks are necessary since we only want to double the current line without the dirichlet nodes
                # (dirichlet nodes can be in current_line_vertices but not in boundary_indices)
                if element[i] in current_line_vertices and element[i] in boundary_indices:
                    index_boundary_list = np.where(boundary_indices == element[i])
                    elements[index][i] = boundary_map[index_boundary_list[0], 1]

    # add new elements with label 0 to connect the two boundaries that are the same
    number_new_elements = len(lines)
    new_elements = []
    for index in range(number_new_elements):
        vertex_1_index = np.where(boundary_indices == lines[index][0])[0]
        if len(vertex_1_index) == 0:
            # vertex_1 is on the Dirichlet boundary and therefore is not doubled
            vertex_1 = lines[index][0]
        else:
            vertex_1 = boundary_map[vertex_1_index[0], 1]

        vertex_2_index = np.where(boundary_indices == lines[index][1])[0]
        if len(vertex_2_index) == 0:
            # vertex_2 is on the Dirichlet boundary and therefore is not doubled
            vertex_2 = lines[index][0]
        else:
            vertex_2 = boundary_map[vertex_2_index[0], 1]

        new_elements.append([vertex_1, vertex_2])
    new_elements = np.array(new_elements, dtype=np.int)
    new_elements = np.concatenate((lines, new_elements), axis=0, dtype=np.int)
    artificial_node_array = np.ones((2*number_new_elements, 1), dtype=np.int) * artificial_node_index
    new_elements = np.concatenate((new_elements, artificial_node_array), axis=1)
    elements = np.concatenate((elements, new_elements), axis=0, dtype=np.int)
    new_elementLabels = np.zeros(len(new_elements), dtype=np.int)
    elementLabels = np.concatenate((elementLabels, new_elementLabels), dtype=np.int)
    vertexLabels = nlfem.get_vertexLabel(elements, elementLabels, vertices)

    mesh = {"elements": elements,
            "vertices": vertices,
            "elementLabels": elementLabels,
            "vertexLabels": vertexLabels,
            "lines": lines,
            "lineLabels": lineLabels,
            "outdim": 1}
    return boundary_map, mesh

