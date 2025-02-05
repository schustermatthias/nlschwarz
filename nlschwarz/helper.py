import copy
import numpy as np

from interpolator import LinearNDInterpolatorExt


def get_conf_dirichlet_from_file(file):
    confs = file.confs
    kernels = file.kernels
    loads = file.loads
    return confs, kernels, loads


def get_conf_neumann_from_file(file):
    confs = file.confs
    kernels = file.kernels
    loads = file.loads
    kappa_list = file.kappa_list
    return confs, kernels, loads, kappa_list


def load_vector(mesh, weights, points, function, subdomainLabel=None):
    vertices = mesh["vertices"]
    number_vertices = len(vertices)
    load = np.zeros(number_vertices)
    if subdomainLabel is None:
        elements = mesh["elements"]
    else:
        elements = mesh["elements"][mesh["elementLabels"] == subdomainLabel]

    basis_function_values = basis_function_value(points.transpose())

    for i in range(len(elements)):
        triangle_i = elements[i]
        triangle_i_vertices = np.array([vertices[triangle_i[0]], vertices[triangle_i[1]], vertices[triangle_i[2]]])
        matrix_i = np.array([triangle_i_vertices[1] - triangle_i_vertices[0],
                             triangle_i_vertices[2] - triangle_i_vertices[0]]).transpose()
        det_matrix_i = matrix_i[0, 0] * matrix_i[1, 1] - matrix_i[1, 0] * matrix_i[0, 1]
        points_transformed = np.array([triangle_i_vertices[0] + matrix_i @ points[i] for i in range(len(points))])
        function_values = function(points_transformed.transpose())
        for a in range(3):
            integral_value = (function_values * basis_function_values[a] * weights).sum()
            load[triangle_i[a]] += abs(det_matrix_i) * integral_value

    return load


def basis_function_value(point):
    return np.array([1. - point[0] - point[1], point[0], point[1]])


def linear_interpolation(points, data, new_points):
    #interpolator = LinearNDInterpolator(points, data)
    interpolator = LinearNDInterpolatorExt(points, data)
    return interpolator(new_points)


def save_dict(dictionary, name):
    file = open(name + ".txt", "w")
    for key in dictionary.keys():
        file.write(str(key) + ": " + str(dictionary[key]) + "\n")
    file.close()


def create_configuration_dict(mesh_file, kernels, confs, loads):
    dictionary = {"mesh": mesh_file,
                  "kernels": kernels,
                  "confs": confs,
                  "loads or source": loads}
    return dictionary


def trim_artificial_node_from_mesh(mesh):
    mesh["elements"] = mesh["elements"][mesh["elementLabels"] != 0]
    mesh["vertices"] = mesh["vertices"][0:-1]
    mesh["vertexLabels"] = mesh["vertexLabels"][0:-1]
    mesh["elementLabels"] = mesh["elementLabels"][mesh["elementLabels"] != 0]

    return mesh

def trim_artificial_node_from_solution_series(solutions):
    iterations = len(solutions)
    subproblems = len(solutions[0])
    trimmed_solutions = [[solutions[i][j][0:-1] for j in range(subproblems)] for i in range(iterations)]
    return trimmed_solutions


def convert_to_cg_function(function, boundary_map):
    number_doubled_vertices = len(boundary_map)
    cg_function = copy.deepcopy(function[0:-number_doubled_vertices])
    for doubled_vertices in boundary_map:
        cg_function[doubled_vertices[0]] = 0.5 * (function[doubled_vertices[0]] + function[doubled_vertices[1]])
    return cg_function


class gmres_counter(object):
    def __init__(self, disp=True):
        self._disp = disp
        self.niter = 0

    def __call__(self, rk=None):
        self.niter += 1
        if self._disp:
            print('iteration %3i\tresidual = %s' % (self.niter, str(rk)))