import numpy as np

import helper
import mesh_data
import plot_results
from schwarz import Schwarz
import schwarz


def neumann_problem(max_iterations, save_history, results_folder, mesh_file, conf_file):
    boundary_map, mesh = mesh_data.get_artificial_node_meshes(mesh_file, problem="neumann")
    neumann_boundary_map = [[2, 1], [4, 3]]
    confs, kernels, loads, kappa_list = helper.get_conf_neumann_from_file(conf_file)

    weights = loads[0]["weights"]
    points = loads[0]["points"]
    load_function_1 = loads[0]["function"]
    load_function_2 = loads[1]["function"]

    load_1 = helper.load_vector(mesh, weights, points, load_function_1, subdomainLabel=1)
    load_2 = helper.load_vector(mesh, weights, points, load_function_2, subdomainLabel=3)

    dictionary = helper.create_configuration_dict(mesh_file, kernels, confs, loads)
    helper.save_dict(dictionary, results_folder + "configuration")

    subproblem_1 = dict(mesh=mesh, kernel=kernels[0], conf=confs[0], load_vector=load_1)
    subproblem_2 = dict(mesh=mesh, kernel=kernels[1], conf=confs[1], load_vector=load_2)
    subdomain_labels = [1, 3]

    problem_data = [subproblem_1, subproblem_2]
    problem = Schwarz(problem_data, cholesky=0, operator="diffusion", boundary_map=boundary_map,
                      nonlocal_neumann_problem=True, neumann_boundary_map=neumann_boundary_map, kappa_list=kappa_list,
                      subdomain_labels=subdomain_labels)

    def boundary_function(x):
        return 0.0

    solutions = problem.block_iterative_method(max_iterations, save_history, method="multiplicative",
                                               boundary_function=boundary_function, print_error=1, residual_tol=1e-9)


    first_mesh = problem_data[0]["mesh"]
    first_mesh["vertexLabels"] = problem.vertexLabels[0]
    first_mesh = helper.trim_artificial_node_from_mesh(first_mesh)
    trimmed_solutions = helper.trim_artificial_node_from_solution_series(solutions)

    plot_results.plot_all(results_folder, first_mesh, trimmed_solutions, problem.energy_error, problem.residual)
    final_solution = helper.convert_to_cg_function(trimmed_solutions[-1][-1], boundary_map)
    return final_solution


def compare_neumann_problem_solutions(max_iterations, save_history, results_folder, conf_file):
    confs, kernels, loads, kappa_list = helper.get_conf_neumann_from_file(conf_file)
    weights = loads[0]["weights"]
    points = loads[0]["points"]
    mesh_file_01 = "mesh/neumann_2_domains_010"
    mesh_file_0075 = "mesh/neumann_2_domains_0075"
    mesh_file_005 = "mesh/neumann_2_domains_005"
    mesh_file_0025 = "mesh/neumann_2_domains_0025"
    mesh_file_0015 = "mesh/neumann_2_domains_0015"
    mesh_file_001 = "mesh/neumann_2_domains_001"

    solution_01 = neumann_problem(max_iterations, save_history, results_folder, mesh_file_01, conf_file)
    solution_0075 = neumann_problem(max_iterations, save_history, results_folder, mesh_file_0075, conf_file)
    solution_005 = neumann_problem(max_iterations, save_history, results_folder, mesh_file_005, conf_file)
    solution_0025 = neumann_problem(max_iterations, save_history, results_folder, mesh_file_0025, conf_file)
    solution_0015 = neumann_problem(max_iterations, save_history, results_folder, mesh_file_0015, conf_file)
    solution_001 = neumann_problem(max_iterations, save_history, results_folder, mesh_file_001, conf_file)

    mesh_01 = mesh_data.get_mesh_dict(mesh_file_01)
    mesh_0075 = mesh_data.get_mesh_dict(mesh_file_0075)
    mesh_005 = mesh_data.get_mesh_dict(mesh_file_005)
    mesh_0025 = mesh_data.get_mesh_dict(mesh_file_0025)
    mesh_0015 = mesh_data.get_mesh_dict(mesh_file_0015)
    mesh_001 = mesh_data.get_mesh_dict(mesh_file_001)

    points_01 = mesh_01["vertices"]
    points_0075 = mesh_0075["vertices"]
    points_005 = mesh_005["vertices"]
    points_0025 = mesh_0025["vertices"]
    points_0015 = mesh_0015["vertices"]
    points_001 = mesh_001["vertices"]

    solution_01_interpolated = helper.linear_interpolation(points_01, solution_01, points_001)
    difference_01 = solution_01_interpolated - solution_001

    print('Compute difference between h=0.1 to h=0.01')
    norm = L2_norm(mesh_001, difference_01, weights, points)
    print('L2 norm:' + str(norm))

    solution_0075_interpolated = helper.linear_interpolation(points_0075, solution_0075, points_001)
    difference_0075 = solution_0075_interpolated - solution_001

    print('Compute difference between h=0.075 to h=0.01')
    norm = L2_norm(mesh_001, difference_0075, weights, points)
    print('L2 norm:' + str(norm))

    solution_005_interpolated = helper.linear_interpolation(points_005, solution_005, points_001)
    difference_005 = solution_005_interpolated - solution_001

    print('Compute difference between h=0.05 to h=0.01')
    norm = L2_norm(mesh_001, difference_005, weights, points)
    print('L2 norm:' + str(norm))

    solution_0025_interpolated = helper.linear_interpolation(points_0025, solution_0025, points_001)
    difference_0025 = solution_0025_interpolated - solution_001

    print('Compute difference between h=0.025 to h=0.01')
    norm = L2_norm(mesh_001, difference_0025, weights, points)
    print('L2 norm:' + str(norm))

    solution_0015_interpolated = helper.linear_interpolation(points_0015, solution_0015, points_001)
    difference_0015 = solution_0015_interpolated - solution_001

    print('Compute difference between h=0.015 to h=0.01')
    norm = L2_norm(mesh_001, difference_0015, weights, points)
    print('L2 norm:' + str(norm))


def compute_condition_numbers_gmres(mesh_file, conf_file):

    boundary_map, mesh = mesh_data.get_artificial_node_meshes(mesh_file)

    confs, kernels, loads_conf = helper.get_conf_dirichlet_from_file(conf_file)

    subdomain_labels = [1, 2]
    load_list = get_load_vector(mesh, loads_conf, subdomain_labels)

    problem_data = create_problem_data(mesh, kernels, confs, load_list, subdomain_labels)

    problem = Schwarz(problem_data, cholesky=0, operator="diffusion", boundary_map=boundary_map,
                      subdomain_labels=subdomain_labels)
    cond_number = problem.get_condition_number()
    print("Condition number without preconditioner: " + str(cond_number))
    cond_number = problem.get_condition_number(method='additive')
    print("Condition number with block-Jacobi-preconditioner: " + str(cond_number))
    cond_number = problem.get_condition_number(method='multiplicative')
    print("Condition number with block-Gauss-Seidel-preconditioner: " + str(cond_number))


def example_gmres(results_folder, mesh_file, conf_file, preconditioner=None):
    boundary_map, mesh = mesh_data.get_artificial_node_meshes(mesh_file)

    confs, kernels, loads_conf = helper.get_conf_dirichlet_from_file(conf_file)

    subdomain_labels = [1, 2]
    load_list = get_load_vector(mesh, loads_conf, subdomain_labels)

    dictionary = helper.create_configuration_dict(mesh_file, kernels, confs, loads_conf)
    helper.save_dict(dictionary, results_folder + "configuration")

    problem_data = create_problem_data(mesh, kernels, confs, load_list, subdomain_labels)

    problem = Schwarz(problem_data, cholesky=0, operator="diffusion", boundary_map=boundary_map,
                      subdomain_labels=subdomain_labels)

    if preconditioner == "jacobi":
        print("GMRES with Jacobi-preconditioner is used.")
        solution = problem.preconditioned_gmres(method="jacobi")
    elif preconditioner == "gauss-seidel":
        print("GMRES with Gauss-Seidel-preconditioner is used.")
        solution = problem.preconditioned_gmres(method="gauss-seidel")
    else:
        print("GMRES without preconditioner is used.")
        solution = problem.preconditioned_gmres()

    first_mesh = problem_data[0]["mesh"]
    first_mesh["vertexLabels"] = problem.vertexLabels[0]
    first_mesh = helper.trim_artificial_node_from_mesh(first_mesh)
    solution = solution[0:-1]

    # plot_results.series_plot(results_folder, trimmed_solutions, first_mesh["vertices"], first_mesh["elements"], "color")
    plot_results.color_plot(results_folder + "color_plot.pdf", solution, first_mesh["vertices"],
                            first_mesh["elements"])
    plot_results.contour_plot(results_folder + "contour_plot.pdf", solution, first_mesh["vertices"],
                              first_mesh["elements"])
    plot_results.surface_plot(results_folder + "surface_plot.pdf", solution, first_mesh["vertices"],
                              first_mesh["elements"])

    plot_results.convert_single_solution_to_xdmf(results_folder + "paraview_plot.xdmf", solution,
                                                 first_mesh)


def dirichlet_5_domains(max_iterations, save_history, results_folder, mesh_file, conf_file):

    def boundary_function(x):
        return 0.0

    subdomain_labels = [1, 2, 3, 4, 5]

    solve_dirichlet_problem(max_iterations, save_history, results_folder, conf_file, mesh_file, boundary_function,
                            subdomain_labels)

def dirichlet_problem(max_iterations, save_history, results_folder, mesh_file, conf_file):

    def boundary_function(x):
        return 0.0

    subdomain_labels = [1, 2, 3]

    solve_dirichlet_problem(max_iterations, save_history, results_folder, conf_file, mesh_file, boundary_function,
                            subdomain_labels)


def compare_dirichlet_problem_solutions(max_iterations, save_history, results_folder, conf_file):
    mesh_file_010 = "mesh/dirichlet_problem_010"
    mesh_file_0075 = "mesh/dirichlet_problem_0075"
    mesh_file_005 = "mesh/dirichlet_problem_005"
    mesh_file_0025 = "mesh/dirichlet_problem_0025"
    mesh_file_0015 = "mesh/dirichlet_problem_0015"
    mesh_file_001 = "mesh/dirichlet_problem_001"

    def boundary_function(x):
        return 0.0 # x[0] + x[1]

    subdomain_labels = [1, 2, 3]

    solution_01 = solve_dirichlet_problem(max_iterations, save_history, results_folder, conf_file, mesh_file_010,
                                          boundary_function, subdomain_labels)
    solution_0075 = solve_dirichlet_problem(max_iterations, save_history, results_folder, conf_file, mesh_file_0075,
                                            boundary_function, subdomain_labels)
    solution_005 = solve_dirichlet_problem(max_iterations, save_history, results_folder, conf_file, mesh_file_005,
                                           boundary_function, subdomain_labels)
    solution_0025 = solve_dirichlet_problem(max_iterations, save_history, results_folder, conf_file, mesh_file_0025,
                                            boundary_function, subdomain_labels)
    solution_0015 = solve_dirichlet_problem(max_iterations, save_history, results_folder, conf_file, mesh_file_0015,
                                            boundary_function, subdomain_labels)
    solution_001 = solve_dirichlet_problem(max_iterations, save_history, results_folder, conf_file, mesh_file_001,
                                           boundary_function, subdomain_labels)

    mesh_01 = mesh_data.get_mesh_dict(mesh_file_010)
    mesh_0075 = mesh_data.get_mesh_dict(mesh_file_0075)
    mesh_005 = mesh_data.get_mesh_dict(mesh_file_005)
    mesh_0025 = mesh_data.get_mesh_dict(mesh_file_0025)
    mesh_0015 = mesh_data.get_mesh_dict(mesh_file_0015)
    mesh_001 = mesh_data.get_mesh_dict(mesh_file_001)

    points_01 = mesh_01["vertices"]
    points_0075 = mesh_0075["vertices"]
    points_005 = mesh_005["vertices"]
    points_0025 = mesh_0025["vertices"]
    points_0015 = mesh_0015["vertices"]
    points_001 = mesh_001["vertices"]
    weights = conf_file.confs[0]["quadrature"]["outer"]["weights"]
    points = conf_file.confs[0]["quadrature"]["outer"]["points"]

    solution_01_interpolated = helper.linear_interpolation(points_01, solution_01, points_001)
    difference_01 = solution_01_interpolated - solution_001

    print('Compute difference between h=0.1 to h=0.01')
    norm = L2_norm(mesh_001, difference_01, weights, points)
    print('L2 norm:' + str(norm))

    solution_0075_interpolated = helper.linear_interpolation(points_0075, solution_0075, points_001)
    difference_0075 = solution_0075_interpolated - solution_001

    print('Compute difference between h=0.075 to h=0.01')
    norm = L2_norm(mesh_001, difference_0075, weights, points)
    print('L2 norm:' + str(norm))

    solution_005_interpolated = helper.linear_interpolation(points_005, solution_005, points_001)
    difference_005 = solution_005_interpolated - solution_001

    print('Compute difference between h=0.05 to h=0.01')
    norm = L2_norm(mesh_001, difference_005, weights, points)
    print('L2 norm:' + str(norm))

    solution_0025_interpolated = helper.linear_interpolation(points_0025, solution_0025, points_001)
    difference_0025 = solution_0025_interpolated - solution_001

    print('Compute difference between h=0.025 to h=0.01')
    norm = L2_norm(mesh_001, difference_0025, weights, points)
    print('L2 norm:' + str(norm))

    solution_0015_interpolated = helper.linear_interpolation(points_0015, solution_0015, points_001)
    difference_0015 = solution_0015_interpolated - solution_001

    print('Compute difference between h=0.015 to h=0.01')
    norm = L2_norm(mesh_001, difference_0015, weights, points)
    print('L2 norm:' + str(norm))


def compare_dirichlet_5_domains_solutions(max_iterations, save_history, results_folder, conf_file):
    mesh_file_010 = "mesh/dirichlet_5_domains_010"
    mesh_file_0075 = "mesh/dirichlet_5_domains_0075"
    mesh_file_005 = "mesh/dirichlet_5_domains_005"
    mesh_file_0025 = "mesh/dirichlet_5_domains_0025"
    mesh_file_0015 = "mesh/dirichlet_5_domains_0015"
    mesh_file_001 = "mesh/dirichlet_5_domains_001"

    def boundary_function(x):
        return 0.0 # x[0] + x[1]

    subdomain_labels = [1, 2, 3, 4, 5]

    solution_01 = solve_dirichlet_problem(max_iterations, save_history, results_folder, conf_file, mesh_file_010,
                                          boundary_function, subdomain_labels)
    solution_0075 = solve_dirichlet_problem(max_iterations, save_history, results_folder, conf_file, mesh_file_0075,
                                            boundary_function, subdomain_labels)
    solution_005 = solve_dirichlet_problem(max_iterations, save_history, results_folder, conf_file, mesh_file_005,
                                           boundary_function, subdomain_labels)
    solution_0025 = solve_dirichlet_problem(max_iterations, save_history, results_folder, conf_file, mesh_file_0025,
                                            boundary_function, subdomain_labels)
    solution_0015 = solve_dirichlet_problem(max_iterations, save_history, results_folder, conf_file, mesh_file_0015,
                                            boundary_function, subdomain_labels)
    solution_001 = solve_dirichlet_problem(max_iterations, save_history, results_folder, conf_file, mesh_file_001,
                                           boundary_function, subdomain_labels)

    mesh_01 = mesh_data.get_mesh_dict(mesh_file_010)
    mesh_0075 = mesh_data.get_mesh_dict(mesh_file_0075)
    mesh_005 = mesh_data.get_mesh_dict(mesh_file_005)
    mesh_0025 = mesh_data.get_mesh_dict(mesh_file_0025)
    mesh_0015 = mesh_data.get_mesh_dict(mesh_file_0015)
    mesh_001 = mesh_data.get_mesh_dict(mesh_file_001)

    points_01 = mesh_01["vertices"]
    points_0075 = mesh_0075["vertices"]
    points_005 = mesh_005["vertices"]
    points_0025 = mesh_0025["vertices"]
    points_0015 = mesh_0015["vertices"]
    points_001 = mesh_001["vertices"]
    weights = conf_file.confs[0]["quadrature"]["outer"]["weights"]
    points = conf_file.confs[0]["quadrature"]["outer"]["points"]

    solution_01_interpolated = helper.linear_interpolation(points_01, solution_01, points_001)
    difference_01 = solution_01_interpolated - solution_001

    print('Compute difference between h=0.1 to h=0.01')
    norm = L2_norm(mesh_001, difference_01, weights, points)
    print('L2 norm:' + str(norm))

    solution_0075_interpolated = helper.linear_interpolation(points_0075, solution_0075, points_001)
    difference_0075 = solution_0075_interpolated - solution_001

    print('Compute difference between h=0.075 to h=0.01')
    norm = L2_norm(mesh_001, difference_0075, weights, points)
    print('L2 norm:' + str(norm))

    solution_005_interpolated = helper.linear_interpolation(points_005, solution_005, points_001)
    difference_005 = solution_005_interpolated - solution_001

    print('Compute difference between h=0.05 to h=0.01')
    norm = L2_norm(mesh_001, difference_005, weights, points)
    print('L2 norm:' + str(norm))

    solution_0025_interpolated = helper.linear_interpolation(points_0025, solution_0025, points_001)
    difference_0025 = solution_0025_interpolated - solution_001

    print('Compute difference between h=0.025 to h=0.01')
    norm = L2_norm(mesh_001, difference_0025, weights, points)
    print('L2 norm:' + str(norm))

    solution_0015_interpolated = helper.linear_interpolation(points_0015, solution_0015, points_001)
    difference_0015 = solution_0015_interpolated - solution_001

    print('Compute difference between h=0.015 to h=0.01')
    norm = L2_norm(mesh_001, difference_0015, weights, points)
    print('L2 norm:' + str(norm))


def patch_test(max_iterations, save_history, results_folder, mesh_file, conf_file):

    def boundary_function(x):
        return x[0] + x[1]

    subdomain_labels = [1, 2]

    solve_dirichlet_problem(max_iterations, save_history, results_folder, conf_file, mesh_file, boundary_function,
                            subdomain_labels)


def compare_patch_test_solutions(max_iterations, save_history, results_folder, conf_file):
    mesh_file_010 = "mesh/init_shape_010"
    mesh_file_0075 = "mesh/init_shape_0075"
    mesh_file_005 = "mesh/init_shape_005"
    mesh_file_0025 = "mesh/init_shape_0025"
    mesh_file_0015 = "mesh/init_shape_0015"
    mesh_file_001 = "mesh/init_shape_001"

    def boundary_function(x):
        return x[0] + x[1]

    subdomain_labels = [1, 2]

    solution_01 = solve_dirichlet_problem(max_iterations, save_history, results_folder, conf_file, mesh_file_010,
                                          boundary_function, subdomain_labels)
    solution_0075 = solve_dirichlet_problem(max_iterations, save_history, results_folder, conf_file, mesh_file_0075,
                                            boundary_function, subdomain_labels)
    solution_005 = solve_dirichlet_problem(max_iterations, save_history, results_folder, conf_file, mesh_file_005,
                                           boundary_function, subdomain_labels)
    solution_0025 = solve_dirichlet_problem(max_iterations, save_history, results_folder, conf_file, mesh_file_0025,
                                            boundary_function, subdomain_labels)
    solution_0015 = solve_dirichlet_problem(max_iterations, save_history, results_folder, conf_file, mesh_file_0015,
                                            boundary_function, subdomain_labels)
    solution_001 = solve_dirichlet_problem(max_iterations, save_history, results_folder, conf_file, mesh_file_001,
                                           boundary_function, subdomain_labels)

    mesh_01 = mesh_data.get_mesh_dict(mesh_file_010)
    mesh_0075 = mesh_data.get_mesh_dict(mesh_file_0075)
    mesh_005 = mesh_data.get_mesh_dict(mesh_file_005)
    mesh_0025 = mesh_data.get_mesh_dict(mesh_file_0025)
    mesh_0015 = mesh_data.get_mesh_dict(mesh_file_0015)
    mesh_001 = mesh_data.get_mesh_dict(mesh_file_001)

    points_01 = mesh_01["vertices"]
    points_0075 = mesh_0075["vertices"]
    points_005 = mesh_005["vertices"]
    points_0025 = mesh_0025["vertices"]
    points_0015 = mesh_0015["vertices"]
    points_001 = mesh_001["vertices"]
    weights = conf_file.confs[0]["quadrature"]["outer"]["weights"]
    points = conf_file.confs[0]["quadrature"]["outer"]["points"]

    g_001 = np.apply_along_axis(boundary_function, axis=1, arr=points_001)

    solution_01_interpolated = helper.linear_interpolation(points_01, solution_01, points_001)
    difference_010 = solution_01_interpolated - g_001

    print('Compute difference between h=0.1 to interpolated exact solution')
    norm = L2_norm(mesh_001, difference_010, weights, points)
    print('L2 norm:' + str(norm))

    solution_0075_interpolated = helper.linear_interpolation(points_0075, solution_0075, points_001)
    difference_0075 = solution_0075_interpolated - g_001

    print('Compute difference between h=0.075 to interpolated exact solution')
    norm = L2_norm(mesh_001, difference_0075, weights, points)
    print('L2 norm:' + str(norm))

    solution_005_interpolated = helper.linear_interpolation(points_005, solution_005, points_001)
    difference_005 = solution_005_interpolated - g_001

    print('Compute difference between h=0.05 to interpolated exact solution')
    norm = L2_norm(mesh_001, difference_005, weights, points)
    print('L2 norm:' + str(norm))

    solution_0025_interpolated = helper.linear_interpolation(points_0025, solution_0025, points_001)
    difference_0025 = solution_0025_interpolated - g_001

    print('Compute difference between h=0.025 to interpolated exact solution')
    norm = L2_norm(mesh_001, difference_0025, weights, points)
    print('L2 norm:' + str(norm))

    solution_0015_interpolated = helper.linear_interpolation(points_0015, solution_0015, points_001)
    difference_0015 = solution_0015_interpolated - g_001

    print('Compute difference between h=0.015 to interpolated exact solution')
    norm = L2_norm(mesh_001, difference_0015, weights, points)
    print('L2 norm:' + str(norm))

    difference_001 = solution_001 - g_001

    print('Compute difference between h=0.01 to interpolated exact solution')
    norm = L2_norm(mesh_001, difference_001, weights, points)
    print('L2 norm:' + str(norm))


def solve_dirichlet_problem(max_iterations, save_history, results_folder, conf_file, mesh_file, boundary_function,
                            subdomain_labels):
    confs, kernels, loads_conf = helper.get_conf_dirichlet_from_file(conf_file)

    boundary_map, mesh = mesh_data.get_artificial_node_meshes(mesh_file)

    load_list = get_load_vector(mesh, loads_conf, subdomain_labels)

    dictionary = helper.create_configuration_dict(mesh_file, kernels, confs, loads_conf)
    helper.save_dict(dictionary, results_folder + "configuration")


    problem_data = create_problem_data(mesh, kernels, confs, load_list, subdomain_labels)

    problem = Schwarz(problem_data, cholesky=0, operator="diffusion", boundary_map=boundary_map,
                      subdomain_labels=subdomain_labels)

    solutions = problem.block_iterative_method(max_iterations, save_history, method="multiplicative",
                                               boundary_function=boundary_function, print_error=1, residual_tol=1e-9)

    first_mesh = problem_data[0]["mesh"]
    first_mesh["vertexLabels"] = problem.vertexLabels[0]
    first_mesh = helper.trim_artificial_node_from_mesh(first_mesh)
    trimmed_solutions = helper.trim_artificial_node_from_solution_series(solutions)

    plot_results.plot_all(results_folder, first_mesh, trimmed_solutions, problem.energy_error, problem.residual)

    final_solution = helper.convert_to_cg_function(trimmed_solutions[-1][-1], boundary_map)
    return final_solution


def L2_norm(mesh, vec, weights, points):
    vec_squared = np.square(vec)
    identity_matrix = schwarz.identity(mesh, weights=weights, points=points)
    norm_squared = np.dot(vec_squared, identity_matrix @ vec_squared)
    return np.sqrt(norm_squared)


def create_problem_data(mesh, kernels, confs, load_list, subdomain_labels):
    problem_data = []
    for i in range(len(subdomain_labels)):
        problem_data.append(dict(mesh=mesh, kernel=kernels[i], conf=confs[i], load_vector=load_list[i]))
    return problem_data

def get_load_vector(mesh, loads_conf, subdomain_labels):
    load_list = []
    for i in range(len(subdomain_labels)):
        weights = loads_conf[i]["weights"]
        points = loads_conf[i]["points"]
        load_function = loads_conf[i]["function"]
        load = helper.load_vector(mesh, weights, points, load_function, subdomainLabel=subdomain_labels[i])
        load_list.append(load)

    return load_list
