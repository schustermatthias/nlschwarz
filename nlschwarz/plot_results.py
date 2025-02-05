import os
import matplotlib
import matplotlib.pyplot as plt

import meshio
import numpy as np

import helper


def write_xdmf_time_series(filename, points, cells, data, time_steps):
    with meshio.xdmf.TimeSeriesWriter(filename) as writer:
        writer.write_points_cells(points, [meshio.CellBlock(type="triangle", data=cells)])
        for t in time_steps:
            writer.write_data(t, point_data={"solution": data[t]})


def convert_solutions_to_xdmf_time_series_old(filename, solutions, problem_data, max_iterations):
    """Create xdmf-time-series by concatenating vertices, cells and solutions for all subdomains."""
    number_subdomains = len(problem_data)
    first_mesh = problem_data[0]["mesh"]
    vertices = np.array(first_mesh["vertices"])
    cells = np.concatenate((np.array(first_mesh["elements"][first_mesh["elementLabels"] < 0]),
                           np.array(first_mesh["elements"][first_mesh["elementLabels"] == 1])))
    data = [solutions[t][0] for t in range(max_iterations)]

    for domain in range(1, number_subdomains):
        domainLabel = domain + 1
        shift = len(vertices)
        mesh = problem_data[domain]["mesh"]
        vertices = np.concatenate((vertices, mesh["vertices"]))
        cells = np.concatenate((cells, mesh["elements"][mesh["elementLabels"] == domainLabel] + shift))
        data = [np.concatenate((data[t], solutions[t][domain])) for t in range(max_iterations)]

    write_xdmf_time_series(filename, vertices, cells, data, range(max_iterations))


def interpolate_and_convert_solutions_to_xdmf_time_series(filename, solutions, vertices_domains, vertexLabels_domains,
                                                          mesh):
    number_solutions = len(solutions)
    number_subdomains = len(vertices_domains)
    number_pictures = number_solutions * number_subdomains
    vertices = mesh["vertices"]
    vertexLabels = mesh["vertexLabels"]

    number_vertices = len(vertices)
    data = [np.zeros(number_vertices) for k in range(number_pictures)]
    k = 0
    for t in range(number_solutions):
        for domain in range(number_subdomains):
            # Interpolate the solution in iteration t regarding the subdomain "domain" for the given mesh
            solution = solutions[t][domain][vertexLabels_domains[domain] > 0]
            solution_interpolated = helper.linear_interpolation(
                                        vertices_domains[domain][vertexLabels_domains[domain] > 0],
                                        solution, vertices[vertexLabels > 0])
            data[k][vertexLabels > 0] = solution_interpolated
            k += 1

    cells = mesh["elements"]
    write_xdmf_time_series(filename, vertices, cells, data, range(number_pictures))


def convert_solutions_to_xdmf_time_series(filename, solutions, mesh):
    # Without interpolation. For every subproblem the same mesh was used.
    number_solutions = len(solutions)
    number_subdomains = len(solutions[0])
    number_pictures = number_solutions * number_subdomains
    number_vertices = len(solutions[0][0])
    data = np.reshape(solutions, (number_pictures, number_vertices))

    vertices = mesh["vertices"]
    cells = mesh["elements"]
    write_xdmf_time_series(filename, vertices, cells, data, range(number_pictures))


def convert_single_solution_to_xdmf(filename, solution, mesh):
    vertices = np.array(mesh["vertices"])
    cells = mesh["elements"]
    data = solution

    write_xdmf_time_series(filename, vertices, cells, [data], range(1))


def series_plot(results_folder, solutions, vertices, cells, plot_type):
    """ plot_type can be chosen from 'color', 'contour' or 'surface' """
    plot_functions = {
        "contour": contour_plot,
        "color": color_plot,
        "surface": surface_plot
    }
    plot_function = plot_functions[plot_type]
    number_iterations = len(solutions)
    number_subproblems = len(solutions[0])
    os.makedirs(results_folder + "/" + plot_type + "/")
    file_start = results_folder + "/" + plot_type + "/" + "plot_"
    for i in range(number_iterations):
        for j in range(number_subproblems):
            filename = file_start + str(i) + "_" + str(j) + ".pdf"
            plot_function(filename, solutions[i][j], vertices, cells)


def contour_plot(filename, solution, vertices, cells):
    plt.tricontourf(vertices[:, 0], vertices[:, 1], cells, solution, cmap=matplotlib.colormaps["turbo"])
    plt.triplot(vertices[:, 0], vertices[:, 1], cells, lw=0.1, color="white", alpha=0.2)
    plt.colorbar()
    plt.savefig(filename)
    plt.show()


def color_plot(filename, solution, vertices, cells):
    plt.tripcolor(vertices[:, 0], vertices[:, 1], cells, solution, cmap=matplotlib.colormaps["turbo"])
    plt.triplot(vertices[:, 0], vertices[:, 1], cells, lw=0.1, color="white", alpha= 0.1)
    plt.colorbar()
    plt.savefig(filename)
    plt.show()


def surface_plot(filename, solution, vertices, cells):
    figure = plt.figure()
    ax = figure.add_subplot(projection='3d', title="Solution u")
    ax.plot_trisurf(vertices[:, 0], vertices[:, 1], cells, solution, cmap=matplotlib.colormaps["turbo"],
                    antialiased=False)
    plt.savefig(filename)
    plt.show()

def plot_error(filename, error_vector, ylabel, x_values = None):
    if x_values is None:
        plt.plot(error_vector)
    else:
        plt.plot(x_values, error_vector)
    plt.xlabel('number of iterations')
    plt.ylabel(ylabel)
    plt.yscale('log')
    # plt.xscale('log')
    plt.savefig(filename)
    plt.show()

def plot_all(results_folder, mesh, solutions, energy_error=None, residual=None):
    color_plot(results_folder + "color_plot.pdf", solutions[-1][-1], mesh["vertices"],
               mesh["elements"])
    contour_plot(results_folder + "contour_plot.pdf", solutions[-1][-1], mesh["vertices"],
                 mesh["elements"])
    surface_plot(results_folder + "surface_plot.pdf", solutions[-1][-1], mesh["vertices"],
                              mesh["elements"])
    #convert_single_solution_to_xdmf(results_folder + "paraview_plot.xdmf", solutions[-1][-1],
    #                                             mesh)
    #convert_solutions_to_xdmf_time_series(results_folder + "solutions.xdmf", solutions,
    #                                                   mesh)

    file_dir = results_folder + "additional_results.txt"
    if energy_error is not None:
        energy_error = energy_error[1:]
        x_values = [i for i in range(1, len(energy_error) + 1)]
        plot_error(results_folder + "energy_error.pdf", energy_error, "energy error", x_values=x_values)
        add_data_to_file(file_dir, "energy_error", energy_error)
    if residual is not None:
        plot_error(results_folder + "residual.pdf", residual, "residual error")
        add_data_to_file(file_dir, "residual", residual)


def add_data_to_file(file_dir, name, data_array):
    file = open(file_dir, "a")
    file.write(name + ": " + str(data_array) + "\n")
    file.close()
