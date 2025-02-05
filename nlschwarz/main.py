import os
import datetime

import examples
import conf_gmres
import conf_neumann_problem
import conf_dirichlet_problem
import conf_dirichlet_5_domains
import conf_patch_test

if __name__ == '__main__':
    # mesh_data.convert_msh_to_xdmf("nlschwarz/mesh/5_domains_001")
    source = [1., 1.]
    max_iterations = 2000
    save_history = 1

    timestamp = datetime.datetime.now().strftime("%Y_%m_%d") + "/" + datetime.datetime.now().strftime("%H_%M")
    results_folder = "results/" + timestamp + "/"
    #cwd = os.getcwd()
    #results_folder = os.path.join(cwd, results_folder)
    if not os.path.exists("./" + results_folder):
        os.makedirs("./" + results_folder)

    # examples.dirichlet_problem(max_iterations, save_history, results_folder, "mesh/dirichlet_problem_005", conf_dirichlet_problem)
    # examples.compare_dirichlet_problem_solutions(max_iterations, save_history, results_folder, conf_dirichlet_problem)

    # examples.dirichlet_5_domains(max_iterations, save_history, results_folder, "mesh/dirichlet_5_domains_010", conf_dirichlet_5_domains)
    # examples.compare_dirichlet_5_domains_solutions(max_iterations, save_history, results_folder, conf_dirichlet_5_domains)

    # examples.neumann_problem(max_iterations, save_history, results_folder, "mesh/neumann_2_domains_0025", conf_neumann_problem)
    # examples.compare_neumann_problem_solutions(max_iterations, save_history, results_folder, conf_neumann_problem)

    # examples.patch_test(max_iterations, save_history, results_folder, "mesh/init_shape_005", conf_patch_test)
    # examples.compare_patch_test_solutions(max_iterations, save_history, results_folder, conf_patch_test)

    examples.example_gmres(results_folder, "mesh/init_shape_005", conf_gmres, preconditioner="jacobi")
    # examples.compute_condition_numbers_gmres("mesh/init_shape_005", conf_gmres)

