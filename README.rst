This is the code for obtaining the results of the paper *Schwarz Methods for Nonlocal Problems* 
by M. Schuster, C. Vollmann and V. Schulz that can be found on https://arxiv.org/abs/2405.01905.

Build and Install on Ubuntu
===========================
In order to clone the project do
::
  git clone https://github.com/schustermatthias/nlschwarz.git path/to/local_folder

| Since this code contains a customized version of **nlfem** the following **basic requirements** of nlfem are needed
| ``gcc, g++, python3-dev, python3-venv, libgmp-dev, libcgal-dev, metis, libmetis-dev, libarmadillo-dev``.
On Ubuntu this can be done via
::
  sudo apt-get install git gcc g++ libarmadillo-dev liblapack-dev libmetis-dev
  sudo apt-get install python3-venv python3-dev libgmp-dev libcgal-dev

| See https://gitlab.uni-trier.de/pde-opt/nonlocal-models/nlfem for more information.
| Moreover to run nlshape **legacy FEniCS(version 2019.1.0)** is required. In order to use FEniCS in a virtual environment, it may has to be 
installed globally and then inherited as a global site package. A virtual environment can be built and activated via
::
  mkdir venv
  python3 -m venv venv/
  source venv/bin/activate

Additionally the packages from the file **requirements.txt** are neccessary and can be installed by
::
  (venv) python3 -m pip install -r requirements.txt

The creation of the virtual environment and the installation of packages from requirements.txt can probably also be done via your IDE.
Finally, nlfem can be installed by
::
  (venv) python3 setup.py build --force install
  
Running the Examples from the Paper
===================================
In order to run an example from the paper just uncomment the corresponding line in main.py and execute main.py. If you would like to produce the results of one single example, execute examples.dirichlet_problem(...), examples.dirichlet_5_domains(...), examples.neumann_problem(...) or examples.patch_test(...).
Here, adjustments can be made in the configuration file that is used as an argument. Moreover, since different mesh sizes were tested, 
there exist a different mesh file for each mesh size, which can be chosen by changing the ending of the utilized mesh file 
(i.e., "mesh/dirichlet_problem_010" for mesh size h=0.1 or "mesh/dirichlet_problem_005" for h=0.05). See inside the folder "mesh" for all available meshes. 
The h-convergence results can be reproduced by running the associated functions starting with examples.compare...(...).

To run one of the (preconditioned) GMRES tests, the function example_gmres(...) has to be executed. Here, add one of the arguments preconditioner="jacobi" 
or preconditioner="gauss-seidel" to chose GMRES with the corresponding preconditioner. Leave this argument out to utilize plain GMRES. Again, the mesh size can be adjusted as described above and
there is also a configuration file (conf_gmres.py) where all remaining parameters are stored and can be changed. 
Lastly, the function examples.compute_condition_numbers_gmres(...) computes the condition number of the system matrix for each GMRES version (with mesh file and conf_gmres as arguments that can again be adjusted). 


Raw Data
========
The raw data of the examples of the Dirichlet and the Neumann problem as well as the patch tests can be found under nlschwarz/results/paper. The results of the GMRES examples can be found in the paper and easily be reproduced as described above.

License
=======
nlschwarz is published under GNU General Public License version 3. Copyright (c) 2025 Matthias Schuster

| Parts of the project are taken from **nlfem** and have been customized.
| nlfem is published under GNU General Public License version 3. Copyright (c) 2021 Manuel Klar, Christian Vollmann
  
