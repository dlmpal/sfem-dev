# SFEM
SFEM is a C++ framework for the solution of PDEs on unstructured meshes, using the Finite Element and Finite Volume methods. It supports distributed memory parallelism via the MPI protocol. SFEM utilizes METIS for mesh partitioning, and, optionally, PETSc and SLEPc can be used for solving linear and eigenproblems respectively. However, an iterative linear solver suite is also included natively with SFEM.

## Requirements
* C++ 23 and above
* Python 3.8 and above
* CMake 3.16 and above 
* MPI (optional)
* METIS (optional if MPI is disabled)
* PETSc (optional)
* SLEPc (optional)

## Installation and usage 
Start by cloning the repository to your machine:

```
git clone https://github.com/dlmpal/sfem.git
```

Then, inside of the SFEM directory, execute the installation script (install.py). The script can be configured with several options, 
which are shown by executing:

```
python install.py --help
```

The installation script is able to download and build all required dependencies of SFEM (**--download-all**). Alternatively,
the user must provide some relevant info (e.g. **--petsc-dir** or **petsc-arch**). The user can also select whether to build the applications (**--with-apps**). The simplest way to install the package is by executing:

```
python install.py --download-all --build-all
```

By default, the downloaded dependencies are installed at **${HOME}/local/sfem-deps** and the library 
is installed at **${HOME}/local/sfem**. These locations can be changed by setting **--third-party-dir** and 
**--install-dir** respectively. 

After successfully running the script, a summary is printed and saved to file (**summary.log**). 
At the end of the summary several lines are printed, which should be copied and appended at the end of the user's **.bashrc** or equivalent file. 

### Using SFEM
The easiest way to build a program with SFEM is by using CMake. Start by creating a CMakeLists.txt as follows:

```
cmake_minimum_required(VERSION 3.16)
project(my_sfem_solver)

set(CMAKE_CXX_COMPILER $ENV{MPICXX})
set(sfem_DIR ${SFEM_INSTALL_DIR}/lib/cmake/sfem)
find_package(sfem REQUIRED)

add_executable(my_sfem_solver /path/to/my_sfem_solver.cc)
target_link_libraries(my_sfem_solver sfem::sfem)
```

The user should verify that the MPI compiler **MPICXX** used to compile my_sfem_solver.cc is the same one that was used to compile SFEM (and PETSc). 

## Non-native mesh formats
SFEM can read and convert meshes in Gmsh format using the gmshToSfem application. SFEM can also produce VTK files for visualization (e.g. via ParaView) with the sfemToVTK application. Application binaries are located under **${SFEM_INSTALL_DIR}/bin**.