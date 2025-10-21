"""
Installation script for SFEM
"""

import os
import pathlib
import zipfile
import tarfile
import argparse
import subprocess
from urllib.request import urlretrieve
from enum import Enum


def run_shell_cmd(cmd: str, cwd: str = None):
    cproc = subprocess.run(cmd, cwd=cwd, shell=True,
                           capture_output=False, text=True)
    if cproc.returncode != 0:
        msg = f'Shell command "{cmd}" exited with nonzero return code {cproc.returncode}'
        raise RuntimeError(msg)


def fetch_dependency_git(git_repo: str, dest_dir: str = None, branch: str = "master"):
    run_shell_cmd(f"git clone -b {branch} {git_repo}", cwd=dest_dir)


def fetch_dependency_url(url: str, dest_dir: str):
    tmp_filename = dest_dir + "_tmp"
    urlretrieve(url, tmp_filename)
    if tarfile.is_tarfile(tmp_filename):
        archive = tarfile.open(tmp_filename, "r")
        archive.extractall(dest_dir)
    elif zipfile.is_zipfile(tmp_filename):
        archive = zipfile.ZipFile(tmp_filename)
        archive.extractall(dest_dir)
    os.remove(tmp_filename)


def n_jobs_make():
    return os.cpu_count()


def build_dependecy(config_cmd: str = "", build_cmd: str = "",
                    install_cmd: str = "", test_cmd: str = "",
                    working_dir: str = None):
    commands = [config_cmd, build_cmd,
                install_cmd, test_cmd]

    for cmd in commands:
        run_shell_cmd(cmd, cwd=working_dir)


def get_openmpi_url_and_version():
    url = "https://download.open-mpi.org/release/open-mpi/v5.0/openmpi-5.0.6.tar.gz"
    version = "5.0.6"
    return url, version


def add_openmpi(openmpi_dir: str, f_compiler: str,
                c_compiler: str, cxx_compiler: str):
    """
    Fetch and build OpenMPI.
    """
    url, _ = get_openmpi_url_and_version()
    fetch_dependency_url(url, openmpi_dir)
    openmpi_config_cmd = f"./configure --prefix={openmpi_dir} --enable-mpi-fortran=yes "
    openmpi_config_cmd += f"FC={f_compiler} CC={c_compiler} CXX={cxx_compiler} FCFLAGS=-O3 CFLAGS=-O3 CXXCFLAGS=-O3"
    openmpi_build_cmd = f"make -j {n_jobs_make()}"
    openmpi_install_cmd = f"make install"
    build_dependecy(config_cmd=openmpi_config_cmd,
                    build_cmd=openmpi_build_cmd,
                    install_cmd=openmpi_install_cmd,
                    working_dir=os.path.join(openmpi_dir, "openmpi-5.0.6"))


def add_petsc4py_deps():
    '''
    Handle petsc4py dependencies
    '''
    # numpy
    try:
        import numpy
    except:
        run_shell_cmd("pip install numpy")

    # setuptools
    try:
        import setuptools
    except:
        run_shell_cmd("pip install setuptools")


def add_petsc(petsc_dir: str, petsc_arch: str,
              f_compiler: str, c_compiler: str, cxx_compiler: str,
              with_debugging: bool, with_petsc4py: bool, with_mpi: bool):
    '''
    Fetch and build PETSc.
    '''
    fetch_dependency_git("https://gitlab.com/petsc/petsc.git",
                         os.path.dirname(petsc_dir),
                         "release")

    if with_petsc4py:
        add_petsc4py_deps()

    petsc_config_cmd = f"./configure PETSC_ARCH={petsc_arch} "
    petsc_config_cmd += f"--with-cc={c_compiler} --with-cxx={cxx_compiler} --with-fc={f_compiler} "
    petsc_config_cmd += f"--with-mpi={int(with_mpi)} "
    if not with_debugging:
        petsc_config_cmd += "COPTFLAGS='-O3' CXXFLAGS='-O3' FOPTFLAGS='-O3' "
    petsc_config_cmd += f"--with-debugging={int(with_debugging)} "
    petsc_config_cmd += f"--with-petsc4py={int(with_petsc4py)} "
    petsc_config_cmd += "--download-f2cblaslapack=1 "
    petsc_config_cmd += "--download-scalapack=1 "
    petsc_config_cmd += "--download-mumps=1"

    petsc_build_cmd = f"make -j {n_jobs_make()} PETSC_DIR={petsc_dir} PETSC_ARCH={petsc_arch} all"
    petsc_test_cmd = f"make -j {n_jobs_make()} PETSC_DIR={petsc_dir} PETSC_ARCH={petsc_arch} check"

    build_dependecy(config_cmd=petsc_config_cmd,
                    build_cmd=petsc_build_cmd,
                    test_cmd=petsc_test_cmd,
                    working_dir=petsc_dir)


def get_petsc4py_dir(petsc_dir: str, petsc_arch: str) -> str:
    return os.path.join(petsc_dir, petsc_arch, "lib")


def add_slepc(slepc_dir: str, petsc_dir: str,
              petsc_arch: str, with_slepc4py: bool):
    '''
    Fetch and build SLEPc.
    '''
    fetch_dependency_git("https://gitlab.com/slepc/slepc.git",
                         os.path.dirname(slepc_dir),
                         "release")

    export_vars = {"SLEPC_DIR": slepc_dir,
                   "PETSC_DIR": petsc_dir,
                   "PETSC_ARCH": petsc_arch,
                   # Required only if PETSc was downloaded and built
                   "PYTHONPATH": get_petsc4py_dir(petsc_dir, petsc_arch)}

    slepc_config_cmd = ""
    for k, v in export_vars.items():
        slepc_config_cmd += f"export {k}={v};"

    slepc_config_cmd += f"./configure --with-slepc4py={int(with_slepc4py)}"
    slepc_build_cmd = f"make -j {n_jobs_make()} SLEPC_DIR={slepc_dir} PETSC_DIR={petsc_dir} PETSC_ARCH={petsc_arch} all"
    slepc_test_cmd = f"make -j {n_jobs_make()} SLEPC_DIR={slepc_dir} PETSC_DIR={petsc_dir} PETSC_ARCH={petsc_arch} check"

    build_dependecy(config_cmd=slepc_config_cmd,
                    build_cmd=slepc_build_cmd,
                    test_cmd=slepc_test_cmd,
                    working_dir=slepc_dir)


def get_slepc4py_dir(slepc_dir: str, petsc_arch: str) -> str:
    return os.path.join(slepc_dir, petsc_arch, "lib")


def add_metis(metis_dir: str, c_compiler: str):
    '''
    Fetch and build METIS, and its dependency GKlib.
    TODO CMake compatibility
    '''
    # Create the sub-directories for GKlib and METIS
    metis_subdirs = {"gklib": os.path.join(metis_dir, "GKlib"),
                     "metis": os.path.join(metis_dir, "METIS")}
    for k, v in metis_subdirs.items():
        os.makedirs(v, exist_ok=True)

    # Fetch and build GKlib
    gklib_git_repo = "https://github.com/KarypisLab/GKlib.git"
    fetch_dependency_git(gklib_git_repo, metis_dir)
    build_dependecy(config_cmd=f"make config prefix={metis_dir}",
                    build_cmd=f"make -j {n_jobs_make()} cc={c_compiler}",
                    install_cmd="make install",
                    working_dir=metis_subdirs["gklib"])

    # Fetch and build METIS
    metis_git_repo = "https://github.com/KarypisLab/METIS.git"
    fetch_dependency_git(metis_git_repo, metis_dir)
    build_dependecy(config_cmd=f"make config shared=1 prefix={metis_dir}",
                    build_cmd=f"make -j {n_jobs_make()} cc={c_compiler}",
                    install_cmd="make install",
                    working_dir=metis_subdirs["metis"])


def add_nanobind():
    '''
    Installs nanobind via pip.
    '''
    pip_cmd = 'pip install nanobind'
    run_shell_cmd(pip_cmd)


def get_pysfem_dir(install_dir: str) -> str:
    return os.path.join(install_dir, "lib", "pysfem")


def get_applications_dir(install_dir: str) -> str:
    return os.path.join(install_dir, "bin")


class FloatingPointPrecision(Enum):
    single = "single"
    double = "double"


# ==============================================================================
# Parse all arguments
parser = argparse.ArgumentParser(
    description='Build and install SFEM',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

# User's home directory
home_dir = pathlib.Path.home()

# Installation directory
parser.add_argument('--install-dir', type=str, metavar="",
                    default=os.path.join(home_dir, "local", "sfem"),
                    help='The installation directory for SFEM')

# Build directory
parser.add_argument('--build-dir', type=str, metavar="",
                    default=os.path.join(os.getcwd(), "build"),
                    help='The build directory for SFEM')

# Fortran compiler
parser.add_argument('--f-compiler', type=str,
                    metavar="", default='gfortran',
                    help='Fortran compiler. Ignored if --with-mpi is provided.')

# C compiler
parser.add_argument('--c-compiler', type=str,
                    metavar="", default='cc',
                    help='C compiler. Ignored if --with-mpi is provided.')

# C++ compiler
parser.add_argument('--cxx-compiler', type=str,
                    metavar="", default='c++',
                    help='C++ compiler. Ignored if --with-mpi is provided.')

# Floating-point arithmetic specification
parser.add_argument('--floating-point-spec', type=FloatingPointPrecision,
                    metavar="", default=FloatingPointPrecision.double,
                    help='The floating-point precision used in SFEM.')

# Build debug
parser.add_argument("--with-debugging", action="store_true",
                    help="Build with debugging enabled.")

# Build Python bindings (pysfem)
parser.add_argument('--with-pysfem', action="store_true",
                    help="Build the Python bindings for SFEM")

# Build the applications
parser.add_argument('--with-apps', action="store_true",
                    help="Build the applications")

# Build and install everything (core library, applications, Python bindings)
parser.add_argument('--build-all', action='store_true',
                    help="Build and install everything (core library, applications, Python bindings)")

# Installation directory for dependencies
parser.add_argument('--third-party-dir', type=str, metavar="",
                    default=os.path.join(home_dir, "local", "sfem-deps"),
                    help="The installation directory for all downloaded dependencies")

# Download and build all dependencies
parser.add_argument('--download-all', action='store_true',
                    help="Download all dependencies")

# MPI
parser.add_argument('--with-mpi', action='store_true',
                    help="Indicates that the C/C++ compilers are MPI-wrapped")

parser.add_argument('--mpi-f-compiler', type=str,
                    metavar="", default='mpif90',
                    help='MPI Fortran compiler')

parser.add_argument('--mpi-c-compiler', type=str,
                    metavar="", default='mpicc',
                    help='MPI C compiler')

parser.add_argument('--mpi-cxx-compiler', type=str,
                    metavar="", default='mpic++',
                    help='MPI C++ compiler')

parser.add_argument('--download-openmpi', action='store_true',
                    help="Download and build OpenMPI")

# PETSc
parser.add_argument('--petsc-dir', type=str,
                    metavar="", default=os.getenv("PETSC_DIR"),
                    help='The installation directory for PETSc')

parser.add_argument('--petsc-arch', type=str,
                    metavar="", default=os.getenv("PETSC_ARCH"),
                    help='The PETSC-ARCH')

parser.add_argument('--download-petsc', action='store_true',
                    help="Download and build PETSc")

# SLEPc
parser.add_argument('--slepc-dir', type=str,
                    metavar="", default=os.getenv("SLEPC_DIR"),
                    help='The installation directory for SLEPc')

parser.add_argument('--download-slepc', action='store_true',
                    help="Download and build SLEPc")

# METIS
parser.add_argument('--metis-dir', type=str,
                    metavar="", default=os.path.join(home_dir, "local"),
                    help='The installation directory for METIS')

parser.add_argument('--download-metis', action='store_true',
                    help="Download and build METIS")

# Parse
args = parser.parse_args()
# ==============================================================================
# Check whether --download-all or --build-all are enabled
if args.download_all:
    args.download_openmpi = True
    args.download_petsc = True
    args.download_slepc = True
    args.download_metis = True

if args.build_all:
    args.with_pysfem = True
    args.with_apps = True
# ==============================================================================
# Remove previous cmake/build/install folders
remove_prev_config = f"rm -rf cmake"
run_shell_cmd(remove_prev_config)

remove_prev_build = f"rm -rf {args.build_dir}"
run_shell_cmd(remove_prev_build)

remove_prev_install = f"rm -rf {args.install_dir}"
run_shell_cmd(remove_prev_install)
# ==============================================================================
# Handle third-party dependencies
# TODO: Check that it does not exist inside cwd
# Perhaps not needed to be created directly
os.makedirs(args.third_party_dir, exist_ok=True)

# OpenMPI
if args.download_openmpi:
    openmpi_dir = os.path.join(args.third_party_dir, "openmpi")
    add_openmpi(openmpi_dir, args.f_compiler,
                args.c_compiler, args.cxx_compiler)
    args.with_mpi = True
    args.mpi_f_compiler = os.path.join(openmpi_dir, "bin", "mpif90")
    args.mpi_c_compiler = os.path.join(openmpi_dir, "bin", "mpicc")
    args.mpi_cxx_compiler = os.path.join(openmpi_dir, "bin", "mpic++")
    args.download_petsc = True
    args.download_slepc = True

# If MPI is enabled (either manually or downloaded)
# the Fortran/C/C++ compilers are set to their MPI-wrapped counterparts
if args.with_mpi is True:
    args.f_compiler = args.mpi_f_compiler
    args.c_compiler = args.mpi_c_compiler
    args.cxx_compiler = args.mpi_cxx_compiler

# PETSc
if args.download_petsc:
    args.petsc_dir = os.path.join(args.third_party_dir, "petsc")
    args.petsc_arch = "arch-sfem"
    os.makedirs(args.petsc_dir, exist_ok=True)
    add_petsc(args.petsc_dir, args.petsc_arch,
              args.f_compiler, args.c_compiler, args.cxx_compiler,
              args.with_debugging, args.with_pysfem,
              args.with_mpi)

# Check petsc-dir
if args.petsc_dir is None or os.path.exists(args.petsc_dir) is False:
    petsc_error = "Variable petsc-dir is either not provided or invalid.\n"
    petsc_error += "Either provide a valid value, or add argument --download-petsc"
    raise RuntimeError(petsc_error)

# Check petsc-arch
if args.petsc_arch is None or os.path.exists(os.path.join(args.petsc_dir, args.petsc_arch)) is False:
    petsc_error = "Variable petsc-arch is either not provided or invalid.\n"
    petsc_error += "Either provide a valid value, or add argument --download-petsc"
    raise RuntimeError(petsc_error)

# SLEPc
if args.download_slepc:
    args.slepc_dir = os.path.join(args.third_party_dir, "slepc")
    os.makedirs(args.slepc_dir, exist_ok=True)
    add_slepc(args.slepc_dir, args.petsc_dir,
              args.petsc_arch, args.with_pysfem)

# Check slepc-dir (only if SLEPc was reqeusted)
if args.slepc_dir is not None or args.download_slepc is True:
    if os.path.exists(args.slepc_dir) is False:
        slepc_error = "Variable slepc-dir is either not provided or invalid.\n"
        slepc_error += "Either provide a valid value, or add argument --download-slepc"
        raise RuntimeError(slepc_error)

# METIS
if args.download_metis:
    args.metis_dir = os.path.join(args.third_party_dir, "metis")
    os.makedirs(args.metis_dir, exist_ok=True)
    add_metis(args.metis_dir, args.c_compiler)

# Check metis-dir
if args.metis_dir is None or os.path.exists(args.metis_dir) is False:
    metis_error = "Variable metis-dir is either not provided or invalid.\n"
    metis_error += "Either provide a valid metis-dir, or add argument --download-metis"
    raise RuntimeError(metis_error)

# nanobind
if args.with_pysfem:
    add_nanobind()
# ==============================================================================
# Create the Config.cmake.in
cmake_config_filepah = os.path.join("cmake", "Config.cmake.in")
os.makedirs(os.path.dirname(cmake_config_filepah), exist_ok=True)
with open(cmake_config_filepah, "w") as file:
    file.writelines(
        "@PACKAGE_INIT@\ninclude(${CMAKE_CURRENT_LIST_DIR}/sfemTargets.cmake)\ncheck_required_components(sfem)")
# ==============================================================================
# Run CMake
cmake_args = {
    "CMAKE_INSTALL_PREFIX": args.install_dir,
    "CMAKE_BUILD_TYPE": "RELEASE" if args.with_debugging is False else "DEBUG",
    "WITH_APPS": "On" if args.with_apps is True else "Off",
    "WITH_PYSFEM": "On" if args.with_pysfem is True else "Off",
    "CMAKE_CXX_COMPILER": args.cxx_compiler,
    "CMAKE_C_COMPILER": args.c_compiler,
    "WITH_MPI": "On" if args.with_mpi is True else "Off",
    "PETSC_DIR":  args.petsc_dir,
    "PETSC_ARCH": args.petsc_arch,
    "SLEPC_DIR": args.slepc_dir,
    "METIS_DIR": args.metis_dir,
}

if args.floating_point_spec == FloatingPointPrecision.double:
    cmake_args["SFEM_USE_DOUBLE_PRECISION"] = "True"
else:
    cmake_args["SFEM_USE_SINGLE_PRECISION"] = "True"

run_cmake = f"cmake -S . -B {args.build_dir} "
for k, v in cmake_args.items():
    if v is not None:
        run_cmake += f"-D{k}={v} "
run_shell_cmd(run_cmake)
# ==============================================================================
# Build
run_make = f"make -j {n_jobs_make()}"
run_shell_cmd(run_make, cwd=args.build_dir)
# ==============================================================================
# Installation
run_make_install = f"cd {args.build_dir}; make install"
run_shell_cmd(run_make_install, cwd=args.build_dir)
# ==============================================================================
# Install the __init__.py file for pysfem
# Create the stub .pyi files for the Python bindings
if args.with_pysfem:
    # __init__.py
    pysfem_src_dir = os.path.join(os.getcwd(), "pysfem")
    pysfem_install_dir = os.path.join(args.install_dir, "lib", "pysfem")
    init_filepath = os.path.join(pysfem_src_dir, "__init__.py")
    install_init_file = f"cp {init_filepath} {pysfem_install_dir}"
    run_shell_cmd(install_init_file)

    # stubs
    create_stubs = f"cd {pysfem_install_dir}; python3 -m nanobind.stubgen -m pysfem --recursive"
    run_shell_cmd(create_stubs)
    move_root_pyi = f"cd {pysfem_install_dir}; mv pysfem.pyi __init__.pyi"
    run_shell_cmd(move_root_pyi)
# ==============================================================================
# Create installation summary file
# Lines to append to .bashrc
bashrc_app = "# SFEM\n"
bashrc_app += f"export SFEM_INSTALL_DIR={args.install_dir}\n"
bashrc_app += f"export SFEM_THIRD_PARTY_DIR={args.third_party_dir}\n"
if args.with_pysfem:
    bashrc_app += "export PYTHONPATH=${PYTHONPATH}:${SFEM_INSTALL_DIR}/lib\n"

if args.with_apps:
    apps_dir = os.path.join(args.install_dir, "bin")
    bashrc_app += f"export PATH=$PATH:{apps_dir}\n"

if args.download_openmpi:
    bashrc_app += "# OpenMPI\n"
    bashrc_app += "export MPICC=${SFEM_THIRD_PARTY_DIR}/openmpi/bin/mpicc\n"
    bashrc_app += "export MPICXX=${SFEM_THIRD_PARTY_DIR}/openmpi/bin/mpic++\n"
    bashrc_app += "export MPIEXEC=${SFEM_THIRD_PARTY_DIR}/openmpi/bin/mpiexec\n"

if args.download_petsc:
    bashrc_app += "# PETSc\n"
    bashrc_app += "export PETSC_DIR=${SFEM_THIRD_PARTY_DIR}/petsc\n"
    bashrc_app += f"export PETSC_ARCH={args.petsc_arch}\n"
    bashrc_app += "export PYTHONPATH=${PYTHONPATH}:${PETSC_DIR}/${PETSC_ARCH}/lib\n"

if args.download_slepc:
    bashrc_app += "# SLEPc\n"
    bashrc_app += "export SLEPC_DIR=${SFEM_THIRD_PARTY_DIR}/slepc\n"
    bashrc_app += "export PYTHONPATH=${PYTHONPATH}:${SLEPC_DIR}/${PETSC_ARCH}/lib\n"

sep_line = "#=================================================#\n"
with open("summary.log", "w") as file:
    lines = sep_line
    lines += "SFEM INSTALLATION SUMMARY\n"

    # SFEM
    lines += sep_line
    lines += f"Installation directory: {args.install_dir}\n"
    lines += f"Third-party dependencies directory: {args.third_party_dir}\n"
    lines += f"pysfem enabled: {args.with_pysfem}\n"
    if args.with_pysfem:
        lines += f"pysfem directory: {get_pysfem_dir(args.install_dir)}\n"
    lines += f"Applications enabled: {args.with_apps}\n"
    if args.with_apps:
        lines += f"Applications directory {get_applications_dir(args.install_dir)}\n"
    lines += f"MPI enabled: {args.with_mpi}\n"
    lines += f"C compiler: {args.c_compiler}\n"
    lines += f"C++ compiler: {args.cxx_compiler}\n"

    # OpenMPI
    if args.download_openmpi:
        lines += sep_line
        lines += "Downloaded OpenMPI\n"
        _, version = get_openmpi_url_and_version()
        lines += f"OpenMPI directory: {openmpi_dir}\n"
        lines += f"OpenMPI version: {version}\n"

    # PETSc
    lines += sep_line
    if args.download_petsc:
        lines += "Downloaded PETSc\n"
    lines += f"PETSc directory: {args.petsc_dir}\n"
    lines += f"PETSc arch: {args.petsc_arch}\n"
    if args.download_petsc and args.with_pysfem:
        lines += f"petsc4py directory: {get_petsc4py_dir(args.petsc_dir, args.petsc_arch)}\n"

    # SLEPc
    lines += sep_line
    if args.download_slepc:
        lines += "Downloaded SLEPc\n"
    lines += f"SLEPc directory: {args.slepc_dir}\n"
    if args.download_slepc and args.with_pysfem:
        lines += f"slepc4py directory: {get_slepc4py_dir(args.slepc_dir, args.petsc_arch)}\n"

    # METIS
    lines += sep_line
    if args.download_metis:
        lines += "Downloaded METIS\n"
    lines += f"METIS directory: {args.metis_dir}\n"

    lines += sep_line
    lines += "INSTALLATION COMPLETE\n"

    # .bashrc
    lines += sep_line
    lines += "Add these lines to the end of your .bashrc (or equivalent) file:\n"
    lines += bashrc_app

    file.writelines(lines)
    print(lines, end='')
