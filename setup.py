import os
from setuptools import setup, Extension
from pybind11.setup_helpers import Pybind11Extension

# Modify these paths to match your setup
path = "/Users/annahart/Documents"
boostdir = "boost_1_74_0"
lemondir = "lemon"

lemon_include_path = f'{path}/{lemondir}/include'
lemon_library_path = f'{path}/{lemondir}/lib'

boost_include_path = f'{path}/{boostdir}/include'
boost_library_path = f'{path}/{boostdir}/lib'
cpp_std_flag = '-std=c++11'

# Define the extension module
ext_modules = [
    Pybind11Extension(
        "clonelib",
        ["src/pybind/gen_trees_pybind.cpp", "src/cnagraph.cpp", "src/cnatree.cpp", "src/basematrix.cpp", "src/basetree.cpp" ,"src/utils.cpp", "src/genotypetree.cpp"],
        include_dirs=[lemon_include_path, boost_include_path, 'src'],
        library_dirs=[lemon_library_path, boost_library_path, 'src'], 
        extra_compile_args=[cpp_std_flag] #,'-g','-O0', '-DDEBUG'],
        # libraries=["boost_1_74_0", "lemon"],  # Optional: Specify additional libraries if needed
    ),
]

# Check if the Lemon Graph Library is installed
if not os.path.exists(lemon_include_path):
    print("Error: Lemon Graph Library not found. Please provide the correct include path.")
    ext_modules = []  # Disable building the extension if Lemon is not found

# Use setuptools to build the extension module
setup(
    name="clonelib",
    ext_modules=ext_modules,
    install_requires=[
        "setuptools",
        "pybind11",  # Make sure pybind11 is installed
    ],
)
