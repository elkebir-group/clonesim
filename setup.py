import os
from setuptools import setup, Extension
from pybind11.setup_helpers import Pybind11Extension

# Modify these paths to match your setup
lemon_include_path = '/Users/melkebir/lemon/include'
lemon_library_path = '/Users/melkebir/lemon/lib'
boost_include_path = '/opt/homebrew/include'
boost_library_path = '/opt/homebrew/lib'
cpp_std_flag = '-std=c++11'

# Define the extension module
ext_modules = [
    Pybind11Extension(
        "cnatrees",
        ["src/gen_trees_pybind.cpp", "src/gencnatrees.cpp", "src/basematrix.cpp", "src/basetree.cpp" ,"src/utils.cpp", "src/genotypetree.cpp", "src/cnatree.cpp"],
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
    name="cnatrees",
    ext_modules=ext_modules,
    install_requires=[
        "setuptools",
        "pybind11",  # Make sure pybind11 is installed
    ],
)
