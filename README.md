# Secure Change-Point Detection

This repository provides a library for performing change-point detection for time series, using the CKKS (Cheon-Kim-Kim-Song) homomorphic encryption scheme.
Our code is built on top of the OpenFHE library.


## Installation (Linux Ubuntu)


### Prerequisites

Install compiler and cmake if needed.

   ```bash
   sudo apt-get install build-essential
   sudo apt-get install cmake
   ```


### OpenFHE Library

To ensure compatibility, we recommend using version **1.1.2** of the OpenFHE library. This version is included within this repository. Follow the steps below to install it:

1. **Navigate to the OpenFHE Directory**

   ```bash
   cd openfhe-development-1.1.2
   ```

2. **Create a Build Directory**

   ```bash
   mkdir build
   ```

3. **Navigate to the Build Directory**

   ```bash
   cd build
   ```

4. **Generate Build Files with CMake**

   ```bash
   cmake ..
   ```

5. **Compile the Library**

   ```bash
   make -j
   ```

6. **Install the Library**

   ```bash
   sudo make install
   ```

7. **Return to the Root Directory**

   ```bash
   cd ../..
   ```

If you do not have sudo access in your machine, specify a different installation path by

   ```bash
   cmake -DCMAKE_INSTALL_PREFIX=~/openfhe ..
   make -j
   make install
   cd ../..
   echo 'export LD_LIBRARY_PATH=~/openfhe/lib:$LD_LIBRARY_PATH' >> ~/.bashrc
   echo 'export LIBRARY_PATH=~/openfhe/lib:$LIBRARY_PATH' >> ~/.bashrc
   echo 'export CMAKE_PREFIX_PATH=~/openfhe:$CMAKE_PREFIX_PATH' >> ~/.bashrc
   echo 'export CPATH=~/openfhe/include:$CPATH' >> ~/.bashrc
   source ~/.bashrc
   ```

For additional support or troubleshooting, refer to the official [OpenFHE Documentation](https://openfhe-development.readthedocs.io/en/latest/).

### Our Library

Follow these steps to compile our library and build the demo and benchmarking executables:

1. **Create a Build Directory**

   Create a directory named `build` in the root of the repository:

   ```bash
   mkdir build
   ```

2. **Navigate to the Build Directory**

   Change into the `build` directory:

   ```bash
   cd build
   ```

3. **Generate Build Files with CMake**

   Use CMake to generate build files based on the configuration specified in `CMakeLists.txt`:

   ```bash
   cmake ..
   ```

4. **Compile the Executables**

   Build the source files to create the executables:

   ```bash
   make -j
   ```

5. **Return to the Root Directory**

   ```bash
   cd ..
   ```

6. **Run the Executables**

   After successful compilation, the `cpd` executable will be available in the `build` directory.
   To run it, use the following syntax:

   ```bash
   ./build/cpd -f|--file <filepath> [-l|--labeled] [-v|--verbose] [-h|--help] [-b|--block-size <size>] [-t|--type <mean|variance|frequency>]
   ```

   Some datasets can be found in the `data` folder. The synthetic ones can be generated with the python scripts in `src/data-utils`.

   For example:

   ```bash
   ./build/cpd -f data/meditation.csv -l -v -t frequency
   ```

---


## Known Issues

**Missing Shared Library**

If you encounter the following error
```error while loading shared libraries: libOPENFHEbinfhe.so.1: cannot open shared object file: No such file or directory```, it means the OpenFHE shared libraries are not in the dynamic linker’s search path.

SOLUTION
1. Temporarily set `LD_LIBRARY_PATH` (valid for the current session only):
```bash
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib
```
2. Make it permanent (applies to future sessions):
```bash
echo 'export LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH' >> ~/.bashrc
source ~/.bashrc
```
If OpenFHE is installed in a different location, replace /usr/local/lib with the correct path.
