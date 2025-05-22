# Privacy-Preserving Vertical K-Means Clustering

This repository provides a library for performing change-point detection for time series, using the CKKS (Cheon-Kim-Kim-Song) homomorphic encryption scheme.
Our code is built on top of the OpenFHE library.


## Installation (Linux)

### OpenFHE Library

To ensure compatibility, we recommend using [version **1.1.2** of the OpenFHE library](https://github.com/openfheorg/openfhe-development/tree/eeb2ca1b88333dccdd96bef9720456710a3facbb).
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

   After successful compilation, the `clustering` executable will be available in the `build` directory.

---
