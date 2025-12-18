# Secure Change-Point Detection

This repository provides a library for performing change-point detection for time series, using the CKKS (Cheon-Kim-Kim-Song) homomorphic encryption scheme.
Our code is built on top of the OpenFHE library.

This repository accompanies the paper *Secure Change-Point Detection for Time Series under Homomorphic Encryption* in publication at PETS 2026.
If you use this work, please cite:
```bibtex
@inproceedings{mazzone2026secure,
  author    = {Mazzone, Federico and Micali, Giorgio and Pronesti, Massimiliano},
  title     = {Secure Change-Point Detection for Time Series under Homomorphic Encryption},
  booktitle = {26th Privacy Enhancing Technologies Symposium (PETS 2026)},
  year      = {2026},
  month     = {jul},
  address   = {Calgary, Canada},
  note      = {To appear}
}
```

## Installation

### Via Docker Image

```bash
docker build -t secure-cpd .
docker run -it secure-cpd
```

### Via Manual Installation

Install prerequisites:
```bash
sudo apt-get install build-essential
sudo apt-get install cmake
sudo apt-get install python3
pip3 install numpy
pip3 install matplotlib
pip3 install scipy
```

Install OpenFHE:
```bash
cd openfhe-development-1.1.2
mkdir build
cd build
cmake ..
make -j
sudo make install
cd ../..
echo 'export LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH' >> ~/.bashrc
source ~/.bashrc
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

Install our library:
```bash
mkdir build
cd build
cmake ..
make -j
cd ..
```



## Usage

After successful compilation, the `cpd` executable will be available in the `build` directory.
To run it, use the following syntax:

```bash
./build/cpd -f|--file <filepath> [-l|--labeled] [-v|--verbose] [-h|--help] [-b|--block-size <size>] [-t|--type <mean|variance|frequency>]
```

Some datasets can be found in the `data` folder. The synthetic ones can be generated with the python scripts in `src/data_utils`.

For example:

```bash
./build/cpd -f data/meditation.csv -l -v -t frequency
```

## Data utilities

The `src/data_utils/` folder contains small scripts to generate synthetic time series with a known change-point and visualize a time series together with detected change-points.

### Output format (generators)

Both generators save a **single CSV row** with:

* `x1,x2,...,xn,t`
* where `t` is the **ground-truth change-point index** (integer, in the original time scale).
* values are rounded to **3 decimals**.
* If `--plot` is enabled, a `.png` with the same basename is also saved next to the CSV.

### Synthetic data generation: distribution change (mean/variance)

Generates i.i.d. samples with a single change in distribution at time `t = floor(tau * n)`. The time series is generated from `dist1(mu1, sigma1)` before the change-point, and from `dist2(mu2, sigma2)` after it.

**Command**

```bash
python3 src/data_utils/ts_gen_distr_change.py \
  --folder <output_folder> \
  --dist1 <normal|laplace|t|uniform|chi2> \
  --dist2 <normal|laplace|t|uniform|chi2> \
  --n <num_points> \
  --mu1 <mean_before> --mu2 <mean_after> \
  --sigma1 <std_before> --sigma2 <std_after> \
  --tau <change_point_fraction> \
  [--plot]
```

**Defaults**

* `--folder data/generated/distr_change`
* `--dist1 normal`, `--dist2 normal`
* `--n 10000`, `--mu1 0.0`, `--mu2 0.0`, `--sigma1 1.0`, `--sigma2 1.0`
* `--tau 0.5`
* `--plot` disabled (flag)

**Example**

```bash
python3 src/data_utils/ts_gen_distr_change.py \
  --dist1 normal --dist2 normal \
  --n 10000 --mu1 0.0 --mu2 1.0 \
  --sigma1 1.0 --sigma2 2.0 \
  --tau 0.5 --plot
```

### Synthetic data generation: frequency change (AR(1) series)

Generates an AR(1) time series with a change in the AR coefficient (from `phi1` to `phi2`) at `t = floor(tau * n)`. A burn-in phase can be used to stabilize the process.

**Command**

```bash
python3 src/data_utils/ts_gen_freq_change.py \
  --folder <output_folder> \
  --distribution <normal|laplace|t> \
  --n <num_points> \
  --phi1 <ar_coeff_before> \
  --phi2 <ar_coeff_after> \
  --tau <change_point_fraction> \
  --burn <burn_in_samples> \
  [--plot]
```

**Defaults**

* `--folder data/generated/frequency-change`
* `--distribution normal`
* `--n 10000`
* `--phi1 0.3`, `--phi2 0.7`
* `--tau 0.5`
* `--burn 500`
* `--plot` disabled (flag)

**Example**

```bash
python3 src/data_utils/ts_gen_freq_change.py \
  --distribution normal \
  --n 10000 --phi1 0.3 --phi2 0.7 \
  --tau 0.5 --burn 500 --plot
```

### Visualization and baseline CPD (plotter)

`plotter.py` loads a CSV in the `x1,...,xn,t` format, runs the CPD algorithm in plaintext, and saves a 3-panel PNG:

1. original time series + vertical lines for ground-truth and computed change-point
2. block statistics `q`
3. CUSUM statistic `|cumsum(q - mean(q))|`

**Command**

```bash
python3 src/data_utils/plotter.py \
  <filepath> \
  --change_type <mean|variance|frequency> \
  [--m <block_size>]
```

**Parameters**

* `<filepath>`: path to a CSV containing `x1,...,xn,t` (required positional argument)
* `--change_type`: statistic used per block (`mean`, `variance`, or `frequency`)
  default: `frequency`
* `--m`: block size. If omitted, it defaults to `floor(sqrt(n))`.

**Example**

```bash
python3 src/data_utils/plotter.py data/meditation.csv --change_type frequency
```

**Output**

* Saves an image next to the CSV, replacing `.csv` with `.png` (e.g., `data/meditation.png`).
