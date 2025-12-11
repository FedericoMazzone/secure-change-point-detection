# Artifact Appendix

Paper title: **Secure Change-Point Detection for Time Series under Homomorphic Encryption**

Requested Badge(s):
  - **Available**
  - **Functional**
  - **Reproduced**

## Description

This artifact relates to the following paper:

```
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

This artifact contains the source code of our algorithm for secure change-point detection of time series under the CKKS encryption scheme. The implementation is built on top of the OpenFHE library, and it provides methods to detect shifts in the mean, variance, and frequency of the input data.

This artifact also comes with the three real-world datasets and the code to generate the synthetic data used in the experimental evaluation of our paper.

Additionally, it provides the code for the local differential privacy (DP) solution against which we compare our own solution.

### Security/Privacy Issues and Ethical Concerns

We do not identify any security or privacy risks associated with installing or using this artifact.

The data used to assess our solution is either synthetically generated or taken from public repositories.

## Basic Requirements

### Hardware Requirements

Minimal hardware requirements: there are no special hardware requirements, our artifact can run on a laptop.

The specifications of the hardware on which the experiments reported in the paper were performed are: AMD EPYC 7763 64-Core Processor running at 2.45GHz, 512 GB RAM.
Although, according to our measurements, 64 GB of RAM should be enough to run our experiments.

### Software Requirements

We have tested our artifact on Linux Ubuntu 24.04 (it should work on 22.04 too). It should theoretically compile and function on Microsoft Windows, but we do not guarantee it.

Building and running our implementation requires a C++ compiler, CMake, and Python. We recommend using GCC v10 or Clang 11 on Linux for best compatibility. We used cmake version 3.28.3, and Python 3.12.3 (with numpy, matplotlib, and scipy).

Our implementation relies on the OpenFHE library, available at https://github.com/openfheorg/openfhe-development. We specifically use version 1.1.2, which is included in
our repository to avoid compatibility issues.

### Estimated Time and Storage Consumption

The estimated overall human time to run the artifact is 2 minutes. While the actual computing to reproduce all the experiments (one time) requires around 2 hours, and uses at most 64 GB of RAM and 0.5 GB of disk space.

## Environment

### Accessibility

This artifact can be accessed through GitHub at https://github.com/FedericoMazzone/secure-change-point-detection.

### Set up the environment

Download the repository:
```
git clone https://github.com/FedericoMazzone/secure-change-point-detection
```

Install prerequisites:
```
sudo apt-get install build-essential
sudo apt-get install cmake
sudo apt-get install python3
pip3 install numpy
pip3 install matplotlib
pip3 install scipy
pip3 install tqdm
```

Install OpenFHE:
```
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

Install our library:
```
mkdir build
cd build
cmake ..
make -j
cd ..
```

If you do not have sudo access in your machine, specify a different installation path for OpenFHE. See more about this in the README file.

### Testing the Environment

To do a basic test run of our solution:

```
./build/cpd -f data/meditation.csv -l -v -t frequency
```

which should start the detection of a frequency-change point for the given dataset. You can expect as output the computation logs together with the runtime measurements.

## Artifact Evaluation

### Main Results and Claims

#### Main Result 1: Runtime/memory scalability

By varying the number of points in the time-series, we expect to see the runtime (Server - Total (s)) and memory usage (Server - Memory (GB)) to roughly match what is reported in Table 1.
This claim is reproducible by executing our
[Experiment 1](#experiment-1-runtimememory-scalability). In this experiment, we test our change-point detection algorithm 1'000, 10'000, 100'000, and 1'000'000 points on synthetic time-series for changes in mean, variance, and frequency.
The runtime and memory consumption of frequency-change should be higher than variance-change, which should be slightly higher than mean-change.

#### Main Result 2: Real-world accuracy

We expect the accuracy of our solution to be very close to the plaintext baseline, in terms of relative error.
This claim is reproducible by executing our
[Experiment 2](#experiment-2-real-world-accuracy). In this experiment, we test our change-point detection algorithm on three real-world datasets for frequency changes.
The accuracy of our approach and of the plaintext baseline should match the one reported in Table 4.

#### Main Result 3: Local DP accuracy

We expect the accuracy of the local DP solution to increase drastically for low privacy budgets, in particular we expect the relative error to be much higher than the baseline is.
This claim is reproducible by executing our
[Experiment 3](#experiment-3-local-dp-accuracy). In this experiment, we test the local DP solution on increasing privacy budget (epsilon) for time-series with 1'000, 10'000, 100'000, and 1'000'000 points. The relative error is averaged over multiple runs.
The relative error trend should reflect the one given in Figure 5 of the paper.


### Experiments

#### Experiment 1: Runtime/memory scalability
- Time: 1 human-minute + 2 compute-hours.
- Memory: 64 GB RAM

This experiment reproduces
[Main Result 1: Runtime/memory scalability](#main-result-1-runtimememory-scalability).
The following script will generate all the synthetic data and run the simulation automatically over it:

```bash
sh ./benchmark-synthetic.sh
```

The resulting runtime and memory results are stored in the CSV file `benchmark-synthetic.out`, one row per experiment. These can be directly compared to the results reported in the paper. Due to the high reliance on CPU power and parallelization, we can only promise the results to be quantitatively in the same order of magnitude of the expected results.

#### Experiment 2: Real-world accuracy
- Time: 1 human-minute + 10 compute-minutes
- Memory: 20 GB RAM

This experiment reproduces
[Main Result 2: Real-world accuracy](#main-result-2-real-world-accuracy).
The following script will run the simulation automatically over the three real-world datasets:

```bash
sh ./benchmark-realworld.sh
```

The resulting accuracy results in terms of relative error are stored in the CSV file `benchmark-realworld.out`, one row per experiment, together with the plaintext baseline. These can be directly compared to the results reported in the paper, and should not quantitatively vary by more than 0.01 percentage point from expected results.


#### Experiment 3: Local DP accuracy
- Time: 1 human-minute + 1 compute-hour

This experiment reproduces
[Main Result 3: Local DP accuracy](#main-result-3-local-dp-accuracy).
The following Python script will run all the experiments a given amount of times and plot the average relative error:

```bash
python3 src/local-dp.py
```

The resulting plot will be stored in the PNG file `benchmark-localDP.png`. This plot can be directly compared to the one in FIgure 5 and should show a similar trend.


## Limitations

All the main results related to our solution can be reproduced through the artifact.

## Notes on Reusability

Our solution has been implemented as the combination of two components: the summarization component (three are available: mean, variance, frequency) and the CUSUM computation component. These components can be modified or optimized indipendently. The CUSUM computation can be used for independent project or works.