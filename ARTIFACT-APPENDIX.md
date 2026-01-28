# Artifact Appendix

Paper title: **Secure Change-Point Detection for Time Series under Homomorphic Encryption**

Requested Badge(s):
  - **Available**
  - **Functional**
  - **Reproduced**

## Description

This artifact relates to the following paper:

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

This artifact contains the source code of our algorithm for secure change-point detection of time series under the CKKS encryption scheme. The implementation is built on top of the OpenFHE library, and it provides methods to detect shifts in the mean, variance, and frequency of the input data.

This artifact also comes with the three real-world datasets and the code to generate the synthetic data used in the experimental evaluation of our paper.

Additionally, it provides the code for the local differential privacy (DP) solution against which we compare our own solution.

### Security/Privacy Issues and Ethical Concerns

We do not identify any security or privacy risks associated with installing or using this artifact.

The data used to assess our solution is either synthetically generated or taken from public repositories.

## Basic Requirements

### Hardware Requirements

Minimal hardware requirements: there are no special hardware requirements, our artifact can run on a laptop.

Our experiments were performed on a machine with:
- CPU: AMD EPYC 7763 64-Core Processor
- RAM: 512 GB (only 64 GB were used)

### Software Requirements

You might want to use `git` to clone the repository: [git 2.39.5 (any version should work)](https://git-scm.com/downloads).

If building through **Docker Image**:
- [Docker Engine 29.1.3](https://docs.docker.com/engine/install/) (should work with any recent version)

If building through **manual installation**:
- OS: Linux Ubuntu 22.04 or 24.04
- C++ compiler: GCC v10 or Clang 11 (other versions might work as well)
- CMake v3.22.1 (other versions might work as well)
- Python 3.10.6 (other versions might work as well)
  - Python packages: numpy, matplotlib, scipy
- [OpenFHE v1.1.2](https://github.com/openfheorg/openfhe-development/releases/tag/v1.1.2) (we include this version in our repository)

### Estimated Time and Storage Consumption

- Human Time: 2 minutes
- Computing Time: 2.5 hours (varies based on hardware)
- Disk Space: 0.5 GB
- Memory: up to 64 GB RAM

## Environment

### Accessibility

The code is publicly available on [GitHub](https://github.com/FedericoMazzone/secure-change-point-detection).

### Set up the environment

Download the repository:
```bash
git clone https://github.com/FedericoMazzone/secure-change-point-detection
cd secure-change-point-detection
```

<details>

<summary>Install via Docker</summary>

Build the Docker image:

```bash
docker build -t secure-cpd .
```

Run the Docker container:

```bash
docker run -it secure-cpd
```

</details>

<details>
  <summary>Install manually</summary>

Install prerequisites:
```bash
sudo apt-get install build-essential cmake python3
pip3 install -r requirements.txt
```

Note: on some systems this may fail with "This environment is externally managed".
In this case, you can use a virtual environment (see below):
```bash
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
```

Build the project:
```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

To clean the build (and remove the local OpenFHE installation):
```bash
rm -rf build
```
</details>

### Testing the Environment

To do a basic test run of our solution:

```bash
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

We expect the error of the local DP solution to increase drastically for low privacy budgets, in particular we expect the relative error to be much higher than the baseline is.
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
bash ./benchmark-synthetic.sh
```

The logs of the individual experiments are stored in the `logs` folder, while the resulting runtime and memory results are extracted and summarized in the file `logs/benchmark-synthetic.out`, one row per experiment. These can be directly compared to the results reported in the paper (Table 1). Due to the high reliance on CPU power and parallelization, we can only promise the runtime results to be quantitatively in the same order of magnitude of the expected results (while the memory results should be reasonably close to the expected ones).

#### Experiment 2: Real-world accuracy
- Time: 1 human-minute + 10 compute-minutes
- Memory: 20 GB RAM

This experiment reproduces
[Main Result 2: Real-world accuracy](#main-result-2-real-world-accuracy).
The following script will run the simulation automatically over the three real-world datasets:

```bash
bash ./benchmark-realworld.sh
```

The logs of the individual experiments are stored in the `logs` folder, while the resulting accuracy results in terms of relative error are stored in the file `logs/benchmark-realworld.out`, one row per experiment, together with the plaintext baseline. These can be directly compared to the results reported in the paper, and should not quantitatively vary by more than 0.01 percentage point from expected results.


#### Experiment 3: Local DP accuracy
- Time: 1 human-minute + 10 compute-minutes

This experiment reproduces
[Main Result 3: Local DP accuracy](#main-result-3-local-dp-accuracy).
The following Python script will run all the experiments a given amount of times and plot the average relative error:

```bash
python3 src/local-dp.py
```

The resulting plot will be stored in the PNG file `logs/local-dp/benchmark-localDP.png`, while the numerical errors will be saved in the CSV file `logs/local-dp/benchmark-localDP.csv`. This benchmark can be directly compared to the one in Figure 5 and should show a similar trend.


## Limitations

All the main results related to our solution can be reproduced through this artifact.

## Notes on Reusability

Our solution has been implemented as the combination of two components: the summarization component (three are available: mean, variance, frequency) and the CUSUM computation component. These components can be modified or optimized independently. The CUSUM computation can be used for independent projects or works.

The precision parameters, approximation degrees, and other HE-related parameters can be easily modified in the code to adapt the solution to different scenarios or requirements.

Since the maximum supported block size is 65,536, for very large time series (e.g., n > 2^32) the default sqrt(n) block size heuristic may exceed this limit and should be adjusted manually.