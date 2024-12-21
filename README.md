# GSRAP
**Geometric Solvers for Reconstruction And Pose estimation**

GSRAP is a library designed to handle problems related to multi-view geometry easily.

Currently, it supports the following problems:
- Relative pose estimation from 2D-2D point correspondences
    - Nisterʼs Five points algorithm
- Absolute pose estimation from 2D-3D point correspondences
    - P3P algorithm, EPnP
- Similarity transformation estimation from 3D-3D point correspondences
    - Umeyama algorithm

It also has the following features:
- RANSAC for outlier rejection
    - A versatile framework for rejecting outliers
- C++17
- Python wrapper
- Tutorial and sample codes

## Tutoral slides
[Slide Link](https://speakerdeck.com/opensourceai/tutorial-of-geometric-solvers-for-reconstruction-and-pose-estimation-gsrap)  
[Japanese Version](https://speakerdeck.com/opensourceai/tutorial-of-geometric-solvers-for-reconstruction-and-pose-estimation-gsrap-b750ec22-e3c1-40f1-bba8-e712e47072c1)  

## Usage
- Relative pose estimation from 2D-2D point correspondences

```cpp
  EssentialSolverPolicy essential_solver_policy     = {};
  essential_solver_policy.ransac_policy.probability = 0.999999;

  essential_solver_policy.flags =
      ESSENTIAL_SOLVER_POLICY_FLAGS_CHECK_SINGULAR_VALUE;
  essential_solver_policy.ransac_policy.flags =
      RANSAC_POLICY_FLAGS_EARLY_STOP |
      RANSAC_POLICY_FLAGS_USE_PROBABILITY_WITHOUT_DUPLICATION_SAMPLE;
  essential_solver_policy.ransac_policy.num_threads =
      std::thread::hardware_concurrency();
  essential_solver_policy.ransac_policy.num_ransac_itr_lower_limit = 30;
  essential_solver_policy.ransac_policy.num_ransac_itr_upper_limit = 100;

  const std::pair<std::optional<EssentialSolverResult<Iterator>>, RansacReport>
      result = ComputeEssentialMatrix(essential_solver_policy, b1s, b2s,
                                      matches.cbegin(), matches.cend());
```

- Absolute pose estimation from 2D-3D point correspondences
```cpp
  PnpSolverPolicy pnp_solver_policy = {};
  pnp_solver_policy.ransac_policy   = {};

  pnp_solver_policy.flags =
      PNP_SOLVER_POLICY_FLAGS_NONE | PNP_SOLVER_POLICY_FLAGS_REFINE;
  pnp_solver_policy.ransac_policy.flags =
      RANSAC_POLICY_FLAGS_EARLY_STOP |
      RANSAC_POLICY_FLAGS_USE_PROBABILITY_WITHOUT_DUPLICATION_SAMPLE;
  pnp_solver_policy.ransac_policy.num_threads =
      std::thread::hardware_concurrency();
  pnp_solver_policy.ransac_policy.num_ransac_itr_lower_limit = 10;
  pnp_solver_policy.ransac_policy.num_ransac_itr_upper_limit = 1000000;

  PnpInlierCheckParamsUsingBearingVector pnp_inlier_check_params = {};
  pnp_inlier_check_params.inlier_thr                                = 1e-6;
  pnp_solver_policy.pnp_inlier_check_params = pnp_inlier_check_params;

  const std::pair<std::optional<PnpSolverResult<Iterator>>, RansacReport>
      result = SolvePnpProblem(pnp_solver_policy, bearings, points,
                               matches.begin(), matches.end());
```

- Similarity transformation estimation from 3D-3D point correspondences
```cpp
 Sim3SolverPolicy sim3_solver_policy = {};

  sim3_solver_policy.flags =
      SIM3_SOLVER_POLICY_FLAGS_NONE | SIM3_SOLVER_POLICY_FLAGS_REFINE;
  sim3_solver_policy.ransac_policy = {};
  sim3_solver_policy.ransac_policy.flags =
      RANSAC_POLICY_FLAGS_EARLY_STOP |
      RANSAC_POLICY_FLAGS_USE_PROBABILITY_WITHOUT_DUPLICATION_SAMPLE;
  sim3_solver_policy.ransac_policy.num_threads =
      std::thread::hardware_concurrency();
  sim3_solver_policy.ransac_policy.num_ransac_itr_lower_limit = 10;
  sim3_solver_policy.ransac_policy.num_ransac_itr_upper_limit = 1000000;

  sim3_solver_policy.inlier_thr = 1e-12;

  const std::pair<std::optional<Sim3SolverResult<Iterator>>, RansacReport>
      result = ComputeSim3Transformation(sim3_solver_policy, points1, points2,
                                         matches.begin(), matches.end());
```

## Build

### Build with Vcpkg
To easily build this project, it is convenient to use Vcpkg to install the required libraries.

```bash
git clone --recursive https://github.com/opensrcai/GSRAP.git
cd GSRAP
git clone https://github.com/Microsoft/vcpkg.git
```

If you're using Windows:

```powershell
.\vcpkg\bootstrap-vcpkg.bat
```

For Mac or Linux:
```
./vcpkg/bootstrap-vcpkg.sh
```

Additionally, to build the Python wrapper, you need to install nanobind.

```
pip install nanobind
```

Next, build this project.

### For Windows
Build with Visual Studio.

### For Mac or Linux

```bash
cmake -H. -Bbuild \
          -DCMAKE_BUILD_TYPE=Release \
          -DCMAKE_EXPORT_COMPILE_COMMANDS=ON \
          -DGSRAP_FIND_EIGEN=ON \
          -DGSRAP_FIND_SPDLOG=ON \
          -DGSRAP_FIND_THREADS=ON \
          -DGSRAP_BUILD_TESTS=ON \
          -DGSRAP_BUILD_EXAMPLES=ON \
          -DGSRAP_BUILD_BENCHMARKS=ON \
          -DGSRAP_ENABLE_NANOBIND=ON \
          -DCMAKE_TOOLCHAIN_FILE=./vcpkg/scripts/buildsystems/vcpkg.cmake \
          -DCMAKE_MODULE_PATH=$(dirname $(which python))/../lib/python$(python -c 'import sys; print(".".join(map(str, sys.version_info[:2])))')/site-packages/nanobind/cmake/ \
          -DCMAKE_POSITION_INDEPENDENT_CODE=ON
cd build
make -j16
```

Then, you can run an example program.
```bash
./example_essential_solver
```

## Python binding

You can use this library from Python.  
There are examples of Python code in `example_py` directory.  
Please also check the demo on Google Colab.  

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/opensrcai/GSRAP/blob/main/example_py/notebooks/PnpRelocalize.ipynb)


## LICENSE

Copyright © 2024 Kai Okawa, Mikiya Shibuya.  
Released under the Pre-Open Source Verification License.  

Please see the [LICENSE](LICENSE) file for details.
