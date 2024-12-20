#include <cstdint>

//
#include <algorithm>
#include <chrono>
#include <iostream>
#include <numeric>
#include <random>
// GSRAP
#include "gsrap/macros.h"
#include "gsrap/pnp_solver.h"

// common
#include "common/camera.h"
#include "common/type.h"
#include "common/utils.h"

GSRAP_IGNORE_STRICT_WARNING_PUSH

#include "Eigen/Core"
#include "Eigen/Dense"

GSRAP_IGNORE_STRICT_WARNING_POP

using namespace gsrap::example;  // NOLINT

namespace gsrap::benchmark {

using Pair     = std::pair<uint32_t, uint32_t>;
using Iterator = std::vector<Pair>::const_iterator;

static Camera CreateCameras() {
  Camera camera = {};
  {
    // Set optional axis and up at the world coordinates.
    // Variable optional_axis_w is equal to the z axis of the coordinates of
    // camera.
    const Vec3d camera_optical_axis_w = Vec3d(-1.0, 1.0, 1.0).normalized();
    const Vec3d camera_up_w           = Vec3d(0.0, -1.0, 0.0);

    // Compute x and y axes of the coordinates of camera.
    const Vec3d camera_xw =
        camera_up_w.cross(camera_optical_axis_w).normalized();
    const Vec3d camera_yw = camera_optical_axis_w.cross(camera_xw).normalized();

    // Create world to camera rotation matrix.
    camera.rotation_cw.row(0) = camera_xw;
    camera.rotation_cw.row(1) = camera_yw;
    camera.rotation_cw.row(2) = camera_optical_axis_w;

    // Set camera to world translation.
    // The distance from camera to camera1 is 1.
    const Vec3d camera_translation_wc =
        Vec3d(1.0 / sqrt(2.0), -1.0 / sqrt(2.0), 0.0);

    // Compute world to camera translation.
    camera.transltaion_cw =
        -camera.rotation_cw.transpose() * camera_translation_wc;
  }

  return camera;
}
static void BenchPnpSolver() {
  const uint32_t num_points = 4;  // Number of points

  const Vec3d center(0.5 / sqrt(2.0), -0.5 / sqrt(2.0),
                     0.5 / sqrt(2.0));    // Center of spherical point cloud.
  const double radius = 0.5 / sqrt(2.0);  // radius of spherical point cloud

  std::mt19937 engine(123456789);

  const size_t num_itr = 100000;

  // To increase the frequency of the CPU
  for (size_t bench_itr = 0; bench_itr < num_itr; ++bench_itr) {
    eigen_vector<Vec3d> points =
        CreatesphericalPointCloud(center, radius, num_points, &engine);

    const auto camera = CreateCameras();

    const auto bearings = ToBearingVectors(points, camera);

    [[maybe_unused]] const auto result =
        gsrap::SolvePnpProblemWithoutRansac(points, points);
  }

  std::vector<int64_t> elapsed_times;
  for (size_t bench_itr = 0; bench_itr < num_itr; ++bench_itr) {
    eigen_vector<Vec3d> points =
        CreatesphericalPointCloud(center, radius, num_points, &engine);

    const auto camera = CreateCameras();

    const auto bearings = ToBearingVectors(points, camera);

    const auto start = std::chrono::high_resolution_clock::now();
    [[maybe_unused]] const auto result = gsrap::SolvePnpProblemWithoutRansac(
        points, points);  //  TODO(kyawakyaw): Check accuracy
    const auto end = std::chrono::high_resolution_clock::now();

    const std::chrono::nanoseconds elapsed_time =
        std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);

    elapsed_times.emplace_back(elapsed_time.count());
  }

  std::sort(elapsed_times.begin(), elapsed_times.end());

  const double mean =
      double(std::accumulate(elapsed_times.begin(), elapsed_times.end(), 0)) /
      double(elapsed_times.size());

  const int64_t median = elapsed_times.at(elapsed_times.size() / 2);

  const int64_t maxi = elapsed_times.back();

  std::cout << "mean: " << mean << " ns, median: " << median
            << " ns, max: " << maxi << " ns" << std::endl;
}

}  // namespace gsrap::benchmark

int main(void) {
  gsrap::benchmark::BenchPnpSolver();
  return 0;
}
