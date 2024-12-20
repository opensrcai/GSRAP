#include <iostream>
#include <random>
#include <thread>

// GSRAP
#include "gsrap/essential_solver.h"
#include "gsrap/macros.h"

// common
#include "common/camera.h"
#include "common/utils.h"

GSRAP_IGNORE_STRICT_WARNING_PUSH

#include "Eigen/Core"
#include "Eigen/Dense"

GSRAP_IGNORE_STRICT_WARNING_POP

// NOLINTNEXTLINE
using namespace gsrap;

// NOLINTNEXTLINE
using namespace gsrap::example;

namespace gsrap::example {

using Pair     = std::pair<uint32_t, uint32_t>;
using Iterator = std::vector<Pair>::const_iterator;

// This function creates two cameras.
static std::pair<Camera, Camera> CreateTwoCameras() {
  // The coordinates of cameara1 are equal to the world coordinates.
  Camera camera1         = {};
  camera1.rotation_cw    = Mat3d::Identity();
  camera1.transltaion_cw = Vec3d::Zero();

  Camera camera2 = {};
  {
    // Set optional axis and up at the world coordinates.
    // Variable optional_axis_w is equal to the z axis of the coordinates of
    // camera2.
    const Vec3d camera2_optical_axis_w = Vec3d(-1.0, 1.0, 1.0).normalized();
    const Vec3d camera2_up_w           = Vec3d(0.0, -1.0, 0.0);

    // Compute x and y axes of the coordinates of camera2.
    const Vec3d camera2_xg =
        camera2_up_w.cross(camera2_optical_axis_w).normalized();
    const Vec3d camera2_yg =
        camera2_optical_axis_w.cross(camera2_xg).normalized();

    // Create world to camera rotation matrix.
    camera2.rotation_cw.row(0) = camera2_xg;
    camera2.rotation_cw.row(1) = camera2_yg;
    camera2.rotation_cw.row(2) = camera2_optical_axis_w;

    // Set camera to world translation.
    // The distance from camera2 to camera1 is 1.
    const Vec3d camera2_translation_wc =
        Vec3d(1.0 / sqrt(2.0), -1.0 / sqrt(2.0), 0.0);

    // Compute world to camera translation.
    camera2.transltaion_cw =
        -camera2.rotation_cw.transpose() * camera2_translation_wc;
  }

  return std::make_pair(camera1, camera2);
}

static void RunEssentialSolver(void) {
  const uint32_t num_points = 500;  // Number of points

  const Vec3d center(0.5 / sqrt(2.0), -0.5 / sqrt(2.0),
                     0.5 / sqrt(2.0));    // Center of spherical point cloud.
  const double radius = 0.5 / sqrt(2.0);  // radius of spherical point cloud

  std::mt19937 engine(123456789);
  eigen_vector<Vec3d> points =
      CreatesphericalPointCloud(center, radius, num_points, &engine);

  const auto [camera1, camera2] = CreateTwoCameras();

  const auto [b1s, b2s] = ToBearingVectors(points, camera1, camera2);

  const std::vector<Pair> matches =
      CreateMatches(num_points, /* inlier_ratio */ 1.0, &engine);

  EssentialSolverPolicy essential_solver_policy = {};

  essential_solver_policy.flags =
      ESSENTIAL_SOLVER_POLICY_FLAGS_CHECK_SINGULAR_VALUE;
  essential_solver_policy.ransac_policy.flags =
      RANSAC_POLICY_FLAGS_EARLY_STOP |
      RANSAC_POLICY_FLAGS_USE_PROBABILITY_WITHOUT_DUPLICATION_SAMPLE;
  essential_solver_policy.ransac_policy.probability = 0.999999;
  essential_solver_policy.ransac_policy.num_threads =
      std::thread::hardware_concurrency();
  essential_solver_policy.ransac_policy.num_ransac_itr_lower_limit = 30;
  essential_solver_policy.ransac_policy.num_ransac_itr_upper_limit = 100;

  const std::pair<std::optional<EssentialSolverResult<Iterator>>, RansacReport>
      result = ComputeEssentialMatrix(essential_solver_policy, b1s, b2s,
                                      matches.cbegin(), matches.cend());

  if (result.first) {
    const Mat3d gt_r = camera2.rotation_cw * camera1.rotation_cw.transpose();
    const Vec3d gt_t = camera2.transltaion_cw - camera1.transltaion_cw;
    Mat3d gt_t_x;
    gt_t_x << 0, -gt_t[2], gt_t[1], gt_t[2], 0, -gt_t[0], -gt_t[1], gt_t[0], 0;

    std::cout << "E :\n" << result.first->essential << "\n" << std::endl;
    std::cout << "gt E :\n" << gt_t_x * gt_r << "\n" << std::endl;

    std::cout << "R :\n" << result.first->rotation << "\n" << std::endl;
    std::cout << "gt R :\n" << gt_r << "\n" << std::endl;

    std::cout << "t :\n" << result.first->translation << "\n" << std::endl;
    std::cout << "gt t :\n" << gt_t << "\n" << std::endl;
  }
}

}  // namespace gsrap::example

int main(void) {
  RunEssentialSolver();
  return 0;
}
