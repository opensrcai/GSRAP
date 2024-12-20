#include <iostream>
#include <random>
#include <thread>

// GSRAP
#include "gsrap/macros.h"
#include "gsrap/pnp_solver.h"

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

static void RunPnpSolver(void) {
  const uint32_t num_points = 500;  // Number of points

  const Vec3d center(0.5 / sqrt(2.0), -0.5 / sqrt(2.0),
                     0.5 / sqrt(2.0));    // Center of spherical point cloud.
  const double radius = 0.5 / sqrt(2.0);  // radius of spherical point cloud

  std::mt19937 engine(123456789);
  eigen_vector<Vec3d> points =
      CreatesphericalPointCloud(center, radius, num_points, &engine);

  const auto camera = CreateCameras();

  const auto bearings = ToBearingVectors(points, camera);

  const std::vector<Pair> matches = CreateMatches(num_points, 1.0, &engine);

  PnpSolverPolicy pnp_solver_policy = {};

  pnp_solver_policy.flags =
      PNP_SOLVER_POLICY_FLAGS_NONE | PNP_SOLVER_POLICY_FLAGS_REFINE;
  pnp_solver_policy.ransac_policy.flags =
      RANSAC_POLICY_FLAGS_EARLY_STOP |
      RANSAC_POLICY_FLAGS_USE_PROBABILITY_WITHOUT_DUPLICATION_SAMPLE;
  pnp_solver_policy.ransac_policy.probability = 0.999999;
  pnp_solver_policy.ransac_policy.num_threads =
      std::thread::hardware_concurrency();
  pnp_solver_policy.ransac_policy.num_ransac_itr_lower_limit = 10;
  pnp_solver_policy.ransac_policy.num_ransac_itr_upper_limit = 1000000;

  PnpInlierCheckParamsUsingBearingVector pnp_inlier_check_params = {};
  pnp_inlier_check_params.inlier_thr                             = 1e-6;
  pnp_solver_policy.pnp_inlier_check_params = pnp_inlier_check_params;

  const std::pair<std::optional<PnpSolverResult<Iterator>>, RansacReport>
      result = SolvePnpProblem(pnp_solver_policy, bearings, points,
                               matches.begin(), matches.end());

  if (result.first) {
    const Mat3d gt_r = camera.rotation_cw;
    const Vec3d gt_t = camera.transltaion_cw;

    std::cout << "R :\n" << result.first->rotation << "\n" << std::endl;
    std::cout << "gt R :\n" << gt_r << "\n" << std::endl;

    std::cout << "t :\n" << result.first->translation << "\n" << std::endl;
    std::cout << "gt t :\n" << gt_t << "\n" << std::endl;
  }
}

}  // namespace gsrap::example

int main(void) {
  RunPnpSolver();
  return 0;
}
