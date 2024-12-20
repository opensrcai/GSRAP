#include <optional>
#include <thread>

// GSRAP
#include "gsrap/macros.h"
#include "gsrap/pnp_solver.h"
// common
#include "common/camera.h"
#include "common/constant.h"
#include "common/type.h"
#include "common/utils.h"

GSRAP_IGNORE_STRICT_WARNING_PUSH

#include "Eigen/Core"
#include "Eigen/Dense"
#include "gtest/gtest.h"

GSRAP_IGNORE_STRICT_WARNING_POP

// NOLINTNEXTLINE
using namespace gsrap::example;

namespace gsrap::test {

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

// This function checks result strictly.
// Check if the error of each element is less than the threshold.
static void CheckResultStrictly(
    const std::pair<std::optional<PnpSolverResult<Iterator>>, RansacReport>&
        result,
    const Camera& camera, const double thr = 1e-13) {
  const Mat3d gt_r = camera.rotation_cw;
  const Vec3d gt_t = camera.transltaion_cw;

  EXPECT_TRUE(result.first.has_value());

  if (!result.first) {
    return;
  }
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      EXPECT_LT(abs(gt_r(i, j) - result.first->rotation(i, j)), thr);
    }
  }
  for (int i = 0; i < 3; ++i) {
    EXPECT_LT(abs(gt_t(i) - result.first->translation(i)), thr);
  }
}

// This function checks result.
static void CheckResult(
    const std::pair<std::optional<PnpSolverResult<Iterator>>, RansacReport>&
        result,
    const Camera& camera,
    const double theta_thr = kPi * 1.0 / 180.0 /* 1.0deg*/,
    const double dist      = 1e-2) {
  const Mat3d gt_r = camera.rotation_cw;
  const Vec3d gt_t = camera.transltaion_cw;
  EXPECT_TRUE(result.first.has_value());

  if (result.first) {
    const Mat3d relative_r = gt_r * result.first->rotation.transpose();
    const double theta     = AxisAngleTheta(relative_r);
    EXPECT_LT(theta, theta_thr);

    EXPECT_LT((gt_t - result.first->translation).norm(), dist);
  }
}

// NOLINTNEXTLINE
TEST(ransac, pnp_solver) {
  const uint32_t num_points = 500;  // Number of points

  const Vec3d center(0.5 / sqrt(2.0), -0.5 / sqrt(2.0),
                     0.5 / sqrt(2.0));    // Center of spherical point cloud.
  const double radius = 0.5 / sqrt(2.0);  // radius of spherical point cloud

  std::mt19937 engine(123456789);
  const uint32_t num_test_iteration = 1000;

  for (uint32_t itr_cnt = 0; itr_cnt < num_test_iteration; ++itr_cnt) {
    eigen_vector<Vec3d> points =
        CreatesphericalPointCloud(center, radius, num_points, &engine);

    const auto camera = CreateCameras();

    const auto bearings = ToBearingVectors(points, camera);

    const std::vector<Pair> matches = CreateMatches(num_points, 1.0, &engine);

    PnpSolverPolicy pnp_solver_policy = {};

    pnp_solver_policy.flags =
        PNP_SOLVER_POLICY_FLAGS_NONE | PNP_SOLVER_POLICY_FLAGS_REFINE;
    pnp_solver_policy.ransac_policy = {};
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

    CheckResultStrictly(result, camera);
  }
}

// NOLINTNEXTLINE
TEST(ransac, pnp_solver_with_miss_matches) {
  const uint32_t num_points = 500;  // Number of points
  const double inlier_ratio = 0.75;

  const Vec3d center(0.5 / sqrt(2.0), -0.5 / sqrt(2.0),
                     0.5 / sqrt(2.0));    // Center of spherical point cloud.
  const double radius = 0.5 / sqrt(2.0);  // radius of spherical point cloud

  std::mt19937 engine(123456789);
  const uint32_t num_test_iteration = 500;

  for (uint32_t itr_cnt = 0; itr_cnt < num_test_iteration; ++itr_cnt) {
    eigen_vector<Vec3d> points =
        CreatesphericalPointCloud(center, radius, num_points, &engine);

    const auto camera = CreateCameras();

    const auto bearings = ToBearingVectors(points, camera);

    const std::vector<Pair> matches =
        CreateMatches(num_points, inlier_ratio, &engine);

    PnpSolverPolicy pnp_solver_policy = {};

    pnp_solver_policy.flags =
        PNP_SOLVER_POLICY_FLAGS_NONE | PNP_SOLVER_POLICY_FLAGS_REFINE;
    pnp_solver_policy.ransac_policy = {};
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

    CheckResultStrictly(result, camera);
  }
}

// NOLINTNEXTLINE
TEST(ransac, pnp_solver_with_noise) {
  const uint32_t num_points = 500;  // Number of points
  const double inlier_ratio = 0.80;

  const Vec3d center(0.5 / sqrt(2.0), -0.5 / sqrt(2.0),
                     0.5 / sqrt(2.0));    // Center of spherical point cloud.
  const double radius = 0.5 / sqrt(2.0);  // Radius of spherical point cloud

  std::mt19937 engine(123456789);
  const uint32_t num_test_iteration = 500;

  for (uint32_t itr_cnt = 0; itr_cnt < num_test_iteration; ++itr_cnt) {
    eigen_vector<Vec3d> points =
        CreatesphericalPointCloud(center, radius, num_points, &engine);

    const auto camera = CreateCameras();

    auto bearings = ToBearingVectors(points, camera);

    RotateBearingVectors(kPi * 0.1 / 180.0, &bearings, &engine);  // 0.1 degree
    RotateSampledBearingVectors(kPi * 10.0 / 180.0, 1.0 - inlier_ratio,
                                &bearings, &engine);  // more than 10 degree

    const std::vector<Pair> matches = CreateMatches(num_points, 1.0, &engine);

    PnpSolverPolicy pnp_solver_policy = {};

    pnp_solver_policy.flags =
        PNP_SOLVER_POLICY_FLAGS_NONE | PNP_SOLVER_POLICY_FLAGS_REFINE;
    pnp_solver_policy.ransac_policy = {};
    pnp_solver_policy.ransac_policy.flags =
        RANSAC_POLICY_FLAGS_EARLY_STOP |
        RANSAC_POLICY_FLAGS_USE_PROBABILITY_WITHOUT_DUPLICATION_SAMPLE;
    pnp_solver_policy.ransac_policy.probability = 0.999;
    pnp_solver_policy.ransac_policy.num_threads =
        std::thread::hardware_concurrency();
    pnp_solver_policy.ransac_policy.num_ransac_itr_lower_limit = 10;
    pnp_solver_policy.ransac_policy.num_ransac_itr_upper_limit = 1000000;

    PnpInlierCheckParamsUsingBearingVector pnp_inlier_check_params = {};
    pnp_inlier_check_params.inlier_thr                             = 5e-4;
    pnp_solver_policy.pnp_inlier_check_params = pnp_inlier_check_params;

    const std::pair<std::optional<PnpSolverResult<Iterator>>, RansacReport>
        result = SolvePnpProblem(pnp_solver_policy, bearings, points,
                                 matches.begin(), matches.end());

    CheckResult(result, camera, kPi * 0.5 / 180.0 /* 0.5 degree */);
  }
}

// NOLINTNEXTLINE
TEST(ransac, pnp_solver_with_miss_matches_and_noise) {
  const uint32_t num_points = 500;  // Number of points

  const Vec3d center(0.5 / sqrt(2.0), -0.5 / sqrt(2.0),
                     0.5 / sqrt(2.0));    // Center of spherical point cloud.
  const double radius = 0.5 / sqrt(2.0);  // Radius of spherical point cloud

  std::mt19937 engine(123456789);
  const uint32_t num_test_iteration = 500;

  for (uint32_t itr_cnt = 0; itr_cnt < num_test_iteration; ++itr_cnt) {
    eigen_vector<Vec3d> points =
        CreatesphericalPointCloud(center, radius, num_points, &engine);

    const auto camera = CreateCameras();

    auto bearings = ToBearingVectors(points, camera);

    RotateBearingVectors(kPi * 0.1 / 180.0, &bearings, &engine);  // 0.1 degree
    RotateSampledBearingVectors(
        kPi * 10.0 / 180.0, 0.20, &bearings,
        &engine);  // 20% of bearing vectors have more than 10 degree errors.

    const std::vector<Pair> matches = CreateMatches(num_points, 0.80, &engine);

    PnpSolverPolicy pnp_solver_policy = {};

    pnp_solver_policy.flags =
        PNP_SOLVER_POLICY_FLAGS_NONE | PNP_SOLVER_POLICY_FLAGS_REFINE;
    pnp_solver_policy.ransac_policy = {};
    pnp_solver_policy.ransac_policy.flags =
        RANSAC_POLICY_FLAGS_EARLY_STOP |
        RANSAC_POLICY_FLAGS_USE_PROBABILITY_WITHOUT_DUPLICATION_SAMPLE;
    pnp_solver_policy.ransac_policy.probability = 0.999;
    pnp_solver_policy.ransac_policy.num_threads =
        std::thread::hardware_concurrency();
    pnp_solver_policy.ransac_policy.num_ransac_itr_lower_limit = 10;
    pnp_solver_policy.ransac_policy.num_ransac_itr_upper_limit = 1000000;

    PnpInlierCheckParamsUsingBearingVector pnp_inlier_check_params = {};
    pnp_inlier_check_params.inlier_thr                             = 2e-3;
    pnp_solver_policy.pnp_inlier_check_params = pnp_inlier_check_params;

    const std::pair<std::optional<PnpSolverResult<Iterator>>, RansacReport>
        result = SolvePnpProblem(pnp_solver_policy, bearings, points,
                                 matches.begin(), matches.end());

    CheckResult(result, camera, kPi * 0.5 / 180.0 /* 0.5 degree */, 1e-3);
  }
}

}  // namespace gsrap::test
