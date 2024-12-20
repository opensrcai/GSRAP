#include <stdint.h>

#include <algorithm>
#include <random>
#include <thread>
#include <unordered_set>
#include <utility>
#include <vector>

// GSRAP
#include "gsrap/essential_solver.h"
#include "gsrap/macros.h"
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
    const Vec3d camera2_xw =
        camera2_up_w.cross(camera2_optical_axis_w).normalized();
    const Vec3d camera2_yw =
        camera2_optical_axis_w.cross(camera2_xw).normalized();

    // Create world to camera rotation matrix.
    camera2.rotation_cw.row(0) = camera2_xw;
    camera2.rotation_cw.row(1) = camera2_yw;
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

// This function checks result strictly.
// Check if the error of each element is less than the threshold.
static void CheckResultStrictly(
    const std::pair<std::optional<EssentialSolverResult<Iterator>>,
                    RansacReport>& result,
    const Camera& camera1, const Camera& camera2, const double thr = 1e-14) {
  const Mat3d gt_r = camera2.rotation_cw * camera1.rotation_cw.transpose();
  const Vec3d gt_t = camera2.transltaion_cw - camera1.transltaion_cw;
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
    const std::pair<std::optional<EssentialSolverResult<Iterator>>,
                    RansacReport>& result,
    const Camera& camera1, const Camera& camera2,
    const double theta_thr = kPi * 1.0 / 180.0 /* 1.0deg*/,
    const double dist      = 1e-2) {
  const Mat3d gt_r = camera2.rotation_cw * camera1.rotation_cw.transpose();
  const Vec3d gt_t = camera2.transltaion_cw - camera1.transltaion_cw;
  EXPECT_TRUE(result.first.has_value());

  if (result.first) {
    const Mat3d relative_r = gt_r * result.first->rotation.transpose();
    const double theta     = AxisAngleTheta(relative_r);
    EXPECT_LT(theta, theta_thr);

    EXPECT_LT((gt_t - result.first->translation).norm(), dist);
  }
}

// NOLINTNEXTLINE
TEST(ransac, essential_solver) {
  const uint32_t num_points = 500;  // Number of points

  const Vec3d center(0.5 / sqrt(2.0), -0.5 / sqrt(2.0),
                     0.5 / sqrt(2.0));    // Center of spherical point cloud.
  const double radius = 0.5 / sqrt(2.0);  // radius of spherical point cloud

  std::mt19937 engine(123456789);
  const uint32_t num_test_iteration = 1000;

  for (uint32_t itr_cnt = 0; itr_cnt < num_test_iteration; ++itr_cnt) {
    eigen_vector<Vec3d> points =
        CreatesphericalPointCloud(center, radius, num_points, &engine);

    const auto [camera1, camera2] = CreateTwoCameras();

    const auto [b1s, b2s] = ToBearingVectors(points, camera1, camera2);

    const std::vector<Pair> matches =
        CreateMatches(num_points, /* inlier_ratio */ 1.0, &engine);

    EssentialSolverPolicy essential_solver_policy = {};

    essential_solver_policy.flags =
        ESSENTIAL_SOLVER_POLICY_FLAGS_NONE |
        ESSENTIAL_SOLVER_POLICY_FLAGS_CHECK_SINGULAR_VALUE;
    essential_solver_policy.ransac_policy.flags =
        RANSAC_POLICY_FLAGS_EARLY_STOP |
        RANSAC_POLICY_FLAGS_USE_PROBABILITY_WITHOUT_DUPLICATION_SAMPLE;
    essential_solver_policy.ransac_policy.probability = 0.999999;
    essential_solver_policy.ransac_policy.num_threads =
        std::thread::hardware_concurrency();
    essential_solver_policy.ransac_policy.num_ransac_itr_lower_limit = 30;
    essential_solver_policy.ransac_policy.num_ransac_itr_upper_limit = 100;

    const std::pair<std::optional<EssentialSolverResult<Iterator>>,
                    RansacReport>
        result = ComputeEssentialMatrix(essential_solver_policy, b1s, b2s,
                                        matches.cbegin(), matches.cend());

    CheckResultStrictly(result, camera1, camera2);
    if (result.first) {
      EXPECT_DOUBLE_EQ(result.first->inlier_ratio, 1.0);
      EXPECT_EQ(result.first->inliers.size(), matches.size());
    }
  }
}

// NOLINTNEXTLINE
TEST(ransac, essential_solver_with_miss_matches) {
  const uint32_t num_points = 500;  // Number of points
  const double inlier_ratio = 0.2;

  const Vec3d center(0.5 / sqrt(2.0), -0.5 / sqrt(2.0),
                     0.5 / sqrt(2.0));    // Center of spherical point cloud.
  const double radius = 0.5 / sqrt(2.0);  // radius of spherical point cloud

  std::mt19937 engine(123456789);
  const uint32_t num_test_iteration = 50;

  for (uint32_t itr_cnt = 0; itr_cnt < num_test_iteration; ++itr_cnt) {
    eigen_vector<Vec3d> points =
        CreatesphericalPointCloud(center, radius, num_points, &engine);
    const auto num_miss_matches =
        uint32_t(std::round(double(num_points) * (1.0 - inlier_ratio)));
    const auto [camera1, camera2] = CreateTwoCameras();

    const auto [b1s, b2s] = ToBearingVectors(points, camera1, camera2);

    const std::vector<Pair> matches =
        CreateMatches(num_points, inlier_ratio, &engine);

    EssentialSolverPolicy essential_solver_policy = {};

    essential_solver_policy.flags =
        ESSENTIAL_SOLVER_POLICY_FLAGS_NONE |
        ESSENTIAL_SOLVER_POLICY_FLAGS_CHECK_SINGULAR_VALUE;
    essential_solver_policy.ransac_policy.flags =
        RANSAC_POLICY_FLAGS_EARLY_STOP |
        RANSAC_POLICY_FLAGS_USE_PROBABILITY_WITHOUT_DUPLICATION_SAMPLE;
    essential_solver_policy.ransac_policy.probability = 0.9999;
    essential_solver_policy.ransac_policy.num_threads =
        std::thread::hardware_concurrency();
    essential_solver_policy.ransac_policy.num_ransac_itr_lower_limit = 10;
    essential_solver_policy.ransac_policy.num_ransac_itr_upper_limit = 1000000;

    const std::pair<std::optional<EssentialSolverResult<Iterator>>,
                    RansacReport>
        result = ComputeEssentialMatrix(essential_solver_policy, b1s, b2s,
                                        matches.cbegin(), matches.cend());

    // printf("result  %s\n", (result.second.ransac_termination_info ==
    //                         gsrap::RansacTerminationInfo::EARLY_STOPED)
    //                            ? "early stoped"
    //                            : "reached limit");

    CheckResultStrictly(result, camera1, camera2);
    EXPECT_DOUBLE_EQ(result.first->inlier_ratio, inlier_ratio);
    EXPECT_EQ(result.first->inliers.size(), matches.size() - num_miss_matches);
  }
}

// NOLINTNEXTLINE
TEST(ransac, essential_solver_with_noise) {
  const uint32_t num_points = 500;  // Number of points

  const Vec3d center(0.5 / sqrt(2.0), -0.5 / sqrt(2.0),
                     0.5 / sqrt(2.0));    // Center of spherical point cloud.
  const double radius = 0.5 / sqrt(2.0);  // radius of spherical point cloud

  std::mt19937 engine(123456789);
  const uint32_t num_test_iteration = 50;

  for (uint32_t itr_cnt = 0; itr_cnt < num_test_iteration; ++itr_cnt) {
    eigen_vector<Vec3d> points =
        CreatesphericalPointCloud(center, radius, num_points, &engine);
    const auto [camera1, camera2] = CreateTwoCameras();

    auto [b1s, b2s] = ToBearingVectors(points, camera1, camera2);

    RotateBearingVectors(kPi * 0.5 / 180.0, &b1s, &engine);  // 0.5 degree
    RotateBearingVectors(kPi * 0.5 / 180.0, &b2s, &engine);  // 0.5 degree

    RotateSampledBearingVectors(kPi * 10.0 / 180.0, 0.25, &b1s, &b2s, &engine);

    const std::vector<Pair> matches =
        CreateMatches(num_points, /* inlier_ratio */ 1.0, &engine);

    EssentialSolverPolicy essential_solver_policy = {};

    essential_solver_policy.flags =
        ESSENTIAL_SOLVER_POLICY_FLAGS_NONE |
        ESSENTIAL_SOLVER_POLICY_FLAGS_CHECK_SINGULAR_VALUE;
    essential_solver_policy.ransac_policy.flags =
        RANSAC_POLICY_FLAGS_EARLY_STOP |
        RANSAC_POLICY_FLAGS_USE_PROBABILITY_WITHOUT_DUPLICATION_SAMPLE;
    essential_solver_policy.ransac_policy.probability = 0.999;
    essential_solver_policy.ransac_policy.num_threads =
        std::thread::hardware_concurrency();
    essential_solver_policy.ransac_policy.num_ransac_itr_lower_limit = 10;
    essential_solver_policy.ransac_policy.num_ransac_itr_upper_limit = 1000000;
    essential_solver_policy.inlier_thr                               = 1e-3;

    const std::pair<std::optional<EssentialSolverResult<Iterator>>,
                    RansacReport>
        result = ComputeEssentialMatrix(essential_solver_policy, b1s, b2s,
                                        matches.cbegin(), matches.cend());

    CheckResult(result, camera1, camera2, kPi * 3.0 / 180.0 /* 3 degree */,
                2.5e-2);
  }
}

// NOLINTNEXTLINE
TEST(ransac, essential_solver_with_miss_matches_and_noise) {
  const uint32_t num_points = 500;  // Number of points

  const Vec3d center(0.5 / sqrt(2.0), -0.5 / sqrt(2.0),
                     0.5 / sqrt(2.0));    // Center of spherical point cloud.
  const double radius = 0.5 / sqrt(2.0);  // radius of spherical point cloud

  std::mt19937 engine(123456789);
  const uint32_t num_test_iteration = 50;

  for (uint32_t itr_cnt = 0; itr_cnt < num_test_iteration; ++itr_cnt) {
    eigen_vector<Vec3d> points =
        CreatesphericalPointCloud(center, radius, num_points, &engine);
    const auto [camera1, camera2] = CreateTwoCameras();

    auto [b1s, b2s] = ToBearingVectors(points, camera1, camera2);

    RotateBearingVectors(kPi * 0.5 / 180.0, &b1s, &engine);  // 0.5 degree
    RotateBearingVectors(kPi * 0.5 / 180.0, &b2s, &engine);  // 0.5 degree

    RotateSampledBearingVectors(kPi * 10.0 / 180.0, 0.25, &b1s, &b2s, &engine);

    const std::vector<Pair> matches =
        CreateMatches(num_points, /* inlier_ratio */ 0.75, &engine);

    EssentialSolverPolicy essential_solver_policy = {};

    essential_solver_policy.flags =
        ESSENTIAL_SOLVER_POLICY_FLAGS_NONE |
        ESSENTIAL_SOLVER_POLICY_FLAGS_CHECK_SINGULAR_VALUE;
    essential_solver_policy.ransac_policy.flags =
        RANSAC_POLICY_FLAGS_EARLY_STOP |
        RANSAC_POLICY_FLAGS_USE_PROBABILITY_WITHOUT_DUPLICATION_SAMPLE;
    essential_solver_policy.ransac_policy.probability = 0.99;
    essential_solver_policy.ransac_policy.num_threads =
        std::thread::hardware_concurrency();
    essential_solver_policy.ransac_policy.num_ransac_itr_lower_limit = 10;
    essential_solver_policy.ransac_policy.num_ransac_itr_upper_limit = 1000000;
    essential_solver_policy.inlier_thr                               = 1e-3;

    const std::pair<std::optional<EssentialSolverResult<Iterator>>,
                    RansacReport>
        result = ComputeEssentialMatrix(essential_solver_policy, b1s, b2s,
                                        matches.cbegin(), matches.cend());

    CheckResult(result, camera1, camera2, kPi * 3.0 / 180.0 /* 3 degree */,
                2.5e-2);
  }
}

}  // namespace gsrap::test
