#include <stdint.h>

#include <thread>

// GSRAP
#include "gsrap/macros.h"
#include "gsrap/sim3_solver.h"
// common
#include "common/utils.h"

GSRAP_IGNORE_STRICT_WARNING_PUSH

#include "gtest/gtest.h"

GSRAP_IGNORE_STRICT_WARNING_POP

// NOLINTNEXTLINE
using namespace gsrap::example;

namespace gsrap::test {

using Pair     = std::pair<uint32_t, uint32_t>;
using Iterator = std::vector<Pair>::const_iterator;

// This function checks result strictly.
// Check if the error of each element is less than the threshold.
static void CheckResultStrictly(
    const std::pair<std::optional<Sim3SolverResult<Iterator>>, RansacReport>&
        result,
    const double gt_scale, const Mat3d& gt_r, const Vec3d& gt_t,
    const double thr = 1e-14) {
  EXPECT_TRUE(result.first.has_value());

  if (!result.first) {
    return;
  }
  EXPECT_LT(abs(gt_scale - result.first->scale), thr);
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
    const std::pair<std::optional<Sim3SolverResult<Iterator>>, RansacReport>&
        result,
    const double gt_scale, const Mat3d& gt_r, const Vec3d& gt_t,
    const double theta_thr = kPi * 1.0 / 180.0 /* 1.0deg*/,
    const double scale_thr = 1e-2, const double dist_thr = 1e-2) {
  EXPECT_TRUE(result.first.has_value());

  if (result.first) {
    EXPECT_LT(abs(gt_scale - result.first->scale) / gt_scale, scale_thr);

    const Mat3d relative_r = gt_r * result.first->rotation.transpose();
    const double theta     = AxisAngleTheta(relative_r);
    EXPECT_LT(theta, theta_thr);

    EXPECT_LT((gt_t - result.first->translation).norm() / gt_t.norm(),
              dist_thr);
  }
}

// NOLINTNEXTLINE
TEST(ransac, sim3_solver) {
  const uint32_t num_points = 500;  // Number of points

  const Vec3d center(0.0, 0.0, 0.0);  // Center of spherical point cloud.
  const double edge_length = 1;       // Edge length of cuboid point cloud

  const double gt_scale = 1.25;
  Mat3d gt_r =
      Eigen::AngleAxis<double>(
          0.25 * kPi,
          (Vec3d::UnitX() + Vec3d::UnitY() + Vec3d::UnitZ()).normalized())
          .toRotationMatrix();
  const Vec3d gt_t(1.0, 0.25, 0.5);

  std::mt19937 engine(123456789);
  const uint32_t num_test_iteration = 10000;

  for (uint32_t itr_cnt = 0; itr_cnt < num_test_iteration; ++itr_cnt) {
    eigen_vector<Vec3d> points1 = CreateCuboidlPointCloud(
        center, edge_length, edge_length, edge_length, num_points, &engine);

    eigen_vector<Vec3d> points2;
    points2.resize(num_points);
    std::transform(points1.begin(), points1.end(), points2.begin(),
                   [&gt_scale, &gt_r, &gt_t](const Vec3d& p1) {
                     return gt_scale * gt_r * p1 + gt_t;
                   });

    const std::vector<Pair> matches = CreateMatches(num_points, 1.0, &engine);

    Sim3SolverPolicy sim3_solver_policy = {};

    sim3_solver_policy.flags =
        SIM3_SOLVER_POLICY_FLAGS_NONE | SIM3_SOLVER_POLICY_FLAGS_REFINE;
    sim3_solver_policy.ransac_policy = {};
    sim3_solver_policy.ransac_policy.flags =
        RANSAC_POLICY_FLAGS_EARLY_STOP |
        RANSAC_POLICY_FLAGS_USE_PROBABILITY_WITHOUT_DUPLICATION_SAMPLE;
    sim3_solver_policy.ransac_policy.probability = 0.999999;
    sim3_solver_policy.ransac_policy.num_threads =
        std::thread::hardware_concurrency();
    sim3_solver_policy.ransac_policy.num_ransac_itr_lower_limit = 10;
    sim3_solver_policy.ransac_policy.num_ransac_itr_upper_limit = 1000000;

    sim3_solver_policy.inlier_thr = 1e-12;

    const std::pair<std::optional<Sim3SolverResult<Iterator>>, RansacReport>
        result = ComputeSim3Transformation(sim3_solver_policy, points1, points2,
                                           matches.begin(), matches.end());

    CheckResultStrictly(result, gt_scale, gt_r, gt_t);
    if (result.first) {
      EXPECT_DOUBLE_EQ(result.first->inlier_ratio, 1.0);
      EXPECT_EQ(result.first->inliers.size(), matches.size());
    }
  }
}

// NOLINTNEXTLINE
TEST(ransac, sim3_solver_with_miss_matches) {
  const uint32_t num_points = 500;  // Number of points
  const double inlier_ratio  = 0.75;

  const Vec3d center(0.0, 0.0, 0.0);  // Center of spherical point cloud.
  const double edge_length = 1;       // Edge length of cuboid point cloud

  const double gt_scale = 1.25;
  Mat3d gt_r =
      Eigen::AngleAxis<double>(
          0.25 * kPi,
          (Vec3d::UnitX() + Vec3d::UnitY() + Vec3d::UnitZ()).normalized())
          .toRotationMatrix();
  const Vec3d gt_t(1.0, 0.25, 0.5);

  std::mt19937 engine(123456789);
  const uint32_t num_test_iteration = 10000;

  for (uint32_t itr_cnt = 0; itr_cnt < num_test_iteration; ++itr_cnt) {
    eigen_vector<Vec3d> points1 = CreateCuboidlPointCloud(
        center, edge_length, edge_length, edge_length, num_points, &engine);

    eigen_vector<Vec3d> points2;
    points2.resize(num_points);
    std::transform(points1.begin(), points1.end(), points2.begin(),
                   [&gt_scale, &gt_r, &gt_t](const Vec3d& p1) {
                     return gt_scale * gt_r * p1 + gt_t;
                   });

    const std::vector<Pair> matches =
        CreateMatches(num_points, inlier_ratio, &engine);

    Sim3SolverPolicy sim3_solver_policy = {};

    sim3_solver_policy.flags =
        SIM3_SOLVER_POLICY_FLAGS_NONE | SIM3_SOLVER_POLICY_FLAGS_REFINE;
    sim3_solver_policy.ransac_policy = {};
    sim3_solver_policy.ransac_policy.flags =
        RANSAC_POLICY_FLAGS_EARLY_STOP |
        RANSAC_POLICY_FLAGS_USE_PROBABILITY_WITHOUT_DUPLICATION_SAMPLE;
    sim3_solver_policy.ransac_policy.probability = 0.999999;
    sim3_solver_policy.ransac_policy.num_threads =
        std::thread::hardware_concurrency();
    sim3_solver_policy.ransac_policy.num_ransac_itr_lower_limit = 10;
    sim3_solver_policy.ransac_policy.num_ransac_itr_upper_limit = 1000000;

    sim3_solver_policy.inlier_thr = 1e-12;

    const std::pair<std::optional<Sim3SolverResult<Iterator>>, RansacReport>
        result = ComputeSim3Transformation(sim3_solver_policy, points1, points2,
                                           matches.begin(), matches.end());

    CheckResultStrictly(result, gt_scale, gt_r, gt_t);
    if (result.first) {
      EXPECT_DOUBLE_EQ(result.first->inlier_ratio, inlier_ratio);
    }
  }
}

// NOLINTNEXTLINE
TEST(ransac, sim3_solver_with_noise) {
  const uint32_t num_points  = 500;  // Number of points
  const double noise_scale   = 0.025;
  const double outlier_scale = 10.0;
  const double inlier_ratio   = 0.75;

  const Vec3d center(0.0, 0.0, 0.0);  // Center of spherical point cloud.
  const double edge_length = 1;       // Edge length of cuboid point cloud

  const double gt_scale = 1.25;
  Mat3d gt_r =
      Eigen::AngleAxis<double>(
          0.25 * kPi,
          (Vec3d::UnitX() + Vec3d::UnitY() + Vec3d::UnitZ()).normalized())
          .toRotationMatrix();
  const Vec3d gt_t(1.0, 0.25, 0.5);

  std::mt19937 engine(123456789);
  const uint32_t num_test_iteration = 5000;

  for (uint32_t itr_cnt = 0; itr_cnt < num_test_iteration; ++itr_cnt) {
    eigen_vector<Vec3d> points1 = CreateCuboidlPointCloud(
        center, edge_length, edge_length, edge_length, num_points, &engine);

    eigen_vector<Vec3d> points2;
    points2.resize(num_points);
    std::transform(points1.begin(), points1.end(), points2.begin(),
                   [&gt_scale, &gt_r, &gt_t](const Vec3d& p1) {
                     return gt_scale * gt_r * p1 + gt_t;
                   });

    AddNoiseToPoints(0.0, noise_scale * edge_length, 1.0, &points2, &engine);

    AddNoiseToPoints(outlier_scale * edge_length, edge_length,
                     1.0 - inlier_ratio, &points2, &engine);

    const std::vector<Pair> matches = CreateMatches(num_points, 1.0, &engine);

    Sim3SolverPolicy sim3_solver_policy = {};

    sim3_solver_policy.flags =
        SIM3_SOLVER_POLICY_FLAGS_NONE | SIM3_SOLVER_POLICY_FLAGS_REFINE;
    sim3_solver_policy.ransac_policy = {};
    sim3_solver_policy.ransac_policy.flags =
        RANSAC_POLICY_FLAGS_EARLY_STOP |
        RANSAC_POLICY_FLAGS_USE_PROBABILITY_WITHOUT_DUPLICATION_SAMPLE;
    sim3_solver_policy.ransac_policy.probability                = 0.999999;
    sim3_solver_policy.ransac_policy.num_threads =
        std::thread::hardware_concurrency();
    sim3_solver_policy.ransac_policy.num_ransac_itr_lower_limit = 10;
    sim3_solver_policy.ransac_policy.num_ransac_itr_upper_limit = 1000000;

    sim3_solver_policy.inlier_thr = noise_scale * 2;
    sim3_solver_policy.num_sample = 3;

    const std::pair<std::optional<Sim3SolverResult<Iterator>>, RansacReport>
        result = ComputeSim3Transformation(sim3_solver_policy, points1, points2,
                                           matches.begin(), matches.end());

    CheckResult(result, gt_scale, gt_r, gt_t, kPi * 1.0 / 180.0);
  }
}

// NOLINTNEXTLINE
TEST(ransac, sim3_solver_with_miss_matches_and_noise) {
  const uint32_t num_points  = 500;  // Number of points
  const double noise_scale   = 0.025;
  const double outlier_scale = 10.0;

  const Vec3d center(0.0, 0.0, 0.0);  // Center of spherical point cloud.
  const double edge_length = 1;       // Edge length of cuboid point cloud

  const double gt_scale = 1.25;
  Mat3d gt_r =
      Eigen::AngleAxis<double>(
          0.25 * kPi,
          (Vec3d::UnitX() + Vec3d::UnitY() + Vec3d::UnitZ()).normalized())
          .toRotationMatrix();
  const Vec3d gt_t(1.0, 0.25, 0.5);

  std::mt19937 engine(123456789);
  const uint32_t num_test_iteration = 5000;

  for (uint32_t itr_cnt = 0; itr_cnt < num_test_iteration; ++itr_cnt) {
    eigen_vector<Vec3d> points1 = CreateCuboidlPointCloud(
        center, edge_length, edge_length, edge_length, num_points, &engine);

    eigen_vector<Vec3d> points2;
    points2.resize(num_points);
    std::transform(points1.begin(), points1.end(), points2.begin(),
                   [&gt_scale, &gt_r, &gt_t](const Vec3d& p1) {
                     return gt_scale * gt_r * p1 + gt_t;
                   });

    AddNoiseToPoints(0.0, noise_scale * edge_length, 1.0, &points2, &engine);

    AddNoiseToPoints(outlier_scale * edge_length, edge_length, 0.1, &points2,
                     &engine);

    const std::vector<Pair> matches = CreateMatches(num_points, 0.75, &engine);

    Sim3SolverPolicy sim3_solver_policy = {};

    sim3_solver_policy.flags =
        SIM3_SOLVER_POLICY_FLAGS_NONE | SIM3_SOLVER_POLICY_FLAGS_REFINE;
    sim3_solver_policy.ransac_policy = {};
    sim3_solver_policy.ransac_policy.flags =
        RANSAC_POLICY_FLAGS_EARLY_STOP |
        RANSAC_POLICY_FLAGS_USE_PROBABILITY_WITHOUT_DUPLICATION_SAMPLE;
    sim3_solver_policy.ransac_policy.probability                = 0.999999;
    sim3_solver_policy.ransac_policy.num_threads =
        std::thread::hardware_concurrency();
    sim3_solver_policy.ransac_policy.num_ransac_itr_lower_limit = 10;
    sim3_solver_policy.ransac_policy.num_ransac_itr_upper_limit = 1000000;

    sim3_solver_policy.inlier_thr = noise_scale * 2;
    sim3_solver_policy.num_sample = 3;

    const std::pair<std::optional<Sim3SolverResult<Iterator>>, RansacReport>
        result = ComputeSim3Transformation(sim3_solver_policy, points1, points2,
                                           matches.begin(), matches.end());

    CheckResult(result, gt_scale, gt_r, gt_t, kPi * 1.0 / 180.0);
  }
}

}  // namespace gsrap::test
