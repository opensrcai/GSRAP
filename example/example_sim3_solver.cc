#include <iostream>
#include <random>
#include <thread>

// GSRAP
#include "gsrap/macros.h"
#include "gsrap/sim3_solver.h"

// common
#include "common/utils.h"

// NOLINTNEXTLINE
using namespace gsrap;

// NOLINTNEXTLINE
using namespace gsrap::example;

namespace gsrap::example {

using Pair     = std::pair<uint32_t, uint32_t>;
using Iterator = std::vector<Pair>::const_iterator;

static void RunSim3Solver(void) {
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

  if (result.first) {
    std::cout << "scale:\n" << result.first->scale << "\n" << std::endl;
    std::cout << "gt scale:\n" << gt_scale << "\n" << std::endl;

    std::cout << "R :\n" << result.first->rotation << "\n" << std::endl;
    std::cout << "gt R :\n" << gt_r << "\n" << std::endl;

    std::cout << "t :\n" << result.first->translation << "\n" << std::endl;
    std::cout << "gt t :\n" << gt_t << "\n" << std::endl;

    std::cout << "inlier rate: " << result.first->inlier_ratio << "\n"
              << std::endl;
  }

  std::cout << "RANSAC itrations: " << result.second.num_iteration << "\n"
            << std::endl;
}

}  // namespace gsrap::example

int main(void) {
  RunSim3Solver();
  return 0;
}
