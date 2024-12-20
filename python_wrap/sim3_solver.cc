#include "sim3_solver.h"

#include <vector>

#include "common.h"
#include "gsrap/sim3_solver.h"

namespace gsrap {

using Pair     = std::pair<uint32_t, uint32_t>;
using Iterator = std::vector<Pair>::const_iterator;

std::pair<std::optional<PySim3SolverResult>, RansacReport>
PyComputeSim3Transformation(
    Sim3SolverPolicy sim3_solver_policy,
    const Eigen::Matrix<double, -1, 3, Eigen::RowMajorBit>& points1,
    const Eigen::Matrix<double, -1, 3, Eigen::RowMajorBit>& points2,
    const RMatrixXu32& pairs) {
  eigen_vector<Vec3d> points1_vec(static_cast<size_t>(points1.rows()));
  eigen_vector<Vec3d> points2_vec(static_cast<size_t>(points2.rows()));

  for (int64_t i = 0; i < points1.rows(); ++i) {
    points1_vec[static_cast<size_t>(i)] = points1.row(i);
  }
  for (int64_t i = 0; i < points2.rows(); ++i) {
    points2_vec[static_cast<size_t>(i)] = points2.row(i);
  }

  std::vector<Pair> pairs_vec(static_cast<size_t>(pairs.rows()));

  for (int64_t i = 0; i < pairs.rows(); ++i) {
    pairs_vec[static_cast<size_t>(i)].first  = pairs.row(i).x();
    pairs_vec[static_cast<size_t>(i)].second = pairs.row(i).y();
  }

  const std::pair<std::optional<Sim3SolverResult<Iterator>>, RansacReport> ret =
      ComputeSim3Transformation(sim3_solver_policy, points1_vec, points2_vec,
                                pairs_vec.cbegin(), pairs_vec.cend());

  if (!ret.first.has_value()) {
    return {std::nullopt, ret.second};
  }

  Sim3SolverResult<Iterator> sim3_solver_result =
      ret.first.value_or(Sim3SolverResult<Iterator>());

  PySim3SolverResult py_ret;
  py_ret.rotation    = sim3_solver_result.rotation;
  py_ret.translation = sim3_solver_result.translation;
  py_ret.scale       = sim3_solver_result.scale;

  py_ret.inliers = std::vector<uint32_t>(sim3_solver_result.inliers.size());
  std::transform(sim3_solver_result.inliers.begin(),
                 sim3_solver_result.inliers.end(), py_ret.inliers.begin(),
                 [&](Iterator itr) {
                   return static_cast<uint32_t>(itr - pairs_vec.begin());
                 });

  py_ret.inlier_ratio = sim3_solver_result.inlier_ratio;

  return {py_ret, ret.second};
}

}  // namespace gsrap
