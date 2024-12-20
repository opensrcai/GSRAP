#pragma once
#include <optional>
#include <vector>

#include "common.h"
#include "gsrap/macros.h"
#include "gsrap/pnp_solver.h"

GSRAP_IGNORE_STRICT_WARNING_PUSH
#include "Eigen/Core"
#include "nanobind/nanobind.h"
GSRAP_IGNORE_STRICT_WARNING_POP

namespace gsrap {

struct PyPnpSolverResult {
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  template <typename Iterator, typename Container>
  void From_(const PnpSolverResult<Iterator>&, const Container&);

  /// Rotation matrix
  Eigen::Matrix3d rotation;
  /// Translation matrix
  Eigen::Vector3d translation;

  /// Iterator to inliers
  std::vector<uint32_t> inliers;
  /// Inlier rate
  double inlier_ratio = 0.0;
};

std::pair<std::optional<PyPnpSolverResult>, RansacReport> PySolvePnpProblem(
    PnpSolverPolicy pnp_solver_policy,
    const Eigen::Matrix<double, -1, 3, Eigen::RowMajorBit>& bearings,
    const Eigen::Matrix<double, -1, 3, Eigen::RowMajorBit>& points,
    const RMatrixXu32& pairs);

void WrapPnpSolver(nanobind::module_& module);

}  // namespace gsrap
