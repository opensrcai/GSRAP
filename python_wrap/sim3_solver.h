#pragma once
#include <optional>
#include <vector>

#include "common.h"
#include "gsrap/macros.h"
#include "gsrap/sim3_solver.h"

GSRAP_IGNORE_STRICT_WARNING_PUSH
#include "Eigen/Core"
GSRAP_IGNORE_STRICT_WARNING_POP

namespace gsrap {

struct PySim3SolverResult {
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  /// Rotation matrix
  RMatrix3d rotation;
  /// Translation vector
  Eigen::Vector3d translation;

  /// scale (scalar)
  double scale = 1.0;

  /// Iterator to inliers
  std::vector<uint32_t> inliers;
  /// Inlier rate
  double inlier_ratio = 0.0;
};

std::pair<std::optional<PySim3SolverResult>, RansacReport>
PyComputeSim3Transformation(
    Sim3SolverPolicy sim3_solver_policy,
    const Eigen::Matrix<double, -1, 3, Eigen::RowMajorBit>& points1,
    const Eigen::Matrix<double, -1, 3, Eigen::RowMajorBit>& points2,
    const RMatrixXu32& pairs);

}  // namespace gsrap
