#pragma once

#include "common.h"
#include "gsrap/essential_solver.h"
#include "gsrap/macros.h"

GSRAP_IGNORE_STRICT_WARNING_PUSH
#include "Eigen/Core"
#include "nanobind/nanobind.h"
GSRAP_IGNORE_STRICT_WARNING_POP

namespace gsrap {

struct PyEssentialSolverResult {
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  template <typename Iterator, typename Container>
  void From_(const EssentialSolverResult<Iterator>&, const Container&);

  /// Rotation matrix
  Eigen::Matrix3d rotation;
  /// Translation matrix
  Eigen::Vector3d translation;

  /// Iterator to inliers
  std::vector<uint32_t> inliers;
  /// Inlier rate
  double inlier_ratio = 0.0;
};

void WrapEssentialSolver(nanobind::module_& module);

}  // namespace gsrap
