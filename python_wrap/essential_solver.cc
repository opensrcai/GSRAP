#include "essential_solver.h"

#include "common.h"
#include "gsrap/essential_solver.h"
#include "gsrap/macros.h"

GSRAP_IGNORE_STRICT_WARNING_PUSH
#include "nanobind/eigen/dense.h"
#include "nanobind/stl/optional.h"
#include "nanobind/stl/pair.h"
#include "nanobind/stl/vector.h"
GSRAP_IGNORE_STRICT_WARNING_POP

namespace gsrap {

using Pair     = std::pair<uint32_t, uint32_t>;
using Iterator = std::vector<Pair>::const_iterator;

template <typename Iterator, typename Container>
void PyEssentialSolverResult::From_(
    const EssentialSolverResult<Iterator>& result, const Container& pairs) {
  rotation    = result.rotation;
  translation = result.translation;

  const auto& _inliers = result.inliers;
  inliers.resize(_inliers.size());
  std::transform(
      _inliers.cbegin(), _inliers.cend(), inliers.begin(),
      [&](Iterator itr) { return static_cast<uint32_t>(itr - pairs.begin()); });

  inlier_ratio = result.inlier_ratio;
}

std::pair<std::optional<PyEssentialSolverResult>, RansacReport>
PyComputeEssentialMatrix(
    EssentialSolverPolicy essential_solver_policy,
    const Eigen::Matrix<double, -1, 3, Eigen::RowMajorBit>& bearings1,
    const Eigen::Matrix<double, -1, 3, Eigen::RowMajorBit>& bearings2,
    const RMatrixXu32& pairs) {
  eigen_vector<Vec3d> bearings1_vec(static_cast<size_t>(bearings1.rows()));
  eigen_vector<Vec3d> bearings2_vec(static_cast<size_t>(bearings2.rows()));

  for (int64_t i = 0; i < bearings1.rows(); ++i) {
    bearings1_vec[static_cast<size_t>(i)] = bearings1.row(i);
  }
  for (int64_t i = 0; i < bearings2.rows(); ++i) {
    bearings2_vec[static_cast<size_t>(i)] = bearings2.row(i);
  }

  std::vector<Pair> pairs_vec(static_cast<size_t>(pairs.rows()));

  for (int64_t i = 0; i < pairs.rows(); ++i) {
    pairs_vec[static_cast<size_t>(i)].first  = pairs.row(i).x();
    pairs_vec[static_cast<size_t>(i)].second = pairs.row(i).y();
  }

  const auto [opt_essential_solver_result, ransac_result] =
      ComputeEssentialMatrix(essential_solver_policy, bearings1_vec,
                             bearings2_vec, pairs_vec.cbegin(),
                             pairs_vec.cend());

  EssentialSolverResult<Iterator> essential_solver_result =
      opt_essential_solver_result.value_or(EssentialSolverResult<Iterator>());

  PyEssentialSolverResult py_ret;
  py_ret.From_(essential_solver_result, pairs_vec);

  return {py_ret, ransac_result};
}

void WrapEssentialSolver(nanobind::module_& module) {
  namespace nb = nanobind;
  using namespace nb::literals;  // NOLINT
                                 //
  // EssentialSolverPolicyFlags
  nb::enum_<EssentialSolverPolicyFlags> essential_solver_policy_flags(
      module, "EssentialSolverPolicyFlags", nb::is_arithmetic());
  essential_solver_policy_flags
      .value("ESSENTIAL_SOLVER_POLICY_FLAGS_NONE",
             EssentialSolverPolicyFlags::ESSENTIAL_SOLVER_POLICY_FLAGS_NONE)
      .value("ESSENTIAL_SOLVER_POLICY_FLAGS_CHECK_SINGULAR_VALUE",
             EssentialSolverPolicyFlags::
                 ESSENTIAL_SOLVER_POLICY_FLAGS_CHECK_SINGULAR_VALUE)
      .export_values();

  nb::class_<EssentialSolverPolicy> essential_solver_policy(
      module, "EssentialSolverPolicy");
  essential_solver_policy.def(nb::init<>())
      .def_rw("flags", &EssentialSolverPolicy::flags)
      .def_rw("ransac_policy", &EssentialSolverPolicy::ransac_policy)
      .def_rw("inlier_thr", &EssentialSolverPolicy::inlier_thr)
      .def_rw("check_singular_value_thr ",
              &EssentialSolverPolicy::check_singular_value_thr);

  nb::class_<PyEssentialSolverResult> essential_solver_result(
      module, "EssentialSolverResult");

  essential_solver_result.def(nb::init<>())
      .def_rw("rotation", &PyEssentialSolverResult::rotation)
      .def_rw("translation", &PyEssentialSolverResult::translation)
      .def_rw("inliers", &PyEssentialSolverResult::inliers)
      .def_rw("inlier_ratio", &PyEssentialSolverResult::inlier_ratio);

  module.def("compute_essential_matrix", &PyComputeEssentialMatrix,
             "essential_solver_policy"_a, "bearings1"_a, "bearings2"_a,
             "pairs"_a,
             "This function compute essentail matrix and execute RANSAC");
}

}  // namespace gsrap
