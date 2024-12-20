#include "pnp_solver.h"

#include "common.h"
#include "gsrap/pnp_solver.h"

GSRAP_IGNORE_STRICT_WARNING_PUSH
#include "nanobind/eigen/dense.h"
#include "nanobind/stl/optional.h"
#include "nanobind/stl/pair.h"
#include "nanobind/stl/vector.h"
#include "nanobind/stl/variant.h"
GSRAP_IGNORE_STRICT_WARNING_POP
namespace gsrap {

using Pair     = std::pair<uint32_t, uint32_t>;
using Iterator = std::vector<Pair>::const_iterator;

template <typename Iterator, typename Container>
void PyPnpSolverResult::From_(const PnpSolverResult<Iterator>& result,
                              const Container& pairs) {
  rotation    = result.rotation;
  translation = result.translation;

  const auto& _inliers = result.inliers;
  inliers.resize(_inliers.size());
  std::transform(
      _inliers.cbegin(), _inliers.cend(), inliers.begin(),
      [&](Iterator itr) { return static_cast<uint32_t>(itr - pairs.begin()); });

  inlier_ratio = result.inlier_ratio;
}

std::pair<std::optional<PyPnpSolverResult>, RansacReport> PySolvePnpProblem(
    PnpSolverPolicy pnp_solver_policy,
    const Eigen::Matrix<double, -1, 3, Eigen::RowMajorBit>& bearings,
    const Eigen::Matrix<double, -1, 3, Eigen::RowMajorBit>& points,
    const RMatrixXu32& pairs) {
  eigen_vector<Vec3d> bearings_vec(static_cast<size_t>(bearings.rows()));
  eigen_vector<Vec3d> points_vec(static_cast<size_t>(points.rows()));

  for (int64_t i = 0; i < bearings.rows(); ++i) {
    bearings_vec[static_cast<size_t>(i)] = bearings.row(i);
  }
  for (int64_t i = 0; i < points.rows(); ++i) {
    points_vec[static_cast<size_t>(i)] = points.row(i);
  }

  std::vector<Pair> pairs_vec(static_cast<size_t>(pairs.rows()));

  for (int64_t i = 0; i < pairs.rows(); ++i) {
    pairs_vec[static_cast<size_t>(i)].first  = pairs.row(i).x();
    pairs_vec[static_cast<size_t>(i)].second = pairs.row(i).y();
  }

  const auto [opt_pnp_solver_result, ransac_result] =
      SolvePnpProblem(pnp_solver_policy, bearings_vec, points_vec,
                      pairs_vec.cbegin(), pairs_vec.cend());

  if (!opt_pnp_solver_result.has_value()) {
    return {std::nullopt, ransac_result};
  }
  PnpSolverResult<Iterator> pnp_solver_result =
      opt_pnp_solver_result.value_or(PnpSolverResult<Iterator>());

  PyPnpSolverResult py_ret;
  py_ret.From_(pnp_solver_result, pairs_vec);

  return {py_ret, ransac_result};
}

void WrapPnpSolver(nanobind::module_& module) {
  namespace nb = nanobind;
  using namespace nb::literals;  // NOLINT

  // PnpSolverPolicyFlags
  nb::enum_<PnpSolverPolicyFlags> pnp_solver_policy_flags(
      module, "PnpSolverPolicyFlags", nb::is_arithmetic());
  pnp_solver_policy_flags
      .value("PNP_SOLVER_POLICY_FLAGS_NONE",
             PnpSolverPolicyFlags::PNP_SOLVER_POLICY_FLAGS_NONE)
      .value("PNP_SOLVER_POLICY_FLAGS_REFINE",
             PnpSolverPolicyFlags::PNP_SOLVER_POLICY_FLAGS_REFINE)
      .export_values();

  nb::class_<PnpSolverPolicy> pnp_solver_policy(module, "PnpSolverPolicy");
  pnp_solver_policy.def(nb::init<>())
      .def_rw("flags", &PnpSolverPolicy::flags)
      .def_rw("ransac_policy", &PnpSolverPolicy::ransac_policy)
      .def_rw("pnp_inlier_check_params", &PnpSolverPolicy::pnp_inlier_check_params);

  nb::class_<PyPnpSolverResult> pnp_solver_result(module, "PnpSolverResult");
  pnp_solver_result.def(nb::init<>())
      .def_rw("rotation", &PyPnpSolverResult::rotation)
      .def_rw("translation", &PyPnpSolverResult::translation)
      .def_rw("inliers", &PyPnpSolverResult::inliers)
      .def_rw("inlier_ratio", &PyPnpSolverResult::inlier_ratio);

  module.def("solve_pnp_problem", &PySolvePnpProblem, "pnp_solver_policy"_a,
             "bearings"_a, "points"_a, "pairs"_a,
             "This function solve PnP prolblem and execute RANSAC");
}

}  // namespace gsrap
