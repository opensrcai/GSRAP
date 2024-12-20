// nanobind
#include <nanobind/eigen/dense.h>
#include <nanobind/nanobind.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/pair.h>
#include <nanobind/stl/vector.h>

#include "essential_solver.h"
#include "pnp_solver.h"
#include "sim3_solver.h"

// GSRAP
#include "gsrap/gsrap-def.h"
#include "gsrap/sim3_solver.h"

namespace nb = nanobind;

NB_MODULE(gsrap_ext, module) {   // NOLINT
  using namespace gsrap;         // NOLINT
  using namespace nb::literals;  // NOLINT

  module.doc() =
      "Geometric Solvers for Reconstruction And Pose estimation(GSRAP)";

  // RansacTerminationInfo
  nb::enum_<RansacTerminationInfo> ransac_termination_info(
      module, "RansacTerminationInfo", nb::is_arithmetic());

  ransac_termination_info
      .value("REACHED_UPPER_LIMIT", RansacTerminationInfo::REACHED_UPPER_LIMIT)
      .value("EARLY_STOPED", RansacTerminationInfo::EARLY_STOPED)
      .export_values();

  // RansacPolicyFlags
  nb::enum_<RansacPolicyFlags> ransac_policy_flags(module, "RansacPolicyFlags",
                                                   nb::is_arithmetic());
  ransac_policy_flags
      .value("RANSAC_POLICY_FLAGS_NONE",
             RansacPolicyFlags::RANSAC_POLICY_FLAGS_NONE)
      .value("RANSAC_POLICY_FLAGS_EARLY_STOP",
             RansacPolicyFlags::RANSAC_POLICY_FLAGS_EARLY_STOP)
      .value("RANSAC_POLICY_FLAGS_USE_PROBABILITY_WITHOUT_DUPLICATION_SAMPLE",
             RansacPolicyFlags::
                 RANSAC_POLICY_FLAGS_USE_PROBABILITY_WITHOUT_DUPLICATION_SAMPLE)
      .export_values();

  // RansacPolicy
  nb::class_<RansacPolicy> ransac_policy(module, "RansacPolicy");

  ransac_policy.def(nb::init<>())
      .def_rw("flags", &RansacPolicy::flags)
      .def_rw("num_ransac_itr_lower_limit",
              &RansacPolicy::num_ransac_itr_lower_limit)
      .def_rw("num_ransac_itr_upper_limit",
              &RansacPolicy::num_ransac_itr_upper_limit)
      .def_rw("probability", &RansacPolicy::probability)
      .def_rw("num_threads", &RansacPolicy::num_threads);

  // RansacReport
  nb::class_<RansacReport> ransac_report(module, "RansacReport");
  ransac_report.def(nb::init<>())
      .def_rw("num_iteration", &RansacReport::num_iteration)
      .def_rw("ransac_termination_info",
              &RansacReport::ransac_termination_info);

  // PnpInlierCheckParamsUsingBearingVector
  nb::class_<PnpInlierCheckParamsUsingBearingVector>
      pnp_inlier_check_params_using_bearing_vetor(
          module, "PnpInlierCheckParamsUsingBearingVector");
  pnp_inlier_check_params_using_bearing_vetor.def(nb::init<>())
      .def_rw("inlier_thr",
              &PnpInlierCheckParamsUsingBearingVector::inlier_thr);

  // PnpInlierCheckParamsUsingProjectedPoint
  nb::class_<PnpInlierCheckParamsUsingProjectedPoint>
      pnp_inlier_check_params_using_projected_point(
          module, "PnpInlierCheckParamsUsingProjectedPoint");
  pnp_inlier_check_params_using_projected_point.def(nb::init<>())
      .def_rw("inlier_thr",
              &PnpInlierCheckParamsUsingProjectedPoint::inlier_thr)
      .def_rw("fx", &PnpInlierCheckParamsUsingProjectedPoint::fx)
      .def_rw("fy", &PnpInlierCheckParamsUsingProjectedPoint::fy);

  // Sim3SolverPolicyFlags
  nb::enum_<Sim3SolverPolicyFlags> sim3_solver_policy_flags(
      module, "Sim3SolverPolicyFlags", nb::is_arithmetic());
  sim3_solver_policy_flags
      .value("SIM3_SOLVER_POLICY_FLAGS_NONE",
             Sim3SolverPolicyFlags::SIM3_SOLVER_POLICY_FLAGS_NONE)
      .value("SIM3_SOLVER_POLICY_FLAGS_REFINE",
             Sim3SolverPolicyFlags::SIM3_SOLVER_POLICY_FLAGS_REFINE)
      .export_values();

  // Sim3SolverPolicy
  nb::class_<Sim3SolverPolicy> sim3_solver_policy(module, "Sim3SolverPolicy");
  sim3_solver_policy.def(nb::init<>())
      .def_rw("flags", &Sim3SolverPolicy::flags)
      .def_rw("ransac_policy", &Sim3SolverPolicy::ransac_policy)
      .def_rw("inlier_thr", &Sim3SolverPolicy::inlier_thr)
      .def_rw("num_sample", &Sim3SolverPolicy::num_sample);

  // Sim3SolverResult
  nb::class_<PySim3SolverResult> sim3_solver_result(module, "Sim3SolverResult");
  sim3_solver_result.def(nb::init<>())
      .def_rw("rotation", &PySim3SolverResult::rotation)
      .def_rw("translation", &PySim3SolverResult::translation)
      .def_rw("scale", &PySim3SolverResult::scale)
      .def_rw("inliers", &PySim3SolverResult::inliers)
      .def_rw("inlier_ratio", &PySim3SolverResult::inlier_ratio);

  module.def("compute_sim3_transformation", &PyComputeSim3Transformation,
             "sim3_solver_policy"_a, "points1"_a, "point2"_a, "pairs"_a,
             "This function computes Sim3 transformation and execute RANSAC");

  WrapPnpSolver(module);

  WrapEssentialSolver(module);
}

#if 0
// pybind11
#include "pybind11/pybind11.h"
#include "sim3_solver.h"

namespace py = pybind11;

// GSRAP
#include "gsrap/gsrap-def.h"
#include "gsrap/sim3_solver.h"

PYBIND11_MODULE(gsrap_py, module) {    // NOLINT
  using namespace gsrap;               // NOLINT
  using namespace pybind11::literals;  // NOLINT

  // RansacTerminationInfo
  py::enum_<RansacTerminationInfo> ransac_termination_info(
      module, "RansacTerminationInfo", py::arithmetic());

  ransac_termination_info
      .value("REACHED_UPPER_LIMIT", RansacTerminationInfo::REACHED_UPPER_LIMIT)
      .value("EARLY_STOPED", RansacTerminationInfo::EARLY_STOPED)
      .export_values();

  // RansacPolicyFlags
  py::enum_<RansacPolicyFlags> ransac_policy_flags(module, "RansacPolicyFlags",
                                                   py::arithmetic());

  ransac_policy_flags
      .value("RANSAC_POLICY_FLAGS_NONE",
             RansacPolicyFlags::RANSAC_POLICY_FLAGS_NONE)
      .value("RANSAC_POLICY_FLAGS_EARLY_STOP",
             RansacPolicyFlags::RANSAC_POLICY_FLAGS_EARLY_STOP)
      .value("RANSAC_POLICY_FLAGS_USE_PROBABILITY_WITHOUT_DUPLICATION_SAMPLE",
             RansacPolicyFlags::
                 RANSAC_POLICY_FLAGS_USE_PROBABILITY_WITHOUT_DUPLICATION_SAMPLE)
      .export_values();

  // RansacPolicy
  py::class_<RansacPolicy> ransac_policy(module, "RansacPolicy");

  ransac_policy.def(py::init<>())
      .def_readwrite("flags", &RansacPolicy::flags)
      .def_readwrite("num_ransac_itr_lower_limit",
                     &RansacPolicy::num_ransac_itr_lower_limit)
      .def_readwrite("num_ransac_itr_upper_limit",
                     &RansacPolicy::num_ransac_itr_upper_limit)
      .def_readwrite("probability", &RansacPolicy::probability)
      .def_readwrite("num_threads", &RansacPolicy::num_threads);

  // RansacReport
  py::class_<RansacReport> ransac_report(module, "RansacReport");
  ransac_report.def(py::init<>())
      .def_readwrite("num_iteration", &RansacReport::num_iteration)
      .def_readwrite("ransac_termination_info",
                     &RansacReport::ransac_termination_info);

  // Sim3SolverPolicyFlags
  py::enum_<Sim3SolverPolicyFlags> sim3_solver_policy_flags(
      module, "Sim3SolverPolicyFlags", py::arithmetic());
  sim3_solver_policy_flags
      .value("SIM3_SOLVER_POLICY_FLAGS_NONE",
             Sim3SolverPolicyFlags::SIM3_SOLVER_POLICY_FLAGS_NONE)
      .value("SIM3_SOLVER_POLICY_FLAGS_REFINE",
             Sim3SolverPolicyFlags::SIM3_SOLVER_POLICY_FLAGS_REFINE)
      .export_values();

  // Sim3SolverPolicy
  py::class_<Sim3SolverPolicy> sim3_solver_policy(module, "Sim3SolverPolicy");
  sim3_solver_policy.def(py::init<>())
      .def_readwrite("flags", &Sim3SolverPolicy::flags)
      .def_readwrite("ransac_policy", &Sim3SolverPolicy::ransac_policy)
      .def_readwrite("inlier_thr", &Sim3SolverPolicy::inlier_thr)
      .def_readwrite("num_sample", &Sim3SolverPolicy::num_sample);

  // Sim3SolverResult
  py::class_<PySim3SolverResult> sim3_solver_result(module, "Sim3SolverResult");
  sim3_solver_result.def(py::init<>())
      .def_readwrite("rotation", &PySim3SolverResult::rotation)
      .def_readwrite("translation", &PySim3SolverResult::translation)
      .def_readwrite("scale", &PySim3SolverResult::scale)
      .def_readwrite("inliers", &PySim3SolverResult::inliers)
      .def_readwrite("inlier_ratio", &PySim3SolverResult::inlier_ratio);

  module.def("compute_sim3_transformation", &PyComputeSim3Transformation,
             "sim3_solver_policy"_a, "points1"_a, "point2"_a, "pairs"_a,
             "This function computes Sim3 transformation and execute RANSAC");
}
#endif
