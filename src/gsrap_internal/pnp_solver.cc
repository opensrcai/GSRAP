#include "gsrap/pnp_solver.h"

#include <assert.h>

#include <array>
#include <memory>
#include <unordered_set>
#include <utility>

#include "gsrap/macros.h"
#include "gsrap_internal/common.h"
#include "gsrap_internal/ransac.h"
#include "gsrap_internal/sample.h"

GSRAP_IGNORE_STRICT_WARNING_PUSH

#include "opengv/absolute_pose/modules/Epnp.hpp"

GSRAP_IGNORE_STRICT_WARNING_POP

namespace gsrap {

// This function checks PnpSolverPolicy.
static PnpSolverPolicy CheckPnpSolverPolicy(PnpSolverPolicy pnp_solver_policy) {
  PnpSolverPolicy ret = pnp_solver_policy;

  if (auto *ptr_hcc = std::get_if<PnpInlierCheckParamsUsingBearingVector>(
          &ret.pnp_inlier_check_params)) {
    // 'inlier_thr' should be greater than zero.
    ptr_hcc->inlier_thr = std::max(0.0, ptr_hcc->inlier_thr);
  } else if (auto *ptr_pc =
                 std::get_if<PnpInlierCheckParamsUsingProjectedPoint>(
                     &ret.pnp_inlier_check_params)) {
    // 'inlier_thr' should be greater than zero.
    ptr_pc->inlier_thr = std::max(0.0, ptr_pc->inlier_thr);
  }
  return ret;
}

template <typename Iterator>
void SolveWithEpnp(const eigen_vector<Eigen::Vector3d>& bearings,
                   const eigen_vector<Eigen::Vector3d>& points,
                   const std::vector<Iterator>& data, Eigen::Matrix3d* R,
                   Eigen::Vector3d* t) {
  const int num_correspondances = int(data.size());

  opengv::absolute_pose::modules::Epnp _epnp_solver = {};

  _epnp_solver.set_maximum_number_of_correspondences(num_correspondances);
  for (const auto& _data : data) {
    assert(0 <= _data->first && _data->first < bearings.size());
    assert(0 <= _data->second && _data->second < points.size());

    const Eigen::Vector3d& bearing = bearings[_data->first];
    const Eigen::Vector3d& point   = points[_data->second];

    _epnp_solver.add_correspondence(point, bearing);
  }
  _epnp_solver.compute_pose(*R, *t);
}

/// Function object to compute  rotation matrix, and translation vector.
template <typename Iterator>
struct Epnp {
  std::vector<std::unique_ptr<Rt>> operator()(
      const eigen_vector<Eigen::Vector3d>& bearings,
      const eigen_vector<Eigen::Vector3d>& points,
      const std::vector<Iterator>& data) const {
    assert(data.size() >= 4);
    if (data.size() < 4) {
      return std::vector<std::unique_ptr<Rt>>();
    }

    const int num_correspondances = int(data.size());

    Mat3d R;
    Vec3d t;
    SolveWithEpnp(bearings, points, data, &R, &t);

    std::vector<std::unique_ptr<Rt>> ret;
    ret.emplace_back(std::make_unique<Rt>(R, t));

    return std::vector<std::unique_ptr<Rt>>(std::move(ret));
  }
};

template <typename Iterator>
std::pair<std::optional<PnpSolverResult<Iterator>>, RansacReport>
SolvePnpProblem(PnpSolverPolicy pnp_solver_policy,
                const eigen_vector<Vec3d>& bearings,
                const eigen_vector<Vec3d>& points, Iterator first,
                Iterator last) {
  pnp_solver_policy = CheckPnpSolverPolicy(pnp_solver_policy);

  Epnp<Iterator> pnp_solver;
  Sample<Iterator> sample;

  const auto num_data = size_t(std::distance(first, last));

  auto result = std::visit(
      [&](const auto& params) {
        IsInlierForPnp<decltype(params), Iterator> is_inlier(params);

        // Execute RANSAC.
        return Ransac(
            pnp_solver_policy.ransac_policy, first, last,
            std::bind(pnp_solver, bearings, points, std::placeholders::_1),
            std::bind(is_inlier, bearings, points, std::placeholders::_1,
                      std::placeholders::_2),
            std::bind(sample, num_data, 4, std::placeholders::_1,
                      std::placeholders::_2, std::placeholders::_3));
      },
      pnp_solver_policy.pnp_inlier_check_params);

  // Refine model by using all inliers after RANSAC
  const bool do_refine =
      result.first && result.first->num_inliers >= 4 &&
      (pnp_solver_policy.flags & PNP_SOLVER_POLICY_FLAGS_REFINE);

  if (do_refine) {
    const std::vector<std::unique_ptr<Rt>> resolved =
        pnp_solver(bearings, points, result.first->inliers);

    if (!resolved.empty()) {
      // Pack the refine results.
      PnpSolverResult<Iterator> ret;

      ret.rotation    = resolved[0]->R;
      ret.translation = resolved[0]->t;

      ret.inliers      = std::move(result.first->inliers);
      ret.inlier_ratio = result.first->inlier_ratio;

      return std::make_pair(ret, result.second);
    }
  }

  // Refinement is not required or failed to that, return the best model of
  // RANSAC Return the best model when it available
  if (result.first) {
    PnpSolverResult<Iterator> ret;
    ret.rotation    = result.first.value().model->R;
    ret.translation = result.first.value().model->t;

    ret.inliers      = std::move(result.first->inliers);
    ret.inlier_ratio = result.first->inlier_ratio;

    return std::make_pair(ret, result.second);
  }

  // Or return empty result
  return std::make_pair(std::nullopt, result.second);
}

template std::pair<std::optional<PnpSolverResult<std::vector<
                       std::pair<uint32_t, uint32_t>>::const_iterator>>,
                   RansacReport>
SolvePnpProblem(
    PnpSolverPolicy pnp_solver_policy, const eigen_vector<Vec3d>& bearings,
    const eigen_vector<Vec3d>& points,
    std::vector<std::pair<uint32_t, uint32_t>>::const_iterator first,
    std::vector<std::pair<uint32_t, uint32_t>>::const_iterator last);

template std::pair<std::optional<PnpSolverResult<
                       std::unordered_map<uint32_t, uint32_t>::const_iterator>>,
                   RansacReport>
SolvePnpProblem(PnpSolverPolicy pnp_solver_policy,
                const eigen_vector<Vec3d>& bearings,
                const eigen_vector<Vec3d>& points,
                std::unordered_map<uint32_t, uint32_t>::const_iterator first,
                std::unordered_map<uint32_t, uint32_t>::const_iterator last);

std::pair<Mat3d, Vec3d> SolvePnpProblemWithoutRansac(
    const eigen_vector<Vec3d>& bearings, const eigen_vector<Vec3d>& points) {
  const size_t num = std::min(bearings.size(), points.size());

  assert(bearings.size() == points.size());
  assert(num >= 4);
  if (num < 4) {
    return std::pair<Eigen::Matrix3d, Eigen::Vector3d>();
  }

  using Iterator = std::vector<std::pair<uint32_t, uint32_t>>::const_iterator;
  std::vector<std::pair<uint32_t, uint32_t>> pairs;
  std::vector<Iterator> data;
  pairs.reserve(num);
  data.reserve(num);
  for (size_t i = 0; i < num; ++i) {
    pairs.emplace_back(i, i);
    data.emplace_back(pairs.end() - 1);
  }

  Mat3d R;
  Vec3d t;
  SolveWithEpnp(bearings, points, data, &R, &t);

  return {R, t};
}

}  // namespace gsrap
