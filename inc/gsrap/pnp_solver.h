#ifndef GSRAP_PNP_SOVLER_H
#define GSRAP_PNP_SOVLER_H

#include <memory>
#include <optional>
#include <unordered_map>
#include <utility>
#include <vector>

#include "gsrap/gsrap-def.h"
#include "gsrap/macros.h"

GSRAP_IGNORE_STRICT_WARNING_PUSH
#include "Eigen/Core"
GSRAP_IGNORE_STRICT_WARNING_POP

namespace gsrap {

/**
 * @enum
 * The flags for PnPSolverPolicy
 */
enum PnpSolverPolicyFlags {
  //! All flags are false,
  PNP_SOLVER_POLICY_FLAGS_NONE = 0,
  //! Refine model by using all inliers after RANSAC
  PNP_SOLVER_POLICY_FLAGS_REFINE = (1 << 0)
};

/// Struct that defines solver runtime policy.
struct PnpSolverPolicy {
  /// Flag
  int flags = PNP_SOLVER_POLICY_FLAGS_NONE;
  /// RANSAC policy
  RansacPolicy ransac_policy;

  PnpInlierCheckParams pnp_inlier_check_params =
      PnpInlierCheckParamsUsingBearingVector{};
};

/// Struct that is packed PnP solver result
template <typename Iterator>
struct PnpSolverResult {
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  /// Rotation matrix
  Eigen::Matrix3d rotation;
  /// Translation matrix
  Eigen::Vector3d translation;

  /// Iterator to inliers
  std::vector<Iterator> inliers;
  /// Inlier rate
  double inlier_ratio = 0.0;
};

// clang-format off
/**
 * @fn
 *
 * @brief This function solves PnP problem and execute RANSAC
 * @param pnp_solver_policy PnP solver policy
 * @param bearings bearing vectors of camera
 * @param points points on world coordinates
 * @param first Interator to the pair of indices of a bearing vectors and a point
 * @param last Interator to the pair of indices of a bearing vector and a point
 * @return Pair of PnP solver result and report of RANSAC
 * @details TODO
 */
// clang-format on
template <typename Iterator>
std::pair<std::optional<PnpSolverResult<Iterator>>, RansacReport>
SolvePnpProblem(PnpSolverPolicy pnp_solver_policy,
                const std::vector<Eigen::Vector3d
#ifdef GSRAP_LESS_EIGEN_3_4
                                  ,
                                  Eigen::aligned_allocator<Eigen::Vector3d>
#endif  // GSRAP_LESS_EIGEN_3_4
                                  >& bearings,
                const std::vector<Eigen::Vector3d
#ifdef GSRAP_LESS_EIGEN_3_4
                                  ,
                                  Eigen::aligned_allocator<Eigen::Vector3d>
#endif  // GSRAP_LESS_EIGEN_3_4
                                  >& points,
                Iterator first, Iterator last);

extern template std::pair<std::optional<PnpSolverResult<std::vector<
                              std::pair<uint32_t, uint32_t>>::const_iterator>>,
                          RansacReport>
SolvePnpProblem(
    PnpSolverPolicy pnp_solver_policy,
    const std::vector<Eigen::Vector3d
#ifdef GSRAP_LESS_EIGEN_3_4
                      ,
                      Eigen::aligned_allocator<Eigen::Vector3d>
#endif  // GSRAP_LESS_EIGEN_3_4
                      >& bearings,
    const std::vector<Eigen::Vector3d
#ifdef GSRAP_LESS_EIGEN_3_4
                      ,
                      Eigen::aligned_allocator<Eigen::Vector3d>
#endif  // GSRAP_LESS_EIGEN_3_4
                      >& points,
    std::vector<std::pair<uint32_t, uint32_t>>::const_iterator first,
    std::vector<std::pair<uint32_t, uint32_t>>::const_iterator last);

extern template std::pair<std::optional<PnpSolverResult<std::unordered_map<
                              uint32_t, uint32_t>::const_iterator>>,
                          RansacReport>
SolvePnpProblem(PnpSolverPolicy pnp_solver_policy,
                const std::vector<Eigen::Vector3d
#ifdef GSRAP_LESS_EIGEN_3_4
                                  ,
                                  Eigen::aligned_allocator<Eigen::Vector3d>
#endif  // GSRAP_LESS_EIGEN_3_4
                                  >& bearings,
                const std::vector<Eigen::Vector3d
#ifdef GSRAP_LESS_EIGEN_3_4
                                  ,
                                  Eigen::aligned_allocator<Eigen::Vector3d>
#endif  // GSRAP_LESS_EIGEN_3_4
                                  >& points,
                std::unordered_map<uint32_t, uint32_t>::const_iterator first,
                std::unordered_map<uint32_t, uint32_t>::const_iterator last);

/**
 * @param bearings bearing vectors of camera
 * @param points points on world coordinates
 */
std::pair<Eigen::Matrix3d, Eigen::Vector3d> SolvePnpProblemWithoutRansac(
    const std::vector<Eigen::Vector3d
#ifdef GSRAP_LESS_EIGEN_3_4
                      ,
                      Eigen::aligned_allocator<Eigen::Vector3d>
#endif  // GSRAP_LESS_EIGEN_3_4
                      >& bearings,
    const std::vector<Eigen::Vector3d
#ifdef GSRAP_LESS_EIGEN_3_4
                      ,
                      Eigen::aligned_allocator<Eigen::Vector3d>
#endif  // GSRAP_LESS_EIGEN_3_4
                      >& points);

}  // namespace gsrap

#endif  // GSRAP_PNP_SOVLER_H
