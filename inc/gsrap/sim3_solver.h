#ifndef GSRAP_SIM3_SOVLER_H
#define GSRAP_SIM3_SOVLER_H

#include <optional>
#include <unordered_map>
#include <vector>

#include "gsrap/gsrap-def.h"
#include "gsrap/macros.h"

GSRAP_IGNORE_STRICT_WARNING_PUSH
#include "Eigen/Core"
GSRAP_IGNORE_STRICT_WARNING_POP

namespace gsrap {
/**
 * @enum
 * The flags for Sim3SolverPolicy
 */
enum Sim3SolverPolicyFlags {
  //! All flags are false,
  SIM3_SOLVER_POLICY_FLAGS_NONE = 0,
  //! Refine model by using all inliers after RANSAC
  SIM3_SOLVER_POLICY_FLAGS_REFINE = (1 << 0),
};

/// Struct that defines solver runtime policy.
struct Sim3SolverPolicy {
  /// Flag
  int flags = SIM3_SOLVER_POLICY_FLAGS_NONE;
  /// RANSAC policy
  RansacPolicy ransac_policy;

  /// Threshold of inlier check
  double inlier_thr = 1e-12;

  size_t num_sample = 3;
};

/// Struct that is packed Sim3 solver result
template <typename Iterator>
struct Sim3SolverResult {
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  /// Rotation matrix
  Eigen::Matrix3d rotation;
  /// Translation vector
  Eigen::Vector3d translation;

  /// scale (scalar)
  double scale = 1.0;

  /// Iterator to inliers
  std::vector<Iterator> inliers;
  /// Inlier rate
  double inlier_ratio = 0.0;
};

/**
 * @fn
 *
 * @brief This function computes Sim3 transformation and execute RANSAC
 * @param sim3_solver_policy Sim3 solver policy
 * @param points1 points on coordinates 1
 * @param points2 points on coordinates 2
 * @param first Interator to the pair of indices of the points
 * @param last Interator to the pair of indices of the points
 * @return Pair of Sim3 solver result and report of RANSAC
 * @details TODO
 */
template <typename Iterator>
std::pair<std::optional<Sim3SolverResult<Iterator>>, RansacReport>
ComputeSim3Transformation(
    Sim3SolverPolicy sim3_solver_policy,
    const std::vector<Eigen::Vector3d
#ifdef GSRAP_LESS_EIGEN_3_4
                      ,
                      Eigen::aligned_allocator<Eigen::Vector3d>
#endif  // GSRAP_LESS_EIGEN_3_4
                      >& points1,
    const std::vector<Eigen::Vector3d
#ifdef GSRAP_LESS_EIGEN_3_4
                      ,
                      Eigen::aligned_allocator<Eigen::Vector3d>
#endif  // GSRAP_LESS_EIGEN_3_4
                      >& points2,
    Iterator first, Iterator last);

extern template std::pair<std::optional<Sim3SolverResult<std::vector<
                              std::pair<uint32_t, uint32_t>>::const_iterator>>,
                          RansacReport>
ComputeSim3Transformation(
    Sim3SolverPolicy sim3_solver_policy,
    const std::vector<Eigen::Vector3d
#ifdef GSRAP_LESS_EIGEN_3_4
                      ,
                      Eigen::aligned_allocator<Eigen::Vector3d>
#endif  // GSRAP_LESS_EIGEN_3_4
                      >& point1,
    const std::vector<Eigen::Vector3d
#ifdef GSRAP_LESS_EIGEN_3_4
                      ,
                      Eigen::aligned_allocator<Eigen::Vector3d>
#endif  // GSRAP_LESS_EIGEN_3_4
                      >& point2,
    std::vector<std::pair<uint32_t, uint32_t>>::const_iterator first,
    std::vector<std::pair<uint32_t, uint32_t>>::const_iterator last);

extern template std::pair<std::optional<Sim3SolverResult<std::unordered_map<
                              uint32_t, uint32_t>::const_iterator>>,
                          RansacReport>
ComputeSim3Transformation(
    Sim3SolverPolicy sim3_solver_policy,
    const std::vector<Eigen::Vector3d

#ifdef GSRAP_LESS_EIGEN_3_4
                      ,
                      Eigen::aligned_allocator<Eigen::Vector3d>
#endif  // GSRAP_LESS_EIGEN_3_4

                      >& points1,
    const std::vector<Eigen::Vector3d
#ifdef GSRAP_LESS_EIGEN_3_4
                      ,
                      Eigen::aligned_allocator<Eigen::Vector3d>
#endif  // GSRAP_LESS_EIGEN_3_4
                      >& points2,
    std::unordered_map<uint32_t, uint32_t>::const_iterator first,
    std::unordered_map<uint32_t, uint32_t>::const_iterator last);

}  // namespace gsrap

#endif  // GSRAP_SIM3_SOVLER_H
