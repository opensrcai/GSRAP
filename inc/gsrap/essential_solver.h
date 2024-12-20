#ifndef GSRAP_ESSENTIAL_SOVLER_H
#define GSRAP_ESSENTIAL_SOVLER_H

#include <stdint.h>

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
 * The flags for EssentialSolverPolicy
 */
enum EssentialSolverPolicyFlags {
  //! All flags are false,
  ESSENTIAL_SOLVER_POLICY_FLAGS_NONE = 0,
  /// If this flag is true, we check the singular value of the essential matrix
  ESSENTIAL_SOLVER_POLICY_FLAGS_CHECK_SINGULAR_VALUE = 1 << 0,
};

/// Struct that defines solver runtime policy.
struct EssentialSolverPolicy {
  /// Flag
  int flags = ESSENTIAL_SOLVER_POLICY_FLAGS_CHECK_SINGULAR_VALUE;
  /// RANSAC policy
  RansacPolicy ransac_policy;

  /// Threshold of inlier check(epipolar constraint)
  double inlier_thr = 1e-15;
  /// Threshold of checking singular value
  double check_singular_value_thr = 1e-15;
};

/// Struct that is packed essential solver result
template <typename Iterator>
struct EssentialSolverResult {
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  /// Essential matrix
  Eigen::Matrix3d essential;
  /// Rotation matrix
  Eigen::Matrix3d rotation;
  /// Translation matrix
  Eigen::Vector3d translation;

  /// Iterator to inliers
  std::vector<Iterator> inliers;
  /// Inlier rate
  double inlier_ratio = 0.0;
};

/**
 * @fn
 *
 * @brief This function computes essential matrix and execute RANSAC
 * @param essential_solver_policy Essential solver policy
 * @param bearings1 bearing vectors of camera 1
 * @param bearings2 bearing vectors of camera 2
 * @param first Interator to the pair of indices of the bearing vectors
 * @param last Interator to the pair of indices of the bearing vectors
 * @return Pair of essential solver result and report of RANSAC
 * @details TODO
 */
template <typename Iterator>
std::pair<std::optional<EssentialSolverResult<Iterator>>, RansacReport>
ComputeEssentialMatrix(
    EssentialSolverPolicy essential_solver_policy,
    const std::vector<Eigen::Vector3d

#ifdef GSRAP_LESS_EIGEN_3_4
                      ,
                      Eigen::aligned_allocator<Eigen::Vector3d>
#endif  // GSRAP_LESS_EIGEN_3_4
                      >& bearings1,
    const std::vector<Eigen::Vector3d
#ifdef GSRAP_LESS_EIGEN_3_4
                      ,
                      Eigen::aligned_allocator<Eigen::Vector3d>
#endif  // GSRAP_LESS_EIGEN_3_4
                      >& bearings2,
    Iterator first, Iterator last);

extern template std::pair<std::optional<EssentialSolverResult<std::vector<
                              std::pair<uint32_t, uint32_t>>::const_iterator>>,
                          RansacReport>
ComputeEssentialMatrix(
    EssentialSolverPolicy essential_solver_policy,
    const std::vector<Eigen::Vector3d
#ifdef GSRAP_LESS_EIGEN_3_4
                      ,
                      Eigen::aligned_allocator<Eigen::Vector3d>
#endif  // GSRAP_LESS_EIGEN_3_4
                      >& bearings1,
    const std::vector<Eigen::Vector3d
#ifdef GSRAP_LESS_EIGEN_3_4
                      ,
                      Eigen::aligned_allocator<Eigen::Vector3d>
#endif  // GSRAP_LESS_EIGEN_3_4
                      >& bearings2,
    std::vector<std::pair<uint32_t, uint32_t>>::const_iterator first,
    std::vector<std::pair<uint32_t, uint32_t>>::const_iterator last);

extern template std::pair<
    std::optional<EssentialSolverResult<
        std::unordered_map<uint32_t, uint32_t>::const_iterator>>,
    RansacReport>
ComputeEssentialMatrix(
    EssentialSolverPolicy essential_solver_policy,
    const std::vector<Eigen::Vector3d
#ifdef GSRAP_LESS_EIGEN_3_4
                      ,
                      Eigen::aligned_allocator<Eigen::Vector3d>
#endif  // GSRAP_LESS_EIGEN_3_4
                      >& bearings1,
    const std::vector<Eigen::Vector3d
#ifdef GSRAP_LESS_EIGEN_3_4
                      ,
                      Eigen::aligned_allocator<Eigen::Vector3d>
#endif  // GSRAP_LESS_EIGEN_3_4
                      >& bearings2,
    std::unordered_map<uint32_t, uint32_t>::const_iterator first,
    std::unordered_map<uint32_t, uint32_t>::const_iterator last);

}  // namespace gsrap

#endif  // GSRAP_ESSENTIAL_SOVLER_H
