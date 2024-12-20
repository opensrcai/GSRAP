#ifndef GSRAP_P3P_SOVLER_H
#define GSRAP_P3P_SOVLER_H

#include <optional>
#include <unordered_map>
#include <variant>
#include <vector>

#include "gsrap/gsrap-def.h"
#include "gsrap/macros.h"

GSRAP_IGNORE_STRICT_WARNING_PUSH
#include "Eigen/Core"
GSRAP_IGNORE_STRICT_WARNING_POP

namespace gsrap {

/**
 * @enum
 * The flags for P3PSolverPolicy
 */
enum P3pSolverPolicyFlags {
  //! All flags are false,
  P3P_SOLVER_POLICY_FLAGS_NONE = 0,
};

/**
 * @enum
 * The algorithm used.
 */
enum class P3pSolverAlgorithm {
  // Tong Ke, Stergios I. Roumeliotis, "An Efficient Algebraic Solution to the
  // Perspective-Three-Point Problem", CVPR2017
  // https://openaccess.thecvf.com/content_cvpr_2017/html/Ke_An_Efficient_Algebraic_CVPR_2017_paper.html
  KE_2017 = 0,
};

/// Struct that defines solver runtime policy.
struct P3pSolverPolicy {
  /// Flag
  int flags = P3P_SOLVER_POLICY_FLAGS_NONE;
  /// Algorithm
  P3pSolverAlgorithm algorithm = P3pSolverAlgorithm::KE_2017;
  /// RANSAC policy
  RansacPolicy ransac_policy;

  PnpInlierCheckParams pnp_inlier_check_params =
      PnpInlierCheckParamsUsingBearingVector{};
};

/// Struct that is packed P3P solver result
template <typename Iterator>
struct P3pSolverResult {
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
 * @param p3p_solver_policy PnP solver policy
 * @param bearings bearing vectors of camera
 * @param points points on world coordinates
 * @param first Interator to the pair of indices of a bearing vectors and a point
 * @param last Interator to the pair of indices of a bearing vector and a point
 * @return Pair of PnP solver result and report of RANSAC
 * @details TODO
 */
// clang-format on
template <typename Iterator>
std::pair<std::optional<P3pSolverResult<Iterator>>, RansacReport>
SolveP3pProblem(P3pSolverPolicy p3p_solver_policy,
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

extern template std::pair<std::optional<P3pSolverResult<std::vector<
                              std::pair<uint32_t, uint32_t>>::const_iterator>>,
                          RansacReport>
SolveP3pProblem(
    P3pSolverPolicy p3p_solver_policy,
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

extern template std::pair<std::optional<P3pSolverResult<std::unordered_map<
                              uint32_t, uint32_t>::const_iterator>>,
                          RansacReport>
SolveP3pProblem(P3pSolverPolicy p3p_solver_policy,
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
std::
    vector<std::pair<Eigen::Matrix3d, Eigen::Vector3d>

#ifdef GSRAP_LESS_EIGEN_3_4
           ,
           Eigen::aligned_allocator<std::pair<Eigen::Matrix3d, Eigen::Vector3d>>
#endif  // GSRAP_LESS_EIGEN_3_4
           >
    SolveP3pProblemWithoutRansac(
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

#endif  // GSRAP_P3P_SOVLER_H
