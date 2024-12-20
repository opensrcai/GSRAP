#ifndef GSRAP_GSRAP_DEF_H_
#define GSRAP_GSRAP_DEF_H_

#include <stdint.h>
#include <variant>

namespace gsrap {

/**
 * @enum RansacTerminationInfo
 * reason for the termination of RANSAC
 */
enum class RansacTerminationInfo {
  //! reached upper limit
  REACHED_UPPER_LIMIT = 0,

  //! early stop
  EARLY_STOPED
};

/**
 * @enum RansacPolicyFlags
 * The flags for RansacPolicy
 */
enum RansacPolicyFlags {
  //! All flags are false.
  RANSAC_POLICY_FLAGS_NONE = 0,
  //! Enable early stopping
  RANSAC_POLICY_FLAGS_EARLY_STOP = 1 << 0,
  //! Use probability without duplication sample
  RANSAC_POLICY_FLAGS_USE_PROBABILITY_WITHOUT_DUPLICATION_SAMPLE = 1 << 1,
};

// Struct that defines RANSAC runtime policy.
struct RansacPolicy {
  //!  Flag
  int flags = RANSAC_POLICY_FLAGS_EARLY_STOP |
              RANSAC_POLICY_FLAGS_USE_PROBABILITY_WITHOUT_DUPLICATION_SAMPLE;
  //! The number of RANSAC iteration's lower limit.
  uint32_t num_ransac_itr_lower_limit = 10;
  //! The number of RANSAC iteration's upper limit.
  uint32_t num_ransac_itr_upper_limit = 100;

  //! The probability that the model could be estimated with only inliers.
  double probability = 0.99;

  //! The number of cpu threads which is used.
  uint32_t num_threads = 1;
};

// Struct that reports about RANSAC
struct RansacReport {
  //! The number of RANSAC iterations.
  uint32_t num_iteration;
  //! Enum of how RANSAC was terminated.
  RansacTerminationInfo ransac_termination_info;
};

struct PnpInlierCheckParamsUsingBearingVector {
  /// Threshold of inlier check(error angle)
  double inlier_thr = 1e-15;  // rad
};

struct PnpInlierCheckParamsUsingProjectedPoint {
  double fx = 0;
  double fy = 0;
  /// Threshold of inlier check(error px)
  double inlier_thr = 1.0;  // px
};

using PnpInlierCheckParams =
    std::variant<PnpInlierCheckParamsUsingBearingVector,
                 PnpInlierCheckParamsUsingProjectedPoint>;

}  // namespace gsrap

#endif  // GSRAP_GSRAP_DEF_H_
