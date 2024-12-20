#ifndef GSRAP_INTERNAL_COMMON_H
#define GSRAP_INTERNAL_COMMON_H

#include "gsrap/gsrap-def.h"
#include "gsrap/macros.h"

GSRAP_IGNORE_STRICT_WARNING_PUSH
#include <Eigen/Core>
GSRAP_IGNORE_STRICT_WARNING_POP

namespace gsrap {

using Mat3d = Eigen::Matrix3d;
using Vec3d = Eigen::Vector3d;

template <typename T>
using eigen_vector = std::vector<T
#ifdef GSRAP_LESS_EIGEN_3_4
                                 ,
                                 Eigen::aligned_allocator<T>
#endif  // GSRAP_LESS_EIGEN_3_4
                                 >;

// Model
struct Rt {
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  inline Rt() = default;
  inline Rt(Mat3d R_, Vec3d t_) : R(std::move(R_)), t(std::move(t_)) {}

  Mat3d R;  // Rotation matrix
  Vec3d t;  // Translation vector
};

/*
/// Function object to determine if a pair of bearing vectors and points is
/// inlier.
template <typename Iterator>
struct IsInlierForPnp {
  bool operator()(const eigen_vector<Eigen::Vector3d>& bearings,
                  const eigen_vector<Eigen::Vector3d>& points,
                  const typename Iterator::value_type& match, const Rt& rt,
                  const double thr = 1e-15) const {
    assert(0 <= match.first && match.first < bearings.size());
    assert(0 <= match.second && match.second < points.size());

    const Vec3d& bearing = bearings[match.first];
    const Vec3d& point   = points[match.second];

    const Vec3d bearing_ = (rt.R * point + rt.t).normalized();

    return acos(std::clamp(bearing.dot(bearing_), -1.0, 1.0)) <
           thr;  // Epipolar constraint
    // Note: acos(x) <- [0, pi]
  }
};
*/

template <typename T>
constexpr bool false_v = false;

/// Function object to determine if a pair of bearing vectors and points is
/// inlier.
template <typename Params, typename Iterator>
struct IsInlierForPnp {
  IsInlierForPnp() = delete;
  explicit IsInlierForPnp(const Params &params) noexcept
      : inlier_check_params_(params) {}

  bool operator()(const eigen_vector<Eigen::Vector3d> &bearings,
                  const eigen_vector<Eigen::Vector3d> &points,
                  const typename Iterator::value_type &match,
                  const Rt &rt) const noexcept {
    assert(0 <= match.first && match.first < bearings.size());
    assert(0 <= match.second && match.second < points.size());

    const Vec3d &bearing = bearings[match.first];
    const Vec3d &point   = points[match.second];

    if constexpr (std::is_same<const PnpInlierCheckParamsUsingProjectedPoint &,
                               Params>::value) {
      const auto *const params =
          std::get_if<PnpInlierCheckParamsUsingProjectedPoint>(
              &inlier_check_params_);

      const double u = params->fx * bearing.x() / bearing.z();
      const double v = params->fy * bearing.y() / bearing.z();

      const Vec3d point_c = rt.R * point + rt.t;
      const double u_     = params->fx * point_c.x() / point_c.z();
      const double v_     = params->fy * point_c.y() / point_c.z();

      const double du = u - u_;
      const double dv = v - v_;

      return du * du + dv * dv < params->inlier_thr * params->inlier_thr;
    } else if constexpr (std::is_same<
                             const PnpInlierCheckParamsUsingBearingVector &,
                             Params>::value) {
      const Vec3d bearing_ = (rt.R * point + rt.t).normalized();

      const double thr = std::get_if<PnpInlierCheckParamsUsingBearingVector>(
                             &inlier_check_params_)
                             ->inlier_thr;

      return acos(std::clamp(bearing.dot(bearing_), -1.0, 1.0)) <
             thr;  // Epipolar constraint
      // Note: acos(x) <- [0, pi]
    } else {
      // https://qiita.com/saka1_p/items/e8c4dfdbfa88449190c5
      static_assert(false_v<Params>);
    }
  }

private:
  PnpInlierCheckParams inlier_check_params_;
};

}  // namespace gsrap

#endif  // GSRAP_INTERNAL_COMMON_H
