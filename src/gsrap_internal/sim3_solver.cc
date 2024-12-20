#include "gsrap/sim3_solver.h"

#include <math.h>
#include <stdint.h>

#include <limits>
#include <unordered_set>
#include <utility>

#include "gsrap/macros.h"
#include "gsrap_internal/common.h"
#include "gsrap_internal/ransac.h"
#include "gsrap_internal/sample.h"

GSRAP_IGNORE_STRICT_WARNING_PUSH

#include "Eigen/Dense"
#include "Eigen/SVD"

GSRAP_IGNORE_STRICT_WARNING_POP

namespace gsrap {

// Model
struct Sim3 {
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  inline Sim3() = default;
  inline Sim3(const double& s_, Mat3d R_, Vec3d t_)
      : R(std::move(R_)), t(std::move(t_)), s(s_) {}

  Mat3d R;         // Rotation matrix
  Vec3d t;         // Translation vector
  double s = 1.0;  // Scale

  double simga_y = 0.0;
};

/// Function object to determine if a pair of points is inlier.
template <typename Iterator>
struct IsInlier {
  bool operator()(const eigen_vector<Eigen::Vector3d>& points1,
                  const eigen_vector<Eigen::Vector3d>& points2,
                  const typename Iterator::value_type& match, const Sim3& srt,
                  const double thr = 1e-12) const {
    assert(0 <= match.first && match.first < points1.size());
    assert(0 <= match.second && match.second < points2.size());

    const Vec3d& p1 = points1[match.first];
    const Vec3d& p2 = points2[match.second];

    const double& s = srt.s;
    const Mat3d& R  = srt.R;
    const Vec3d& t  = srt.t;

    return (s * R * p1 + t - p2).norm() < thr * srt.simga_y;
  }
};

/// Function object to compute Sim3 transformation.
template <typename Iterator>
struct Umeyama {
  std::vector<std::unique_ptr<Sim3>> operator()(
      const eigen_vector<Eigen::Vector3d>& points1,
      const eigen_vector<Eigen::Vector3d>& points2,
      const std::vector<Iterator>& data) const {
    const int n                    = int(data.size());
    const double inv_n             = 1.0 / double(n);
    Eigen::Matrix<double, 3, -1> x = Eigen::Matrix<double, 3, -1>::Zero(3, n);
    Eigen::Matrix<double, 3, -1> y = Eigen::Matrix<double, 3, -1>::Zero(3, n);

    for (int i = 0; i < n; ++i) {
      x.col(i) = points1[data[size_t(i)]->first];
    }
    for (int i = 0; i < n; ++i) {
      y.col(i) = points2[data[size_t(i)]->second];
    }

    const Vec3d mu_x = inv_n * x.rowwise().sum();
    const Vec3d mu_y = inv_n * y.rowwise().sum();

    const Eigen::Matrix<double, 3, -1> x_m_mu_x = x.colwise() - mu_x;
    const Eigen::Matrix<double, 3, -1> y_m_mu_y = y.colwise() - mu_y;

    Mat3d sigma_xy = Mat3d::Zero();
    for (int i = 0; i < n; ++i) {
      sigma_xy += y_m_mu_y.col(i) * x_m_mu_x.col(i).transpose();
    }
    sigma_xy *= inv_n;

    Eigen::JacobiSVD<Mat3d> svd(sigma_xy,
                                Eigen::ComputeFullU | Eigen::ComputeFullV);
    if (svd.rank() < 2) {
      return std::vector<std::unique_ptr<Sim3>>();
    }
    const Mat3d S = [&svd](void) {
      const auto IsDetOne = [](const Mat3d& M) {
        return abs(M.determinant() - 1.0) < 1e-8;
      };

      Mat3d ret = Mat3d::Identity();
      if (IsDetOne(svd.matrixU()) ^ /*xor*/ IsDetOne(svd.matrixV())) {
        ret(2, 2) = -1;
      }
      return ret;
    }();

    std::unique_ptr<Sim3> srt(new Sim3);
    srt->R = svd.matrixU() * S * svd.matrixV().transpose();
    const double inv_sigma2_x =
        1.0 / (inv_n * x_m_mu_x.colwise().squaredNorm().sum());
    const double trDS =
        (Eigen::DiagonalMatrix<double, 3>(svd.singularValues()) * S).trace();
    srt->s = inv_sigma2_x * trDS;
    srt->t = mu_y - srt->s * srt->R * mu_x;

    const double sigma2_y = inv_n * y_m_mu_y.colwise().squaredNorm().sum();

    srt->simga_y = sqrt(sigma2_y);

    std::vector<std::unique_ptr<Sim3>> ret;
    ret.emplace_back(std::move(srt));

    return ret;
  }
};

template <typename Iterator>
std::pair<std::optional<Sim3SolverResult<Iterator>>, RansacReport>
ComputeSim3Transformation(Sim3SolverPolicy sim3_solver_policy,
                          const eigen_vector<Vec3d>& points1,
                          const eigen_vector<Vec3d>& points2, Iterator first,
                          Iterator last) {
  Umeyama<Iterator> umeyama;
  IsInlier<Iterator> is_inlier;
  Sample<Iterator> sample;

  const auto num_data = size_t(std::distance(first, last));

  // Execute RANSAC.
  auto result =
      Ransac(sim3_solver_policy.ransac_policy, first, last,
             std::bind(umeyama, points1, points2, std::placeholders::_1),
             std::bind(is_inlier, points1, points2, std::placeholders::_1,
                       std::placeholders::_2, sim3_solver_policy.inlier_thr),
             std::bind(sample, num_data, sim3_solver_policy.num_sample,
                       std::placeholders::_1, std::placeholders::_2,
                       std::placeholders::_3));

  // Refine model by using all inliers after RANSAC
  const bool do_refine =
      result.first && result.first->num_inliers >= 3 &&
      (sim3_solver_policy.flags & SIM3_SOLVER_POLICY_FLAGS_REFINE);
  if (do_refine) {
    const std::vector<std::unique_ptr<Sim3>> resolved =
        umeyama(points1, points2, result.first->inliers);

    if (!resolved.empty()) {
      Sim3SolverResult<Iterator> ret;

      ret.scale       = resolved[0]->s;
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
    Sim3SolverResult<Iterator> ret;

    ret.scale       = result.first.value().model->s;
    ret.rotation    = result.first.value().model->R;
    ret.translation = result.first.value().model->t;

    ret.inliers      = std::move(result.first->inliers);
    ret.inlier_ratio = result.first->inlier_ratio;

    return std::make_pair(ret, result.second);
  }

  // Or return empty resulim
  return std::make_pair(std::nullopt, result.second);
}

template std::pair<std::optional<Sim3SolverResult<std::vector<
                       std::pair<uint32_t, uint32_t>>::const_iterator>>,
                   RansacReport>
ComputeSim3Transformation(
    Sim3SolverPolicy sim3_solver_policy, const eigen_vector<Vec3d>& point1,
    const eigen_vector<Vec3d>& point2,
    std::vector<std::pair<uint32_t, uint32_t>>::const_iterator first,
    std::vector<std::pair<uint32_t, uint32_t>>::const_iterator last);

template std::pair<std::optional<Sim3SolverResult<
                       std::unordered_map<uint32_t, uint32_t>::const_iterator>>,
                   RansacReport>
ComputeSim3Transformation(
    Sim3SolverPolicy sim3_solver_policy, const eigen_vector<Vec3d>& points1,
    const eigen_vector<Vec3d>& points2,
    std::unordered_map<uint32_t, uint32_t>::const_iterator first,
    std::unordered_map<uint32_t, uint32_t>::const_iterator last);

}  // namespace gsrap
