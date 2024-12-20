#include "gsrap/essential_solver.h"

#include <assert.h>

#include <algorithm>
#include <memory>
#include <random>
#include <unordered_set>
#include <utility>

#include "gsrap/macros.h"
#include "gsrap_internal/common.h"
#include "gsrap_internal/essential_solver/fivept_nister.h"
#include "gsrap_internal/ransac.h"
#include "gsrap_internal/sample.h"

GSRAP_IGNORE_STRICT_WARNING_PUSH

#include "Eigen/SVD"

// opengv
#include "opengv/relative_pose/modules/main.hpp"

GSRAP_IGNORE_STRICT_WARNING_POP

namespace gsrap {

// Model
struct ERt {
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  ERt() = default;
  ERt(Mat3d E_, Mat3d R_, Vec3d t_)
      : E(std::move(E_)), R(std::move(R_)), t(std::move(t_)) {}

  Mat3d E;  // Essential matrix
  Mat3d R;  // Rotation matrix
  Vec3d t;  // Translation vector
};

// This function determines if two bearing vectors intersect in the forward.
static bool CheckFront(Vec3d b1, const Vec3d& b2, const Mat3d& R21,
                       const Vec3d& t21) {
  b1 = R21 * b1;

  return t21.cross(b2).dot(t21.cross(b1)) > 0 && t21.dot(b2) > t21.dot(b1);
}

// This function checks EssentialSolverPolicy.
static EssentialSolverPolicy CheckEssentialSolverPolicy(
    EssentialSolverPolicy essential_solver_policy) {
  EssentialSolverPolicy ret = essential_solver_policy;

  // 'inlier_thr' should be greater than zero.
  ret.inlier_thr = std::max(0.0, ret.inlier_thr);
  // 'check_singular_value_thr' should be greater than zero.
  ret.check_singular_value_thr = std::max(0.0, ret.check_singular_value_thr);

  return ret;
}

/// Function object to compute essential matrix, rotation matrix, and
/// translation vector.
template <typename Iterator>
struct FivePt {
  std::vector<std::unique_ptr<ERt>> operator()(
      const eigen_vector<Eigen::Vector3d>& bearings1,
      const eigen_vector<Eigen::Vector3d>& bearings2,
      const std::vector<Iterator>& data, const bool check_singular_value = true,
      const double singular_value_thr = 1e15) const {
    assert(data.size() == 5);
    if (data.size() < 5) {
      return {};
    }

    /*
     * 1. Compute essential matrix from 5 pairs of bearing vectors.
     */
    std::vector<Mat3d, Eigen::aligned_allocator<Mat3d>> Es =
        FivePtNister(bearings1, bearings2, data);

    std::vector<std::unique_ptr<ERt>> ret;

    /*
     * 2. Essential matrices computed in the previous process
     */
    for (auto& E : Es) {
      /// 2.1 Perform singular value decomposition.
      Eigen::JacobiSVD<Mat3d> SVD_E(E,
                                    Eigen::ComputeFullU | Eigen::ComputeFullV);
      const double scale = SVD_E.singularValues()[0];
      if (check_singular_value) {
        /// Check if the three singular values are {|t|, |t|, 0}.
        if (abs(1.0 - SVD_E.singularValues()[1] / scale) > singular_value_thr) {
          continue;
        }
        if (abs(SVD_E.singularValues()[2] / scale) > singular_value_thr) {
          continue;
        }
      }
      // Normalize essential matrix.
      E /= scale;

      const Mat3d& U  = SVD_E.matrixU();
      const Mat3d& Vt = SVD_E.matrixV().transpose();

      const auto IsDetOne = [](const Mat3d& M) {
        assert((M.row(0).cross(M.row(1)) - M.row(2)).squaredNorm() < 1e-10 ||
               abs((M.row(0).cross(M.row(1)) - M.row(2)).squaredNorm() - 4.0) <
                   1e-10);
        return (M.row(0).cross(M.row(1)) - M.row(2)).squaredNorm() < 1.0;
      };

      // Adjust the determinant of the matrices computed later to be 1.
      Mat3d W;
      W << 0.0, -1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
          (IsDetOne(U) ^ /*xor*/ IsDetOne(Vt) ? -1.0 : 1.0);
      const Mat3d Wt = W.transpose();

      const Vec3d u = U.col(2);

      const std::vector<std::pair<Mat3d, Vec3d>
#ifdef GSRAP_LESS_EIGEN_3_4
                        ,
                        Eigen::aligned_allocator<std::pair<Mat3d, Vec3d>>
#endif  // GSRAP_LESS_EIGEN_3_4
                        >
          Rts = {{U * W * Vt, u},
                 {U * W * Vt, -u},
                 {U * Wt * Vt, u},
                 {U * Wt * Vt, -u}};

      // Determines if pairs of bearing vectors intersect in the forward.
      for (const auto& Rt : Rts) {
        const Mat3d& R = Rt.first;
        const Vec3d& t = Rt.second;

        if (std::all_of(data.begin(), data.end(),
                        [&R, &t, &bearings1, &bearings2](Iterator itr) {
                          return CheckFront(bearings1.at(itr->first),
                                            bearings2.at(itr->second), R, t);
                        })) {
          ret.emplace_back(std::make_unique<ERt>(E, R, t));
        }
      }
    }
    return ret;
  }
};

/// Function object to determine if a pair of bearing vectors is inlier.
template <typename Iterator>
struct IsInlier {
  bool operator()(const eigen_vector<Eigen::Vector3d>& bearings1,
                  const eigen_vector<Eigen::Vector3d>& bearings2,
                  const typename Iterator::value_type& match, const ERt& ert,
                  const double thr = 1e-15) const {
    assert(0 <= match.first && match.first < bearings1.size());
    assert(0 <= match.second && match.second < bearings2.size());

    const Vec3d& b1 = bearings1[match.first];
    const Vec3d& b2 = bearings2[match.second];

    return CheckFront(b1, b2, ert.R, ert.t) &&
           abs(b2.transpose() * ert.E * b1) < thr;  // Epipolar constraint
  }
};

template <typename Iterator>
std::pair<std::optional<EssentialSolverResult<Iterator>>, RansacReport>
ComputeEssentialMatrix(EssentialSolverPolicy essential_solver_policy,
                       const eigen_vector<Eigen::Vector3d>& bearings1,
                       const eigen_vector<Eigen::Vector3d>& bearings2,
                       Iterator first, Iterator last) {
  essential_solver_policy = CheckEssentialSolverPolicy(essential_solver_policy);

  FivePt<Iterator> five_pt;
  IsInlier<Iterator> is_inlier;
  SampleN<Iterator, 5> sample;

  const auto num_data = size_t(std::distance(first, last));

  // Execute RANSAC.
  auto result = Ransac(
      essential_solver_policy.ransac_policy, first, last,
      std::bind(five_pt, bearings1, bearings2, std::placeholders::_1,
                bool(essential_solver_policy.flags &
                     ESSENTIAL_SOLVER_POLICY_FLAGS_CHECK_SINGULAR_VALUE),
                essential_solver_policy.check_singular_value_thr),
      std::bind(is_inlier, bearings1, bearings2, std::placeholders::_1,
                std::placeholders::_2, essential_solver_policy.inlier_thr),
      std::bind(sample, num_data, std::placeholders::_1, std::placeholders::_2,
                std::placeholders::_3));

  // Pack the calculation results.
  std::optional<EssentialSolverResult<Iterator>> ret;
  if (result.first) {
    EssentialSolverResult<Iterator> tmp;
    tmp.essential   = result.first.value().model->E;
    tmp.rotation    = result.first.value().model->R;
    tmp.translation = result.first.value().model->t;

    tmp.inliers      = std::move(result.first->inliers);
    tmp.inlier_ratio = result.first->inlier_ratio;

    ret = tmp;
  }

  return std::make_pair(ret, result.second);
}

template std::pair<std::optional<EssentialSolverResult<std::vector<
                       std::pair<uint32_t, uint32_t>>::const_iterator>>,
                   RansacReport>
ComputeEssentialMatrix(
    EssentialSolverPolicy essential_solver_policy,
    const eigen_vector<Eigen::Vector3d>& bearings1,
    const eigen_vector<Eigen::Vector3d>& bearings2,
    std::vector<std::pair<uint32_t, uint32_t>>::const_iterator first,
    std::vector<std::pair<uint32_t, uint32_t>>::const_iterator last);

template std::pair<std::optional<EssentialSolverResult<
                       std::unordered_map<uint32_t, uint32_t>::const_iterator>>,
                   RansacReport>
ComputeEssentialMatrix(
    EssentialSolverPolicy essential_solver_policy,
    const eigen_vector<Eigen::Vector3d>& bearings1,
    const eigen_vector<Eigen::Vector3d>& bearings2,
    std::unordered_map<uint32_t, uint32_t>::const_iterator first,
    std::unordered_map<uint32_t, uint32_t>::const_iterator last);

}  // namespace gsrap
