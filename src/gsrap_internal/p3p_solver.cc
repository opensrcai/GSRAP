#include "gsrap/p3p_solver.h"

#include <cfloat>
//
#include <complex>  // FIXME(kyawakyawa)
#include <memory>
#include <type_traits>
#include <unordered_map>
#include <vector>

#include "gsrap/gsrap-def.h"
#include "gsrap/macros.h"
#include "gsrap_internal/common.h"
#include "gsrap_internal/p3p_solver/compute_poses.h"
#include "gsrap_internal/ransac.h"
#include "gsrap_internal/sample.h"

GSRAP_IGNORE_STRICT_WARNING_PUSH
#include "Eigen/Core"
GSRAP_IGNORE_STRICT_WARNING_POP

#define KAHAN_ADD(sum, add)    \
  diff   = add - remain;       \
  tmp    = sum + diff;         \
  remain = (tmp - sum) - diff; \
  sum    = tmp

#define KAHAN_ADD_LAST(sum, add) \
  diff = add - remain;           \
  tmp  = sum + diff;             \
  sum  = tmp

namespace gsrap {

static void refine_cubic_equation_solution(const double b, const double c,
                                           const double d, double *x) {
  const double x2 = *x * *x;

  double numerator   = 0;
  double denominator = 0;
  double diff        = 0.0;
  double tmp         = 0.0;
  double remain      = 0.0;

  remain = 0;
  KAHAN_ADD(numerator, d);
  KAHAN_ADD(numerator, c * *x);
  KAHAN_ADD(numerator, b * x2);
  KAHAN_ADD_LAST(numerator, *x * x2);

  remain = 0;
  KAHAN_ADD(denominator, c);
  KAHAN_ADD(denominator, 2.0 * b * *x);
  KAHAN_ADD_LAST(denominator, 3.0 * x2);

  if (fabs(denominator) >
      DBL_EPSILON) {  // TODO(kyawakyawa) Use the better way with Ofast option
    *x -= numerator / denominator;
  }
}

static void compute_real_solution_of_quadratic_equation(double b,
                                                        const double c,
                                                        double *solutions,
                                                        size_t *num_solutions,
                                                        const double D_thr) {
  // ref http://www.math.twcu.ac.jp/ogita/lec/na_basic.pdf
  // x^2  + b x + c =  0;
  //     |
  //     v
  // x^2 + 2b x + c =  0;
  //
  b *= 0.5;

  const double D = b * b - c;
  if (D > D_thr) {
    *num_solutions        = 2;
    const double sqD      = sqrt(D);
    const bool signbit_b  = signbit(b);  // if b >= 0 => 0 else b < 0 => 1
    const double sign_b   = signbit_b ? -1.0 : 1.0;
    solutions[signbit_b]  = (-b - sign_b * sqD);
    solutions[!signbit_b] = c / solutions[signbit_b];
  } else if (D > -D_thr) {
    *num_solutions = 1;
    solutions[0]   = -b;
  } else {
    *num_solutions = 0;
  }
}

static void compute_one_of_real_non_negative_solution_of_cubic_equation(
    double b, const double c, const double d, double D_thr, double *solution,
    bool *find_root) {
  // x^3 + b x^2  + c x + d =  0;
  //             |
  //             v
  // x^3 + 3b x^2  + c x + d =  0;

  b *= 1.0 / 3.0;

  const double b2 = b * b;

  const double p = b2 - c * (1.0 / 3.0);
  const double q = (b * c - 2.0 * b * b2 - d) * 0.5;

  // y^3 - 3 p y - 2 q = 0
  //       where y = x + b

  const double p3 = p * p * p;
  const double D  = q * q - p3;

  double x_candidate = -1.0;
  if (D < -D_thr) {
    assert(p > 0.0);
    // This solution is larger than the other two real-valued solutions.
    x_candidate = 2.0 * sqrt(p) * cos(acos(q * sqrt(1.0 / p3)) / 3.0) - b;

  } else {
    const double sqD = (D < D_thr ? 0.0 : sqrt(D));
    const double alpha_3 =
        (!signbit(q) ? cbrt(q + sqD)
                     : cbrt(q - sqD));  //  Make sure that the absolute value of
                                        //  alpha3 is large
    const double beta_3 =
        (D < D_thr
             ? alpha_3
             : p / alpha_3);  // When it is multiple solutions, it is close to
                              // zero division, so it is divided into cases.

    const double y1 = alpha_3 + beta_3;

    x_candidate = y1 - b;

    if (signbit(x_candidate) && D < D_thr) {
      x_candidate = y1 * -0.5 - b;
    }
  }
  if (signbit(x_candidate)) {
    *find_root = false;
  } else {
    *solution  = x_candidate;
    *find_root = true;
  }
}

// NOLINTNEXTLINE
void compute_real_solution_of_quartic_equation(
    double b, const double c, const double d, const double e,
    double solutions[4], size_t *num_solutions, const double D_thr,
    const double biquadratic_equation_thr, const double multiple_solution_thr);

void compute_real_solution_of_quartic_equation(
    double b, const double c, const double d, const double e,
    double solutions[4], size_t *num_solutions, const double D_thr,
    const double biquadratic_equation_thr, const double multiple_solution_thr) {
  b *= 0.25;

  const double b2 = b * b;

  const double p = c - 6.0 * b2;
  const double q = d - 2.0 * c * b + 8.0 * b2 * b;
  const double r = e - d * b + c * b2 - 3.0 * b2 * b2;

  if (fabs(q) < biquadratic_equation_thr) {
    double quadratic_equation_solutions[2];
    size_t num_quadratic_equation_solutions = 0;
    compute_real_solution_of_quadratic_equation(
        p, r, quadratic_equation_solutions, &num_quadratic_equation_solutions,
        D_thr);
    *num_solutions = 0;
    if (num_quadratic_equation_solutions >= 1) {
      bool zero_flag = false;  // To deal with floating point arithmetic errors.
      if (quadratic_equation_solutions[0] > multiple_solution_thr) {
        const double sq               = sqrt(quadratic_equation_solutions[0]);
        solutions[(*num_solutions)++] = -sq;
        solutions[(*num_solutions)++] = sq;
      } else if (quadratic_equation_solutions[0] > -multiple_solution_thr) {
        solutions[(*num_solutions)++] = 0.0;
        zero_flag                     = true;
      }

      if (num_quadratic_equation_solutions == 2) {
        if (quadratic_equation_solutions[1] > multiple_solution_thr) {
          const double sq               = sqrt(quadratic_equation_solutions[1]);
          solutions[(*num_solutions)++] = -sq;
          solutions[(*num_solutions)++] = sq;
        } else if (!zero_flag &&
                   quadratic_equation_solutions[1] > -multiple_solution_thr) {
          solutions[(*num_solutions)++] = 0.0;
        }
      }
    }
  } else {
    double u        = 0.0;
    bool find_root  = false;
    const double c1 = 2.0 * p;
    const double c2 = p * p - 4.0 * r;
    const double c3 = -q * q;

    compute_one_of_real_non_negative_solution_of_cubic_equation(
        c1, c2, c3, D_thr, &u, &find_root);

    if (!find_root) {
      *num_solutions = 0;
      return;
    }

    refine_cubic_equation_solution(c1, c2, c3, &u);

    const double sq_u      = sqrt(u);
    const double alpha     = (p + u) * 0.5;
    const double beta      = (q * 0.5) / u;
    const double sq_u_beta = sq_u * beta;

    size_t num_quadratic_equation_solutions0 = 0;

    compute_real_solution_of_quadratic_equation(
        sq_u, alpha - sq_u_beta, solutions, &num_quadratic_equation_solutions0,
        D_thr);

#if 0
    double quadratic_equation_solutions1[2];
    size_t num_quadratic_equation_solutions1 = 0;

    compute_real_solution_of_quadratic_equation(
        -sq_u, alpha + sq_u_beta, quadratic_equation_solutions1,
        &num_quadratic_equation_solutions1, D_thr);

    char duplicate_with_solution0[2] = {0, 0};
    for (size_t i = 0; i < num_quadratic_equation_solutions0; ++i) {
      for (size_t j = 0; j < num_quadratic_equation_solutions1; ++j) {
        if (fabs(solutions[i] - quadratic_equation_solutions1[j]) <
            multiple_solution_thr) {
          duplicate_with_solution0[j] = 1;
        }
      }
    }

    *num_solutions = num_quadratic_equation_solutions0;

    if (num_quadratic_equation_solutions1 >= 1 &&
        (duplicate_with_solution0[0] == 0)) {
      solutions[(*num_solutions)++] = quadratic_equation_solutions1[0];
    }
    if (num_quadratic_equation_solutions1 == 2 &&
        (duplicate_with_solution0[1] == 0)) {
      solutions[(*num_solutions)++] = quadratic_equation_solutions1[1];
    }
#else
    size_t num_quadratic_equation_solutions1 = 0;

    compute_real_solution_of_quadratic_equation(
        -sq_u, alpha + sq_u_beta, solutions + num_quadratic_equation_solutions0,
        &num_quadratic_equation_solutions1, D_thr);
    *num_solutions =
        num_quadratic_equation_solutions0 + num_quadratic_equation_solutions1;
#endif
  }

  for (size_t i = 0; i < *num_solutions; ++i) {
    // TODO(kyawakyawa) SIMD
    solutions[i] -= b;  // NOLINT
  }
}
template <typename Iterator>
struct P3pKe2017 {
  std::vector<std::unique_ptr<Rt>> operator()(

      const eigen_vector<Eigen::Vector3d> &bearings,
      const eigen_vector<Eigen::Vector3d> &points,
      const std::vector<Iterator> &data) const {
    assert(data.size() == 3);
    if (data.size() != 3) {
      // TODO(kyawakyawa) error handling
      return std::vector<std::unique_ptr<Rt>>();
    }

    double feature_vectors[3][4];  // TODO(kyawakyawa) Delete P4p mode
    double world_points[3][4];     // TODO(kyawakyawa) Delete P4p mode

    const auto df0        = data[0]->first;
    feature_vectors[0][0] = bearings[df0][0];
    feature_vectors[1][0] = bearings[df0][1];
    feature_vectors[2][0] = bearings[df0][2];

    const auto df1        = data[1]->first;
    feature_vectors[0][1] = bearings[df1][0];
    feature_vectors[1][1] = bearings[df1][1];
    feature_vectors[2][1] = bearings[df1][2];

    const auto df2        = data[2]->first;
    feature_vectors[0][2] = bearings[df2][0];
    feature_vectors[1][2] = bearings[df2][1];
    feature_vectors[2][2] = bearings[df2][2];

    feature_vectors[0][3] = 0.0;  // Not used
    feature_vectors[1][3] = 0.0;  // Not used
    feature_vectors[2][3] = 0.0;  // Not used

    const auto ds0     = data[0]->second;
    world_points[0][0] = points[ds0][0];
    world_points[1][0] = points[ds0][1];
    world_points[2][0] = points[ds0][2];

    const auto ds1     = data[1]->second;
    world_points[0][1] = points[ds1][0];
    world_points[1][1] = points[ds1][1];
    world_points[2][1] = points[ds1][2];

    const auto ds2     = data[2]->second;
    world_points[0][2] = points[ds2][0];
    world_points[1][2] = points[ds2][1];
    world_points[2][2] = points[ds2][2];

    world_points[0][3] = 0.0;  // Not used
    world_points[1][3] = 0.0;  // Not used
    world_points[2][3] = 0.0;  // Not used

    double R[4][3][3];
    double t[4][3];

    const int num_solutions = computePoses(feature_vectors, world_points, R, t);

    std::vector<std::unique_ptr<Rt>> ret;
    ret.reserve(size_t(num_solutions));
    for (size_t i = 0; i < size_t(num_solutions); ++i) {
      ret.emplace_back(std::make_unique<Rt>());

      Mat3d &ref_R = ret.back()->R;
      ref_R(0, 0)  = R[i][0][0];
      ref_R(0, 1)  = R[i][0][1];
      ref_R(0, 2)  = R[i][0][2];

      ref_R(1, 0) = R[i][1][0];
      ref_R(1, 1) = R[i][1][1];
      ref_R(1, 2) = R[i][1][2];

      ref_R(2, 0) = R[i][2][0];
      ref_R(2, 1) = R[i][2][1];
      ref_R(2, 2) = R[i][2][2];

      Vec3d &ref_t = ret.back()->t;
      ref_t[0]     = t[i][0];
      ref_t[1]     = t[i][1];
      ref_t[2]     = t[i][2];
    }

    return ret;
  }
};

// This function checks P3pSolverPolicy.
static P3pSolverPolicy CheckP3pSolverPolicy(P3pSolverPolicy p3p_solver_policy) {
  P3pSolverPolicy &ret = p3p_solver_policy;

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
std::pair<std::optional<P3pSolverResult<Iterator>>, RansacReport>
SolveP3pProblem(P3pSolverPolicy p3p_solver_policy,
                const eigen_vector<Vec3d> &bearings,
                const eigen_vector<Vec3d> &points, Iterator first,
                Iterator last) {
  p3p_solver_policy = CheckP3pSolverPolicy(p3p_solver_policy);

  P3pKe2017<Iterator> p3p_solver;
  SampleN<Iterator, 3> sample;
  const auto num_data = size_t(std::distance(first, last));

  auto result = std::visit(
      [&](const auto &params) {
        IsInlierForPnp<decltype(params), Iterator> is_inlier(params);

        return Ransac(
            p3p_solver_policy.ransac_policy, first, last,
            std::bind(p3p_solver, bearings, points, std::placeholders::_1),
            std::bind(is_inlier, bearings, points, std::placeholders::_1,
                      std::placeholders::_2),
            std::bind(sample, num_data, std::placeholders::_1,
                      std::placeholders::_2, std::placeholders::_3));
      },
      p3p_solver_policy.pnp_inlier_check_params);

  // Refinement is not required or failed to that, return the best model of
  // RANSAC Return the best model when it available
  if (result.first) {
    P3pSolverResult<Iterator> ret;
    ret.rotation    = result.first.value().model->R;
    ret.translation = result.first.value().model->t;

    ret.inliers      = std::move(result.first->inliers);
    ret.inlier_ratio = result.first->inlier_ratio;

    return std::make_pair(ret, result.second);
  }

  return std::make_pair(std::nullopt, RansacReport());
}

template std::pair<std::optional<P3pSolverResult<std::vector<
                       std::pair<uint32_t, uint32_t>>::const_iterator>>,
                   RansacReport>
SolveP3pProblem(
    P3pSolverPolicy p3p_solver_policy, const eigen_vector<Vec3d> &bearings,
    const eigen_vector<Vec3d> &points,
    std::vector<std::pair<uint32_t, uint32_t>>::const_iterator first,
    std::vector<std::pair<uint32_t, uint32_t>>::const_iterator last);

template std::pair<std::optional<P3pSolverResult<
                       std::unordered_map<uint32_t, uint32_t>::const_iterator>>,
                   RansacReport>
SolveP3pProblem(P3pSolverPolicy p3p_solver_policy,
                const eigen_vector<Vec3d> &bearings,
                const eigen_vector<Vec3d> &points,
                std::unordered_map<uint32_t, uint32_t>::const_iterator first,
                std::unordered_map<uint32_t, uint32_t>::const_iterator last);

eigen_vector<std::pair<Mat3d, Vec3d>> SolveP3pProblemWithoutRansac(
    const eigen_vector<Vec3d> &bearings, const eigen_vector<Vec3d> &points) {
  const size_t num = std::min(bearings.size(), points.size());

  assert(bearings.size() == points.size());
  assert(num == 3);
  if (num != 3) {
    return {};
  }

  double feature_vectors[3][4];  // TODO(kyawakyawa) Delete P4p mode
  double world_points[3][4];     // TODO(kyawakyawa) Delete P4p mode

  feature_vectors[0][0] = bearings[0][0];
  feature_vectors[1][0] = bearings[0][1];
  feature_vectors[2][0] = bearings[0][2];

  feature_vectors[0][1] = bearings[1][0];
  feature_vectors[1][1] = bearings[1][1];
  feature_vectors[2][1] = bearings[1][2];

  feature_vectors[0][2] = bearings[2][0];
  feature_vectors[1][2] = bearings[2][1];
  feature_vectors[2][2] = bearings[2][2];

  feature_vectors[0][3] = 0.0;  // Not used
  feature_vectors[1][3] = 0.0;  // Not used
  feature_vectors[2][3] = 0.0;  // Not used

  world_points[0][0] = points[0][0];
  world_points[1][0] = points[0][1];
  world_points[2][0] = points[0][2];

  world_points[0][1] = points[1][0];
  world_points[1][1] = points[1][1];
  world_points[2][1] = points[1][2];

  world_points[0][2] = points[2][0];
  world_points[1][2] = points[2][1];
  world_points[2][2] = points[2][2];

  world_points[0][3] = 0.0;  // Not used
  world_points[1][3] = 0.0;  // Not used
  world_points[2][3] = 0.0;  // Not used

  double R[4][3][3];
  double t[4][3];

  const int num_solutions = computePoses(feature_vectors, world_points, R, t);

  eigen_vector<std::pair<Mat3d, Vec3d>> ret;
  ret.reserve(size_t(num_solutions));
  for (size_t i = 0; i < size_t(num_solutions); ++i) {
    ret.emplace_back();

    Mat3d &ref_R = ret.back().first;
    ref_R(0, 0)  = R[i][0][0];
    ref_R(0, 1)  = R[i][0][1];
    ref_R(0, 2)  = R[i][0][2];

    ref_R(1, 0) = R[i][1][0];
    ref_R(1, 1) = R[i][1][1];
    ref_R(1, 2) = R[i][1][2];

    ref_R(2, 0) = R[i][2][0];
    ref_R(2, 1) = R[i][2][1];
    ref_R(2, 2) = R[i][2][2];

    Vec3d &ref_t = ret.back().second;
    ref_t[0]     = t[i][0];
    ref_t[1]     = t[i][1];
    ref_t[2]     = t[i][2];
  }
  return ret;
}

}  // namespace gsrap
