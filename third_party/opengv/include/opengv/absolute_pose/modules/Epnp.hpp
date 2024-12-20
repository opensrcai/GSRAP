/******************************************************************************
 * Author:   Laurent Kneip                                                    *
 * Contact:  kneip.laurent@gmail.com                                          *
 * License:  Copyright (c) 2013 Laurent Kneip, ANU. All rights reserved.      *
 *                                                                            *
 * Redistribution and use in source and binary forms, with or without         *
 * modification, are permitted provided that the following conditions         *
 * are met:                                                                   *
 * * Redistributions of source code must retain the above copyright           *
 *   notice, this list of conditions and the following disclaimer.            *
 * * Redistributions in binary form must reproduce the above copyright        *
 *   notice, this list of conditions and the following disclaimer in the      *
 *   documentation and/or other materials provided with the distribution.     *
 * * Neither the name of ANU nor the names of its contributors may be         *
 *   used to endorse or promote products derived from this software without   *
 *   specific prior written permission.                                       *
 *                                                                            *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"*
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE  *
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE *
 * ARE DISCLAIMED. IN NO EVENT SHALL ANU OR THE CONTRIBUTORS BE LIABLE        *
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL *
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR *
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER *
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT         *
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY  *
 * OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF     *
 * SUCH DAMAGE.                                                               *
 ******************************************************************************/

// Note: this code has been downloaded from the homepage of the "Computer
// Vision Laboratory" at EPFL Lausanne, and was originally developped by the
// authors of [4]. I only adapted it to Eigen.

#ifndef OPENGV_ABSOLUTE_POSE_MODULES_EPNP_HPP_
#define OPENGV_ABSOLUTE_POSE_MODULES_EPNP_HPP_

#include <Eigen/src/Core/util/DisableStupidWarnings.h>
#include <stdlib.h>

#include <Eigen/Eigen>

namespace opengv {
namespace absolute_pose {
namespace modules {

class Epnp {
  using Mat3d = Eigen::Matrix3d;
  using Vec3d = Eigen::Vector3d;
  using Vec4d = Eigen::Vector4d;

public:
  Epnp(void);
  ~Epnp();

  void set_maximum_number_of_correspondences(const int n);
  void reset_correspondences();
  void add_correspondence(const Vec3d& point, const Vec3d& bearing);

  double compute_pose(Mat3d& R, Vec3d& t);

  double reprojection_error(const Mat3d& R, const Vec3d& t);

private:
  void choose_control_points();
  void compute_barycentric_coordinates();
  void fill_M(Eigen::MatrixXd& M, const int row, const Vec4d& alpha,
              const double u, const double v);
  void compute_ccs(const double* betas, const Eigen::MatrixXd& ut);
  void compute_pcs();

  void solve_for_sign();

  void find_betas_approx_1(const Eigen::Matrix<double, 6, 10>& L_6x10,
                           const Eigen::Matrix<double, 6, 1>& Rho,
                           double* betas);
  void find_betas_approx_2(const Eigen::Matrix<double, 6, 10>& L_6x10,
                           const Eigen::Matrix<double, 6, 1>& Rho,
                           double* betas);
  void find_betas_approx_3(const Eigen::Matrix<double, 6, 10>& L_6x10,
                           const Eigen::Matrix<double, 6, 1>& Rho,
                           double* betas);
  void qr_solve(Eigen::Matrix<double, 6, 4>& A, Eigen::Matrix<double, 6, 1>& b,
                Eigen::Matrix<double, 4, 1>& X);

  void compute_rho(Eigen::Matrix<double, 6, 1>& Rho);
  void compute_L_6x10(const Eigen::MatrixXd& Ut,
                      Eigen::Matrix<double, 6, 10>& L_6x10);

  void gauss_newton(const Eigen::Matrix<double, 6, 10>& L_6x10,
                    const Eigen::Matrix<double, 6, 1>& Rho,
                    double current_betas[4]);
  void compute_A_and_b_gauss_newton(const Eigen::Matrix<double, 6, 10>& L_6x10,
                                    const Eigen::Matrix<double, 6, 1>& Rho,
                                    const double cb[4],
                                    Eigen::Matrix<double, 6, 4>& A,
                                    Eigen::Matrix<double, 6, 1>& b);

  double compute_R_and_t(const Eigen::MatrixXd& Ut, const double* betas,
                         Mat3d& R, Vec3d& t);

  void estimate_R_and_t(Mat3d& R, Vec3d& t);

  static constexpr double uc = 0.0;
  static constexpr double vc = 0.0;
  static constexpr double fu = 1.0;
  static constexpr double fv = 1.0;

  std::vector<double> us;
  std::vector<double> vs;
  std::vector<Vec3d, Eigen::aligned_allocator<Vec3d>> pws;
  std::vector<Vec3d, Eigen::aligned_allocator<Vec3d>> pcs;
  std::vector<Vec4d, Eigen::aligned_allocator<Vec4d>> alphas;
  std::vector<int> signs;  // added!
  int maximum_number_of_correspondences;
  int number_of_correspondences;

  std::vector<Vec3d, Eigen::aligned_allocator<Vec3d>> cws;
  std::vector<Vec3d, Eigen::aligned_allocator<Vec3d>> ccs;
};

}  // namespace modules
}  // namespace absolute_pose
}  // namespace opengv

#endif /* OPENGV_ABSOLUTE_POSE_MODULES_EPNP_HPP_ */
