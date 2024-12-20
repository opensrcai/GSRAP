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

#include <iostream>
using namespace std;

#include <opengv/absolute_pose/modules/Epnp.hpp>

opengv::absolute_pose::modules::Epnp::Epnp() {
  maximum_number_of_correspondences = 0;
  number_of_correspondences         = 0;
}

opengv::absolute_pose::modules::Epnp::~Epnp() = default;

void opengv::absolute_pose::modules::Epnp::
    set_maximum_number_of_correspondences(int n) {
  if (maximum_number_of_correspondences < n) {
    maximum_number_of_correspondences = n;
    us.resize(maximum_number_of_correspondences);
    vs.resize(maximum_number_of_correspondences);
    pws.resize(maximum_number_of_correspondences);
    pcs.resize(maximum_number_of_correspondences);
    alphas.resize(maximum_number_of_correspondences);
    signs.resize(maximum_number_of_correspondences);
  }
}

void opengv::absolute_pose::modules::Epnp::reset_correspondences() {
  number_of_correspondences = 0;
}

void opengv::absolute_pose::modules::Epnp::add_correspondence(
    const Vec3d& point, const Vec3d& bearing) {
  pws.at(number_of_correspondences) = point;

  us.at(number_of_correspondences) = bearing(0) / bearing(2);
  vs.at(number_of_correspondences) = bearing(1) / bearing(2);

  // added the following
  if (bearing(2) > 0.0)
    signs[number_of_correspondences] = 1;
  else
    signs[number_of_correspondences] = -1;

  number_of_correspondences++;
}

void opengv::absolute_pose::modules::Epnp::choose_control_points() {
  // Take C0 as the reference points centroid:
  cws.resize(4);
  cws.at(0) = Vec3d::Zero();
  for (int i = 0; i < number_of_correspondences; i++)
    for (int j = 0; j < 3; j++) cws[0](j) += pws[i][j];  // TODO

  cws[0] /= number_of_correspondences;

  // Take C1, C2, and C3 from PCA on the reference points:
  Eigen::MatrixXd PW0(number_of_correspondences, 3);

  for (int i = 0; i < number_of_correspondences; i++)
    PW0.row(i) = pws.at(i) - cws.at(0);

  Eigen::MatrixXd PW0tPW0 = PW0.transpose() * PW0;
  Eigen::JacobiSVD<Eigen::MatrixXd> SVD(
      PW0tPW0, Eigen::ComputeFullV | Eigen::ComputeFullU);
  Eigen::MatrixXd D  = SVD.singularValues();
  Eigen::MatrixXd Ut = SVD.matrixU();

  for (int i = 0; i < 3; i++) {
    double k      = sqrt(D(i, 0) / number_of_correspondences);
    cws.at(i + 1) = cws.at(0) + k * Ut.col(i);
  }
}

void opengv::absolute_pose::modules::Epnp::compute_barycentric_coordinates() {
  Eigen::Matrix3d CC;

  for (int j = 0; j < 3; j++) CC.col(j) = cws.at(j + 1) - cws.at(0);

  Eigen::Matrix3d CC_inv = CC.fullPivLu().inverse();

  for (int i = 0; i < number_of_correspondences; i++) {
    const Vec3d& pi = pws.at(i);
    Vec4d& a        = alphas.at(i);

    a.bottomRows(3) = CC_inv * (pi - cws.at(0));
    a(0)            = 1.0f - a(1) - a(2) - a(3);
  }
}

void opengv::absolute_pose::modules::Epnp::fill_M(Eigen::MatrixXd& M,
                                                  const int row,
                                                  const Vec4d& alpha,
                                                  const double u,
                                                  const double v) {
  for (int i = 0; i < 4; i++) {
    M(row, 3 * i)     = alpha(i) * fu;
    M(row, 3 * i + 1) = 0.0;
    M(row, 3 * i + 2) = alpha(i) * (uc - u);

    M(row + 1, 3 * i)     = 0.0;
    M(row + 1, 3 * i + 1) = alpha(i) * fv;
    M(row + 1, 3 * i + 2) = alpha(i) * (vc - v);
  }
}

void opengv::absolute_pose::modules::Epnp::compute_ccs(
    const double* betas, const Eigen::MatrixXd& ut) {
  ccs.resize(4);
  for (int i = 0; i < 4; i++) ccs.at(i) = Vec3d::Zero();

  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++)
      ccs.at(j) += betas[i] * ut.block(11 - i, 3 * j, 1, 3).transpose();
  }
}

void opengv::absolute_pose::modules::Epnp::compute_pcs() {
  for (int i = 0; i < number_of_correspondences; i++) {
    const Vec4d& a = alphas.at(i);
    pcs.at(i)      = a(0) * ccs.at(0) + a(1) * ccs.at(1) + a(2) * ccs.at(2) +
                a(3) * ccs.at(3);
  }
}

double opengv::absolute_pose::modules::Epnp::compute_pose(Mat3d& R, Vec3d& t) {
  choose_control_points();
  compute_barycentric_coordinates();

  Eigen::MatrixXd M(2 * number_of_correspondences, 12);

  for (int i = 0; i < number_of_correspondences; i++)
    fill_M(M, 2 * i, alphas.at(i), us.at(i), vs.at(i));

  Eigen::MatrixXd MtM = M.transpose() * M;
  Eigen::JacobiSVD<Eigen::MatrixXd> SVD(MtM, Eigen::ComputeThinU);
  Eigen::MatrixXd Ut = SVD.matrixU().transpose();

  // Fill elements in null-space column via deterministic way
  // (SVD with ComputeFullU computes different results every time even if the
  // input is same)
  if (2 * number_of_correspondences < 12) {
    for (int L = 2 * number_of_correspondences; L < 12; L++) {
      Eigen::MatrixXd A      = Ut.topLeftCorner(L, L);
      Eigen::MatrixXd b      = Ut.block(0, L, L, 1);
      Eigen::MatrixXd x_part = -A.fullPivLu().inverse() * b;
      Eigen::MatrixXd x      = Eigen::MatrixXd::Zero(1, 12);
      x.leftCols(L)          = x_part.transpose();
      x(0, L)                = 1.0;
      Ut.block(L, 0, 1, 12)  = x.normalized();
    }
  }

  Eigen::Matrix<double, 6, 10> L_6x10;
  Eigen::Matrix<double, 6, 1> Rho;

  compute_L_6x10(Ut, L_6x10);
  compute_rho(Rho);

  double rep_errors[4];
  Mat3d Rs[4];
  Vec3d ts[4];
  double beta[4];

  find_betas_approx_1(L_6x10, Rho, beta);
  gauss_newton(L_6x10, Rho, beta);
  rep_errors[1] = compute_R_and_t(Ut, beta, Rs[1], ts[1]);

  find_betas_approx_2(L_6x10, Rho, beta);
  gauss_newton(L_6x10, Rho, beta);
  rep_errors[2] = compute_R_and_t(Ut, beta, Rs[2], ts[2]);

  find_betas_approx_3(L_6x10, Rho, beta);
  gauss_newton(L_6x10, Rho, beta);
  rep_errors[3] = compute_R_and_t(Ut, beta, Rs[3], ts[3]);

  int N = 1;
  if (rep_errors[2] < rep_errors[1]) N = 2;
  if (rep_errors[3] < rep_errors[N]) N = 3;

  R = Rs[N];
  t = ts[N];

  return rep_errors[N];
}

double opengv::absolute_pose::modules::Epnp::reprojection_error(
    const Mat3d& R, const Vec3d& t) {
  double sum2 = 0.0;

  for (int i = 0; i < number_of_correspondences; i++) {
    const Vec3d& pw = pws.at(i);
    double Xc       = R.row(0).dot(pw) + t[0];  // TODO
    double Yc       = R.row(1).dot(pw) + t[1];
    double inv_Zc   = 1.0 / (R.row(2).dot(pw) + t[2]);
    double ue       = uc + fu * Xc * inv_Zc;
    double ve       = vc + fv * Yc * inv_Zc;
    double u        = us.at(i);
    double v        = vs.at(i);

    sum2 += sqrt((u - ue) * (u - ue) + (v - ve) * (v - ve));
  }

  return sum2 / number_of_correspondences;
}

void opengv::absolute_pose::modules::Epnp::estimate_R_and_t(Mat3d& R,
                                                            Vec3d& t) {
  Vec3d pc0 = Vec3d::Zero();
  Vec3d pw0 = Vec3d::Zero();

  for (int i = 0; i < number_of_correspondences; i++) {
    const Vec3d& pc = pcs.at(i);
    const Vec3d& pw = pws.at(i);
    pc0 += pc;
    pw0 += pw;
  }
  pc0 /= number_of_correspondences;
  pw0 /= number_of_correspondences;

  Eigen::MatrixXd Abt(3, 3);

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) Abt(i, j) = 0.0;
  }

  for (int i = 0; i < number_of_correspondences; i++) {
    const Vec3d& pc = pcs.at(i);
    const Vec3d& pw = pws.at(i);
    Abt += (pc - pc0) * (pw - pw0).transpose();
  }

  Eigen::JacobiSVD<Eigen::MatrixXd> SVD(
      Abt, Eigen::ComputeFullV | Eigen::ComputeFullU);
  Eigen::MatrixXd Abt_u = SVD.matrixU();
  Eigen::MatrixXd Abt_v = SVD.matrixV();

  R = Abt_u * Abt_v.transpose();

  const double det = R.determinant();

  // change 1: negative determinant problem is solved by changing Abt_v, not R

  if (det < 0) {
    Eigen::MatrixXd Abt_v_prime = Abt_v;
    Abt_v_prime.col(2)          = -Abt_v.col(2);
    R                           = Abt_u * Abt_v_prime.transpose();
  }

  t = pc0 - R * pw0;
}

void opengv::absolute_pose::modules::Epnp::solve_for_sign() {
  // change to this (using original depths)
  if ((pcs.at(0)(2) < 0.0 && signs[0] > 0) ||
      (pcs.at(0)(2) > 0.0 && signs[0] < 0)) {
    for (int i = 0; i < 4; i++)
      for (int j = 0; j < 3; j++) ccs[i](j) = -ccs[i](j);

    for (int i = 0; i < number_of_correspondences; i++) {
      pcs.at(i) *= -1;
    }
  }
}

double opengv::absolute_pose::modules::Epnp::compute_R_and_t(
    const Eigen::MatrixXd& Ut, const double* betas, Mat3d& R, Vec3d& t) {
  compute_ccs(betas, Ut);
  compute_pcs();

  solve_for_sign();

  estimate_R_and_t(R, t);

  return reprojection_error(R, t);
}

// betas10        = [B11 B12 B22 B13 B23 B33 B14 B24 B34 B44]
// betas_approx_1 = [B11 B12     B13         B14]

void opengv::absolute_pose::modules::Epnp::find_betas_approx_1(
    const Eigen::Matrix<double, 6, 10>& L_6x10,
    const Eigen::Matrix<double, 6, 1>& Rho, double* betas) {
  Eigen::MatrixXd L_6x4(6, 4);

  for (int i = 0; i < 6; i++) {
    L_6x4(i, 0) = L_6x10(i, 0);
    L_6x4(i, 1) = L_6x10(i, 1);
    L_6x4(i, 2) = L_6x10(i, 3);
    L_6x4(i, 3) = L_6x10(i, 6);
  }

  Eigen::JacobiSVD<Eigen::MatrixXd> SVD(
      L_6x4, Eigen::ComputeFullV | Eigen::ComputeFullU);
  Eigen::VectorXd Rho_temp = Rho;
  Eigen::VectorXd b4       = SVD.solve(Rho_temp);

  if (b4[0] < 0) {
    betas[0] = sqrt(-b4[0]);
    betas[1] = -b4[1] / betas[0];
    betas[2] = -b4[2] / betas[0];
    betas[3] = -b4[3] / betas[0];
  } else {
    betas[0] = sqrt(b4[0]);
    betas[1] = b4[1] / betas[0];
    betas[2] = b4[2] / betas[0];
    betas[3] = b4[3] / betas[0];
  }
}

// betas10        = [B11 B12 B22 B13 B23 B33 B14 B24 B34 B44]
// betas_approx_2 = [B11 B12 B22                            ]

void opengv::absolute_pose::modules::Epnp::find_betas_approx_2(
    const Eigen::Matrix<double, 6, 10>& L_6x10,
    const Eigen::Matrix<double, 6, 1>& Rho, double* betas) {
  Eigen::MatrixXd L_6x3(6, 3);

  for (int i = 0; i < 6; i++) {
    L_6x3(i, 0) = L_6x10(i, 0);
    L_6x3(i, 1) = L_6x10(i, 1);
    L_6x3(i, 2) = L_6x10(i, 2);
  }

  Eigen::JacobiSVD<Eigen::MatrixXd> SVD(
      L_6x3, Eigen::ComputeFullV | Eigen::ComputeFullU);
  Eigen::VectorXd Rho_temp = Rho;
  Eigen::VectorXd b3       = SVD.solve(Rho_temp);

  if (b3[0] < 0) {
    betas[0] = sqrt(-b3[0]);
    betas[1] = (b3[2] < 0) ? sqrt(-b3[2]) : 0.0;
  } else {
    betas[0] = sqrt(b3[0]);
    betas[1] = (b3[2] > 0) ? sqrt(b3[2]) : 0.0;
  }

  if (b3[1] < 0) betas[0] = -betas[0];

  betas[2] = 0.0;
  betas[3] = 0.0;
}

// betas10        = [B11 B12 B22 B13 B23 B33 B14 B24 B34 B44]
// betas_approx_3 = [B11 B12 B22 B13 B23                    ]

void opengv::absolute_pose::modules::Epnp::find_betas_approx_3(
    const Eigen::Matrix<double, 6, 10>& L_6x10,
    const Eigen::Matrix<double, 6, 1>& Rho, double* betas) {
  Eigen::MatrixXd L_6x5(6, 5);

  for (int i = 0; i < 6; i++) {
    L_6x5(i, 0) = L_6x10(i, 0);
    L_6x5(i, 1) = L_6x10(i, 1);
    L_6x5(i, 2) = L_6x10(i, 2);
    L_6x5(i, 3) = L_6x10(i, 3);
    L_6x5(i, 4) = L_6x10(i, 4);
  }

  Eigen::JacobiSVD<Eigen::MatrixXd> SVD(
      L_6x5, Eigen::ComputeFullV | Eigen::ComputeFullU);
  Eigen::VectorXd Rho_temp = Rho;
  Eigen::VectorXd b5       = SVD.solve(Rho_temp);

  if (b5[0] < 0) {
    betas[0] = sqrt(-b5[0]);
    betas[1] = (b5[2] < 0) ? sqrt(-b5[2]) : 0.0;
  } else {
    betas[0] = sqrt(b5[0]);
    betas[1] = (b5[2] > 0) ? sqrt(b5[2]) : 0.0;
  }
  if (b5[1] < 0) betas[0] = -betas[0];
  betas[2] = b5[3] / betas[0];
  betas[3] = 0.0;
}

void opengv::absolute_pose::modules::Epnp::compute_L_6x10(
    const Eigen::MatrixXd& Ut, Eigen::Matrix<double, 6, 10>& L_6x10) {
  Vec3d dv[4][6];

  for (int i = 0; i < 4; i++) {
    int a = 0, b = 1;
    for (int j = 0; j < 6; j++) {
      dv[i][j] = (Ut.block(11 - i, 3 * a, 1, 3) - Ut.block(11 - i, 3 * b, 1, 3))
                     .transpose();

      b++;
      if (b > 3) {
        a++;
        b = a + 1;
      }
    }
  }

  for (int i = 0; i < 6; i++) {
    L_6x10(i, 0) = dv[0][i].dot(dv[0][i]);
    L_6x10(i, 1) = 2.0f * dv[0][i].dot(dv[1][i]);
    L_6x10(i, 2) = dv[1][i].dot(dv[1][i]);
    L_6x10(i, 3) = 2.0f * dv[0][i].dot(dv[2][i]);
    L_6x10(i, 4) = 2.0f * dv[1][i].dot(dv[2][i]);
    L_6x10(i, 5) = dv[2][i].dot(dv[2][i]);
    L_6x10(i, 6) = 2.0f * dv[0][i].dot(dv[3][i]);
    L_6x10(i, 7) = 2.0f * dv[1][i].dot(dv[3][i]);
    L_6x10(i, 8) = 2.0f * dv[2][i].dot(dv[3][i]);
    L_6x10(i, 9) = dv[3][i].dot(dv[3][i]);
  }
}

void opengv::absolute_pose::modules::Epnp::compute_rho(
    Eigen::Matrix<double, 6, 1>& Rho) {
  Rho[0] = (cws.at(0) - cws.at(1)).squaredNorm();
  Rho[1] = (cws.at(0) - cws.at(2)).squaredNorm();
  Rho[2] = (cws.at(0) - cws.at(3)).squaredNorm();
  Rho[3] = (cws.at(1) - cws.at(2)).squaredNorm();
  Rho[4] = (cws.at(1) - cws.at(3)).squaredNorm();
  Rho[5] = (cws.at(2) - cws.at(3)).squaredNorm();
}

void opengv::absolute_pose::modules::Epnp::compute_A_and_b_gauss_newton(
    const Eigen::Matrix<double, 6, 10>& L_6x10,
    const Eigen::Matrix<double, 6, 1>& Rho, const double betas[4],
    Eigen::Matrix<double, 6, 4>& A, Eigen::Matrix<double, 6, 1>& b) {
  for (int i = 0; i < 6; i++) {
    A(i, 0) = 2 * L_6x10(i, 0) * betas[0] + L_6x10(i, 1) * betas[1] +
              L_6x10(i, 3) * betas[2] + L_6x10(i, 6) * betas[3];
    A(i, 1) = L_6x10(i, 1) * betas[0] + 2 * L_6x10(i, 2) * betas[1] +
              L_6x10(i, 4) * betas[2] + L_6x10(i, 7) * betas[3];
    A(i, 2) = L_6x10(i, 3) * betas[0] + L_6x10(i, 4) * betas[1] +
              2 * L_6x10(i, 5) * betas[2] + L_6x10(i, 8) * betas[3];
    A(i, 3) = L_6x10(i, 6) * betas[0] + L_6x10(i, 7) * betas[1] +
              L_6x10(i, 8) * betas[2] + 2 * L_6x10(i, 9) * betas[3];

    b(i, 0) = Rho[i] - (L_6x10(i, 0) * betas[0] * betas[0] +
                        L_6x10(i, 1) * betas[0] * betas[1] +
                        L_6x10(i, 2) * betas[1] * betas[1] +
                        L_6x10(i, 3) * betas[0] * betas[2] +
                        L_6x10(i, 4) * betas[1] * betas[2] +
                        L_6x10(i, 5) * betas[2] * betas[2] +
                        L_6x10(i, 6) * betas[0] * betas[3] +
                        L_6x10(i, 7) * betas[1] * betas[3] +
                        L_6x10(i, 8) * betas[2] * betas[3] +
                        L_6x10(i, 9) * betas[3] * betas[3]);
  }
}

void opengv::absolute_pose::modules::Epnp::gauss_newton(
    const Eigen::Matrix<double, 6, 10>& L_6x10,
    const Eigen::Matrix<double, 6, 1>& Rho, double betas[4]) {
  const int iterations_number = 5;

  Eigen::Matrix<double, 6, 4> A;
  Eigen::Matrix<double, 6, 1> B;
  Eigen::Matrix<double, 4, 1> X;

  for (int k = 0; k < iterations_number; k++) {
    compute_A_and_b_gauss_newton(L_6x10, Rho, betas, A, B);
    qr_solve(A, B, X);

    for (int i = 0; i < 4; i++) betas[i] += X[i];
  }
}

void opengv::absolute_pose::modules::Epnp::qr_solve(
    Eigen::Matrix<double, 6, 4>& A_orig, Eigen::Matrix<double, 6, 1>& b,
    Eigen::Matrix<double, 4, 1>& X) {
  Eigen::Matrix<double, 4, 6> A = A_orig.transpose();

  double A1[4];
  double A2[4];

  const int nr = A_orig.rows();
  const int nc = A_orig.cols();

  double *pA = A.data(), *ppAkk = pA;
  for (int k = 0; k < nc; k++) {
    double *ppAik = ppAkk, eta = fabs(*ppAik);
    for (int i = k + 1; i < nr; i++) {
      double elt = fabs(*ppAik);
      if (eta < elt) eta = elt;
      ppAik += nc;
    }

    if (eta == 0) {
      A1[k] = A2[k] = 0.0;
      cerr << "God damnit, A is singular, this shouldn't happen." << endl;
      return;
    } else {
      double *ppAik = ppAkk, sum = 0.0, inv_eta = 1. / eta;
      for (int i = k; i < nr; i++) {
        *ppAik *= inv_eta;
        sum += *ppAik * *ppAik;
        ppAik += nc;
      }
      double sigma = sqrt(sum);
      if (*ppAkk < 0) sigma = -sigma;
      *ppAkk += sigma;
      A1[k] = sigma * *ppAkk;
      A2[k] = -eta * sigma;
      for (int j = k + 1; j < nc; j++) {
        double *ppAik = ppAkk, sum = 0;
        for (int i = k; i < nr; i++) {
          sum += *ppAik * ppAik[j - k];
          ppAik += nc;
        }
        double tau = sum / A1[k];
        ppAik      = ppAkk;
        for (int i = k; i < nr; i++) {
          ppAik[j - k] -= tau * *ppAik;
          ppAik += nc;
        }
      }
    }
    ppAkk += nc + 1;
  }

  // b <- Qt b
  double *ppAjj = pA, *pb = b.data();
  for (int j = 0; j < nc; j++) {
    double *ppAij = ppAjj, tau = 0;
    for (int i = j; i < nr; i++) {
      tau += *ppAij * pb[i];
      ppAij += nc;
    }
    tau /= A1[j];
    ppAij = ppAjj;
    for (int i = j; i < nr; i++) {
      pb[i] -= tau * *ppAij;
      ppAij += nc;
    }
    ppAjj += nc + 1;
  }

  // X = R-1 b
  double* pX = X.data();
  pX[nc - 1] = pb[nc - 1] / A2[nc - 1];
  for (int i = nc - 2; i >= 0; i--) {
    double *ppAij = pA + i * nc + (i + 1), sum = 0;

    for (int j = i + 1; j < nc; j++) {
      sum += *ppAij * pX[j];
      ppAij++;
    }
    pX[i] = (pb[i] - sum) / A2[i];
  }
}
