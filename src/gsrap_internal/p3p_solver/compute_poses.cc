#include "gsrap_internal/p3p_solver/compute_poses.h"

#include <math.h>

namespace gsrap {

void compute_real_solution_of_quartic_equation(
    double b, const double c, const double d, const double e,
    double solutions[4], size_t *num_solutions, const double D_thr,
    const double biquadratic_equation_thr, const double multiple_solution_thr);

// from OpenCV start ///////////////////////////////////////////////////////////

// OpenCV LICENSE
/*
                                 Apache License
                           Version 2.0, January 2004
                        http://www.apache.org/licenses/

   TERMS AND CONDITIONS FOR USE, REPRODUCTION, AND DISTRIBUTION

   1. Definitions.

      "License" shall mean the terms and conditions for use, reproduction,
      and distribution as defined by Sections 1 through 9 of this document.

      "Licensor" shall mean the copyright owner or entity authorized by
      the copyright owner that is granting the License.

      "Legal Entity" shall mean the union of the acting entity and all
      other entities that control, are controlled by, or are under common
      control with that entity. For the purposes of this definition,
      "control" means (i) the power, direct or indirect, to cause the
      direction or management of such entity, whether by contract or
      otherwise, or (ii) ownership of fifty percent (50%) or more of the
      outstanding shares, or (iii) beneficial ownership of such entity.

      "You" (or "Your") shall mean an individual or Legal Entity
      exercising permissions granted by this License.

      "Source" form shall mean the preferred form for making modifications,
      including but not limited to software source code, documentation
      source, and configuration files.

      "Object" form shall mean any form resulting from mechanical
      transformation or translation of a Source form, including but
      not limited to compiled object code, generated documentation,
      and conversions to other media types.

      "Work" shall mean the work of authorship, whether in Source or
      Object form, made available under the License, as indicated by a
      copyright notice that is included in or attached to the work
      (an example is provided in the Appendix below).

      "Derivative Works" shall mean any work, whether in Source or Object
      form, that is based on (or derived from) the Work and for which the
      editorial revisions, annotations, elaborations, or other modifications
      represent, as a whole, an original work of authorship. For the purposes
      of this License, Derivative Works shall not include works that remain
      separable from, or merely link (or bind by name) to the interfaces of,
      the Work and Derivative Works thereof.

      "Contribution" shall mean any work of authorship, including
      the original version of the Work and any modifications or additions
      to that Work or Derivative Works thereof, that is intentionally
      submitted to Licensor for inclusion in the Work by the copyright owner
      or by an individual or Legal Entity authorized to submit on behalf of
      the copyright owner. For the purposes of this definition, "submitted"
      means any form of electronic, verbal, or written communication sent
      to the Licensor or its representatives, including but not limited to
      communication on electronic mailing lists, source code control systems,
      and issue tracking systems that are managed by, or on behalf of, the
      Licensor for the purpose of discussing and improving the Work, but
      excluding communication that is conspicuously marked or otherwise
      designated in writing by the copyright owner as "Not a Contribution."

      "Contributor" shall mean Licensor and any individual or Legal Entity
      on behalf of whom a Contribution has been received by Licensor and
      subsequently incorporated within the Work.

   2. Grant of Copyright License. Subject to the terms and conditions of
      this License, each Contributor hereby grants to You a perpetual,
      worldwide, non-exclusive, no-charge, royalty-free, irrevocable
      copyright license to reproduce, prepare Derivative Works of,
      publicly display, publicly perform, sublicense, and distribute the
      Work and such Derivative Works in Source or Object form.

   3. Grant of Patent License. Subject to the terms and conditions of
      this License, each Contributor hereby grants to You a perpetual,
      worldwide, non-exclusive, no-charge, royalty-free, irrevocable
      (except as stated in this section) patent license to make, have made,
      use, offer to sell, sell, import, and otherwise transfer the Work,
      where such license applies only to those patent claims licensable
      by such Contributor that are necessarily infringed by their
      Contribution(s) alone or by combination of their Contribution(s)
      with the Work to which such Contribution(s) was submitted. If You
      institute patent litigation against any entity (including a
      cross-claim or counterclaim in a lawsuit) alleging that the Work
      or a Contribution incorporated within the Work constitutes direct
      or contributory patent infringement, then any patent licenses
      granted to You under this License for that Work shall terminate
      as of the date such litigation is filed.

   4. Redistribution. You may reproduce and distribute copies of the
      Work or Derivative Works thereof in any medium, with or without
      modifications, and in Source or Object form, provided that You
      meet the following conditions:

      (a) You must give any other recipients of the Work or
          Derivative Works a copy of this License; and

      (b) You must cause any modified files to carry prominent notices
          stating that You changed the files; and

      (c) You must retain, in the Source form of any Derivative Works
          that You distribute, all copyright, patent, trademark, and
          attribution notices from the Source form of the Work,
          excluding those notices that do not pertain to any part of
          the Derivative Works; and

      (d) If the Work includes a "NOTICE" text file as part of its
          distribution, then any Derivative Works that You distribute must
          include a readable copy of the attribution notices contained
          within such NOTICE file, excluding those notices that do not
          pertain to any part of the Derivative Works, in at least one
          of the following places: within a NOTICE text file distributed
          as part of the Derivative Works; within the Source form or
          documentation, if provided along with the Derivative Works; or,
          within a display generated by the Derivative Works, if and
          wherever such third-party notices normally appear. The contents
          of the NOTICE file are for informational purposes only and
          do not modify the License. You may add Your own attribution
          notices within Derivative Works that You distribute, alongside
          or as an addendum to the NOTICE text from the Work, provided
          that such additional attribution notices cannot be construed
          as modifying the License.

      You may add Your own copyright statement to Your modifications and
      may provide additional or different license terms and conditions
      for use, reproduction, or distribution of Your modifications, or
      for any such Derivative Works as a whole, provided Your use,
      reproduction, and distribution of the Work otherwise complies with
      the conditions stated in this License.

   5. Submission of Contributions. Unless You explicitly state otherwise,
      any Contribution intentionally submitted for inclusion in the Work
      by You to the Licensor shall be under the terms and conditions of
      this License, without any additional terms or conditions.
      Notwithstanding the above, nothing herein shall supersede or modify
      the terms of any separate license agreement you may have executed
      with Licensor regarding such Contributions.

   6. Trademarks. This License does not grant permission to use the trade
      names, trademarks, service marks, or product names of the Licensor,
      except as required for reasonable and customary use in describing the
      origin of the Work and reproducing the content of the NOTICE file.

   7. Disclaimer of Warranty. Unless required by applicable law or
      agreed to in writing, Licensor provides the Work (and each
      Contributor provides its Contributions) on an "AS IS" BASIS,
      WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
      implied, including, without limitation, any warranties or conditions
      of TITLE, NON-INFRINGEMENT, MERCHANTABILITY, or FITNESS FOR A
      PARTICULAR PURPOSE. You are solely responsible for determining the
      appropriateness of using or redistributing the Work and assume any
      risks associated with Your exercise of permissions under this License.

   8. Limitation of Liability. In no event and under no legal theory,
      whether in tort (including negligence), contract, or otherwise,
      unless required by applicable law (such as deliberate and grossly
      negligent acts) or agreed to in writing, shall any Contributor be
      liable to You for damages, including any direct, indirect, special,
      incidental, or consequential damages of any character arising as a
      result of this License or out of the use or inability to use the
      Work (including but not limited to damages for loss of goodwill,
      work stoppage, computer failure or malfunction, or any and all
      other commercial damages or losses), even if such Contributor
      has been advised of the possibility of such damages.

   9. Accepting Warranty or Additional Liability. While redistributing
      the Work or Derivative Works thereof, You may choose to offer,
      and charge a fee for, acceptance of support, warranty, indemnity,
      or other liability obligations and/or rights consistent with this
      License. However, in accepting such obligations, You may act only
      on Your own behalf and on Your sole responsibility, not on behalf
      of any other Contributor, and only if You agree to indemnify,
      defend, and hold each Contributor harmless for any liability
      incurred by, or claims asserted against, such Contributor by reason
      of your accepting any such warranty or additional liability.

   END OF TERMS AND CONDITIONS

   APPENDIX: How to apply the Apache License to your work.

      To apply the Apache License to your work, attach the following
      boilerplate notice, with the fields enclosed by brackets "[]"
      replaced with your own identifying information. (Don't include
      the brackets!)  The text should be enclosed in the appropriate
      comment syntax for the file format. We also recommend that a
      file or class name and description of purpose be included on the
      same "printed page" as the copyright notice for easier
      identification within third-party archives.

   Copyright [yyyy] [name of copyright owner]

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
   */


static void vect_cross(const double *a, const double *b, double *result) {
  result[0] = a[1] * b[2] - a[2] * b[1];
  result[1] = -(a[0] * b[2] - a[2] * b[0]);
  result[2] = a[0] * b[1] - a[1] * b[0];
}

static double vect_dot(const double *a, const double *b) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

static double vect_norm(const double *a) {
  return sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
}

inline void vect_scale(const double s, const double *a, double *result) {
  result[0] = a[0] * s;
  result[1] = a[1] * s;
  result[2] = a[2] * s;
}

inline void vect_sub(const double *a, const double *b, double *result) {
  result[0] = a[0] - b[0];
  result[1] = a[1] - b[1];
  result[2] = a[2] - b[2];
}

inline void vect_divide(const double *a, const double d, double *result) {
  result[0] = a[0] / d;
  result[1] = a[1] / d;
  result[2] = a[2] / d;
}

inline void mat_mult(const double a[3][3], const double b[3][3],
                     double result[3][3]) {
  result[0][0] = a[0][0] * b[0][0] + a[0][1] * b[1][0] + a[0][2] * b[2][0];
  result[0][1] = a[0][0] * b[0][1] + a[0][1] * b[1][1] + a[0][2] * b[2][1];
  result[0][2] = a[0][0] * b[0][2] + a[0][1] * b[1][2] + a[0][2] * b[2][2];

  result[1][0] = a[1][0] * b[0][0] + a[1][1] * b[1][0] + a[1][2] * b[2][0];
  result[1][1] = a[1][0] * b[0][1] + a[1][1] * b[1][1] + a[1][2] * b[2][1];
  result[1][2] = a[1][0] * b[0][2] + a[1][1] * b[1][2] + a[1][2] * b[2][2];

  result[2][0] = a[2][0] * b[0][0] + a[2][1] * b[1][0] + a[2][2] * b[2][0];
  result[2][1] = a[2][0] * b[0][1] + a[2][1] * b[1][1] + a[2][2] * b[2][1];
  result[2][2] = a[2][0] * b[0][2] + a[2][1] * b[1][2] + a[2][2] * b[2][2];
}

// This algorithm is from "Tong Ke, Stergios Roumeliotis, An Efficient Algebraic
// Solution to the Perspective-Three-Point Problem" (Accepted by CVPR 2017) See
// https://arxiv.org/pdf/1701.08237.pdf featureVectors: The 3 bearing
// measurements (normalized) stored as column vectors worldPoints: The positions
// of the 3 feature points stored as column vectors solutionsR: 4 possible
// solutions of rotation matrix of the world w.r.t the camera frame solutionsT:
// 4 possible solutions of translation of the world origin w.r.t the camera
// frame
int computePoses(const double featureVectors[3][4],
                 const double worldPoints[3][4], double solutionsR[4][3][3],
                 double solutionsT[4][3]) {
  // world point vectors
  double w1[3] = {worldPoints[0][0], worldPoints[1][0], worldPoints[2][0]};
  double w2[3] = {worldPoints[0][1], worldPoints[1][1], worldPoints[2][1]};
  double w3[3] = {worldPoints[0][2], worldPoints[1][2], worldPoints[2][2]};
  // k1
  double u0[3];
  vect_sub(w1, w2, u0);

  double nu0 = vect_norm(u0);
  double k1[3];
  vect_divide(u0, nu0, k1);
  // bi
  double b1[3] = {featureVectors[0][0], featureVectors[1][0],
                  featureVectors[2][0]};
  double b2[3] = {featureVectors[0][1], featureVectors[1][1],
                  featureVectors[2][1]};
  double b3[3] = {featureVectors[0][2], featureVectors[1][2],
                  featureVectors[2][2]};
  // k3,tz
  double k3[3];
  vect_cross(b1, b2, k3);
  double nk3 = vect_norm(k3);
  vect_divide(k3, nk3, k3);

  double tz[3];
  vect_cross(b1, k3, tz);
  // ui,vi
  double v1[3];
  vect_cross(b1, b3, v1);
  double v2[3];
  vect_cross(b2, b3, v2);

  double u1[3];
  vect_sub(w1, w3, u1);
  // coefficients related terms
  double u1k1 = vect_dot(u1, k1);
  double k3b3 = vect_dot(k3, b3);
  // f1i
  double f11 = k3b3;
  double f13 = vect_dot(k3, v1);
  double f15 = -u1k1 * f11;
  // delta
  double nl[3];
  vect_cross(u1, k1, nl);
  double delta = vect_norm(nl);
  vect_divide(nl, delta, nl);
  f11 *= delta;
  f13 *= delta;
  // f2i
  double u2k1 = u1k1 - nu0;
  double f21  = vect_dot(tz, v2);
  double f22  = nk3 * k3b3;
  double f23  = vect_dot(k3, v2);
  double f24  = u2k1 * f22;
  double f25  = -u2k1 * f21;
  f21 *= delta;
  f22 *= delta;
  f23 *= delta;
  double g1        = f13 * f22;
  double g2        = f13 * f25 - f15 * f23;
  double g3        = f11 * f23 - f13 * f21;
  double g4        = -f13 * f24;
  double g5        = f11 * f22;
  double g6        = f11 * f25 - f15 * f21;
  double g7        = -f15 * f24;
  double coeffs[5] = {
      g5 * g5 + g1 * g1 + g3 * g3, 2 * (g5 * g6 + g1 * g2 + g3 * g4),
      g6 * g6 + 2 * g5 * g7 + g2 * g2 + g4 * g4 - g1 * g1 - g3 * g3,
      2 * (g6 * g7 - g1 * g2 - g3 * g4), g7 * g7 - g2 * g2 - g4 * g4};

  double s[4];
  size_t num_solutions_of_quartic_equation = 0;
  const double inv_coeff0                  = 1.0 / coeffs[0];
  compute_real_solution_of_quartic_equation(
      coeffs[1] * inv_coeff0, coeffs[2] * inv_coeff0, coeffs[3] * inv_coeff0,
      coeffs[4] * inv_coeff0, s, &num_solutions_of_quartic_equation, 1e-14,
      1e-14, 1e-14);

  // TODO(kyawakyawa): Refine s[]

  double temp[3];
  vect_cross(k1, nl, temp);

  double Ck1nl[3][3] = {{k1[0], nl[0], temp[0]},
                        {k1[1], nl[1], temp[1]},
                        {k1[2], nl[2], temp[2]}};

  double Cb1k3tzT[3][3] = {
      {b1[0], b1[1], b1[2]}, {k3[0], k3[1], k3[2]}, {tz[0], tz[1], tz[2]}};

  double b3p[3];
  vect_scale((delta / k3b3), b3, b3p);

  int nb_solutions = 0;
  for (int i = 0; i < int(num_solutions_of_quartic_equation); ++i) {
    double ctheta1p = s[i];
    if (abs(ctheta1p) > 1) {
      continue;
    }
    double stheta1p = sqrt(1 - ctheta1p * ctheta1p);
    stheta1p        = (k3b3 > 0) ? stheta1p : -stheta1p;
    double ctheta3  = g1 * ctheta1p + g2;
    double stheta3  = g3 * ctheta1p + g4;
    double ntheta3  = stheta1p / ((g5 * ctheta1p + g6) * ctheta1p + g7);
    ctheta3 *= ntheta3;
    stheta3 *= ntheta3;

    double C13[3][3] = {{ctheta3, 0, -stheta3},
                        {stheta1p * stheta3, ctheta1p, stheta1p * ctheta3},
                        {ctheta1p * stheta3, -stheta1p, ctheta1p * ctheta3}};

    double temp_matrix[3][3];
    double R[3][3];
    mat_mult(Ck1nl, C13, temp_matrix);
    mat_mult(temp_matrix, Cb1k3tzT, R);

    // R' * p3
    double rp3[3] = {w3[0] * R[0][0] + w3[1] * R[1][0] + w3[2] * R[2][0],
                     w3[0] * R[0][1] + w3[1] * R[1][1] + w3[2] * R[2][1],
                     w3[0] * R[0][2] + w3[1] * R[1][2] + w3[2] * R[2][2]};

    double pxstheta1p[3];
    vect_scale(stheta1p, b3p, pxstheta1p);

    vect_sub(pxstheta1p, rp3, solutionsT[nb_solutions]);

    solutionsR[nb_solutions][0][0] = R[0][0];
    solutionsR[nb_solutions][1][0] = R[0][1];
    solutionsR[nb_solutions][2][0] = R[0][2];
    solutionsR[nb_solutions][0][1] = R[1][0];
    solutionsR[nb_solutions][1][1] = R[1][1];
    solutionsR[nb_solutions][2][1] = R[1][2];
    solutionsR[nb_solutions][0][2] = R[2][0];
    solutionsR[nb_solutions][1][2] = R[2][1];
    solutionsR[nb_solutions][2][2] = R[2][2];

    nb_solutions++;
  }
  return nb_solutions;
}

// from OpenCV end /////////////////////////////////////////////////////////////

}  // namespace gsrap
