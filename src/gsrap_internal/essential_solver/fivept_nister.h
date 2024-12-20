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

#ifndef GSRAP_INTERNAL_ESSENTIAL_SOLVER_FIVEPT_NISTER_H
#define GSRAP_INTERNAL_ESSENTIAL_SOLVER_FIVEPT_NISTER_H

#include "Eigen/SVD"
#include "gsrap/macros.h"
#include "gsrap_internal/common.h"

GSRAP_IGNORE_STRICT_WARNING_PUSH

// opengv
#include "opengv/relative_pose/modules/main.hpp"

GSRAP_IGNORE_STRICT_WARNING_POP

namespace gsrap {

template <typename Iterator>
std::vector<Mat3d, Eigen::aligned_allocator<Mat3d>> FivePtNister(
    const eigen_vector<Eigen::Vector3d>& bearings1,
    const eigen_vector<Eigen::Vector3d>& bearings2,
    const std::vector<Iterator>& data) {
  Eigen::Matrix<double, 5, 9> Q;
  for (size_t i = 0; i < 5; ++i) {
    assert(0 <= data[i]->first && data[i]->first < bearings1.size());
    assert(0 <= data[i]->second && data[i]->second < bearings2.size());

    const Vec3d& f       = bearings1[data[i]->first];
    const Vec3d& f_prime = bearings2[data[i]->second];

    Mat3d tmp                  = f * f_prime.transpose();
    Q.row(static_cast<int>(i)) = Eigen::Matrix<double, 1, 9>(tmp.data());
  }
  Eigen::JacobiSVD<Eigen::Matrix<double, 5, 9>> SVD_Q(Q, Eigen::ComputeFullV);
  Eigen::Matrix<double, 9, 4> EE = SVD_Q.matrixV().block(0, 5, 9, 4);
  std::vector<Mat3d, Eigen::aligned_allocator<Mat3d>> Es;
  opengv::relative_pose::modules::fivept_nister_main(EE, Es);

  return Es;
}

}  // namespace gsrap

#endif  // GSRAP_INTERNAL_ESSENTIAL_SOLVER_FIVEPT_NISTER_H
