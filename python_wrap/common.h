#pragma once

#include "gsrap/macros.h"

GSRAP_IGNORE_STRICT_WARNING_PUSH
#include "Eigen/Core"
GSRAP_IGNORE_STRICT_WARNING_POP

namespace gsrap {

using RMatrix3d = Eigen::Matrix<double, 3, 3, Eigen::RowMajorBit>;
using RMatrixXd =
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajorBit>;

using RMatrixXu32 =
    Eigen::Matrix<uint32_t, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajorBit>;

using Vec3d = Eigen::Vector3d;

template <typename T>
using eigen_vector = std::vector<T
#ifdef GSRAP_LESS_EIGEN_3_4
                                 ,
                                 Eigen::aligned_allocator<T>
#endif  // GSRAP_LESS_EIGEN_3_4
                                 >;

}  // namespace gsrap
