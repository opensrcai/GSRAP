#pragma once

#include "gsrap/macros.h"

GSRAP_IGNORE_STRICT_WARNING_PUSH

#include "Eigen/Core"

GSRAP_IGNORE_STRICT_WARNING_POP

namespace gsrap::example {

using Mat3d = Eigen::Matrix3d;
using Vec2d = Eigen::Vector2d;
using Vec3d = Eigen::Vector3d;

template <typename T>
using eigen_vector = std::vector<T
#ifdef GSRAP_LESS_EIGEN_3_4
                                 ,
                                 Eigen::aligned_allocator<T>
#endif  // GSRAP_LESS_EIGEN_3_4
                                 >;

}  // namespace gsrap::example
