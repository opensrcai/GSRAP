#pragma once

#include "common/type.h"

namespace gsrap::example {
struct Camera {
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  Mat3d rotation_cw;     // world to camera
  Vec3d transltaion_cw;  // world to camera
};

}  // namespace gsrap::example
