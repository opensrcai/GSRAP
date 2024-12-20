#pragma once

#include <stdint.h>

#include <random>
#include <unordered_set>
#include <utility>

#include "common/camera.h"
#include "common/constant.h"
#include "common/type.h"
#include "gsrap/macros.h"

GSRAP_IGNORE_STRICT_WARNING_PUSH

#include "Eigen/Core"
#include "Eigen/Geometry"

GSRAP_IGNORE_STRICT_WARNING_POP
namespace gsrap::example {

// This function calculates angle of rotaion.
inline double AxisAngleTheta(const Mat3d& R) {
  const double cos_theta = (R.trace() - 1.0) / 2.0;
  return acos(cos_theta);  // FIXME:(kyawakyawa): Range check
}

// This function creates orthonormal basis its z axis is n.
inline std::pair</*x*/ Vec3d, /*y*/ Vec3d> CreateOnb(const Vec3d& n) {
  const double sign = (n.z() >= 0.0 ? 1.0 : -1.0);
  const double a    = -1.0 / (sign + n.z());
  const double b    = n.x() * n.y() * a;
  return std::make_pair(
      Vec3d(1.0 + sign * n.x() * n.x() * a, sign * b, -sign * n.x()),
      Vec3d(b, sign + n.y() * n.y() * a, -n.y()));
}

// This function samples a 3D point in the unit sphere uniformly.
inline Vec3d SampleInUnitSphere(double u, double v, double w) {
  u               = pow(u, 1.0 / 3.0);
  const double _v = sqrt(1 - v * v);
  w               = 2.0 * kPi * w;

  Vec3d ret;
  ret.x() = u * _v * cos(w);
  ret.y() = u * _v * sin(w);
  ret.z() = u * v;

  return ret;
}

// This function creates spherical point cloud randomly.
inline eigen_vector<Vec3d> CreatesphericalPointCloud(const Vec3d& center,
                                                     const double radius,
                                                     const uint32_t num_point,
                                                     std::mt19937* engine) {
  eigen_vector<Vec3d> ret(num_point);

  std::uniform_real_distribution<double> dist(0.0, 1.0);

  for (size_t i = 0; i < num_point; ++i) {
    Vec3d p = SampleInUnitSphere(dist(*engine), dist(*engine), dist(*engine));

    p = radius * p + center;

    ret[i] = p;
  }
  return ret;
}

// This function creates cuboid point cloud randomly.
inline eigen_vector<Vec3d> CreateCuboidlPointCloud(
    const Vec3d& center, const double x, const double y, const double z,
    const uint32_t num_point, std::mt19937* engine) {
  eigen_vector<Vec3d> ret(num_point);

  std::uniform_real_distribution<double> dist(0.0, 1.0);

  for (size_t i = 0; i < num_point; ++i) {
    const double u = dist(*engine) - 0.5;
    const double v = dist(*engine) - 0.5;
    const double w = dist(*engine) - 0.5;

    ret[i] = Vec3d(center[0] + u * x, center[1] + v * y, center[2] + w * z);
  }
  return ret;
}

// This function converts the points to the bearing vectors.
inline eigen_vector<Vec3d> ToBearingVectors(const eigen_vector<Vec3d>& points,
                                            const Camera& camera) {
  const size_t num_points = points.size();
  eigen_vector<Vec3d> bearings;
  bearings.reserve(num_points);
  for (uint32_t m_id = 0; m_id < num_points; ++m_id) {
    const Vec3d& point = points.at(m_id);

    const Vec3d point_c = camera.rotation_cw * point + camera.transltaion_cw;

    bearings.emplace_back(point_c.normalized());
  }

  return bearings;
}
// This function converts the points to the bearing vectors.
inline std::pair<eigen_vector<Vec3d>, eigen_vector<Vec3d>> ToBearingVectors(
    const eigen_vector<Vec3d>& points, const Camera& camera1,
    const Camera& camera2) {
  return std::make_pair(ToBearingVectors(points, camera1),
                        ToBearingVectors(points, camera2));
}

// This function create matches.
inline std::vector<std::pair<uint32_t, uint32_t>> CreateMatches(
    const size_t num_points, const double inlier_ratio, std::mt19937* engine) {
  const size_t num_matches = num_points;
  const auto num_miss_matches =
      uint32_t(std::round(double(num_points) * (1.0 - inlier_ratio)));

  std::unordered_set<int64_t> miss_match_ids;
  std::uniform_int_distribution<int64_t> dist(0, int64_t(num_matches) - 1);

  while (miss_match_ids.size() < num_miss_matches) {
    miss_match_ids.emplace(dist(*engine));
  }

  const std::vector<size_t> ids1(miss_match_ids.begin(), miss_match_ids.end());
  std::vector<size_t> ids2;
  while (true) {
    ids2 = ids1;
    std::shuffle(ids2.begin(), ids2.end(), *engine);
    bool flag = true;
    for (size_t i = 0; i < ids1.size(); ++i) {
      if (ids1[i] == ids2[i]) {
        flag = false;
        break;
      }
    }
    if (flag) {
      break;
    }
  }

  std::vector<std::pair<uint32_t, uint32_t>> matches;
  matches.reserve(num_matches);

  for (size_t i = 0; i < ids1.size(); ++i) {
    matches.emplace_back(uint32_t(ids1[i]), uint32_t(ids2[i]));
  }

  for (uint32_t i = 0; i < num_matches; ++i) {
    if (miss_match_ids.count(i) == 0) {
      matches.emplace_back(uint32_t(i), uint32_t(i));
    }
  }

  return matches;
}
// This function rotates bearing vector.
inline Vec3d RotateBearingVector(Vec3d n, const double theta /* rad */,
                                 const double phi /* rad */) {
  n.normalize();
  auto [x, y] = CreateOnb(n);

  const Vec3d axis = cos(phi) * x + sin(phi) * y;

  Eigen::Affine3d rot = {};
  rot                 = Eigen::AngleAxisd(theta, axis);

  n = rot * n;

  return n;
}

inline void RotateBearingVectors(const double angle /* rad */,
                                 eigen_vector<Vec3d>* ns,
                                 std::mt19937* engine) {
  std::uniform_real_distribution<double> real_dist(0.0, 1.0);

  // More than 99% of vectors have an error less than angle
  std::normal_distribution<double> normal_dist(0.0, angle / 3.0);

  for (Vec3d& v : *ns) {
    const double phi   = 2.0 * kPi * real_dist(*engine);
    const double theta = normal_dist(*engine);

    v = RotateBearingVector(v, theta, phi);
  }
}

inline void RotateSampledBearingVectors(const double angle, double rate,
                                        eigen_vector<Vec3d>* ns,
                                        std::mt19937* engine) {
  rate                  = std::clamp(rate, 0.0, 1.0);
  const auto num_sample = size_t(std::round(double(ns->size()) * rate));

  std::unordered_set<int64_t> sampled_ids;
  std::uniform_int_distribution<int64_t> dist(0, int64_t(ns->size()) - 1);

  while (sampled_ids.size() < num_sample) {
    sampled_ids.emplace(dist(*engine));
  }

  std::uniform_real_distribution<double> real_dist(0.0, 1.0);

  std::normal_distribution<double> normal_dist(0.0, angle * 0.1);

  for (const auto sampled_id : sampled_ids) {
    {
      const double phi   = 2.0 * kPi * real_dist(*engine);
      const double theta = angle + normal_dist(*engine);
      ns->at(size_t(sampled_id)) =
          RotateBearingVector(ns->at(size_t(sampled_id)), theta, phi);
    }
  }
}

inline void RotateSampledBearingVectors(const double angle, double rate,
                                        eigen_vector<Vec3d>* n1s,
                                        eigen_vector<Vec3d>* n2s,
                                        std::mt19937* engine) {
  assert(n1s->size() == n2s->size());
  if (n1s->size() != n2s->size()) {
    return;
  }
  rate                  = std::clamp(rate, 0.0, 1.0);
  const auto num_sample = size_t(std::round(double(n1s->size()) * rate));

  std::unordered_set<int64_t> sampled_ids;
  std::uniform_int_distribution<int64_t> dist(0, int64_t(n1s->size()) - 1);

  while (sampled_ids.size() < num_sample) {
    sampled_ids.emplace(dist(*engine));
  }

  std::uniform_real_distribution<double> real_dist(0.0, 1.0);

  std::normal_distribution<double> normal_dist(0.0, angle * 0.1);

  for (const auto sampled_id : sampled_ids) {
    {
      const double phi   = 2.0 * kPi * real_dist(*engine);
      const double theta = angle + normal_dist(*engine);
      n1s->at(size_t(sampled_id)) =
          RotateBearingVector(n1s->at(size_t(sampled_id)), theta, phi);
    }
    {
      const double phi   = 2.0 * kPi * real_dist(*engine);
      const double theta = angle + normal_dist(*engine);
      n2s->at(size_t(sampled_id)) =
          RotateBearingVector(n2s->at(size_t(sampled_id)), theta, phi);
    }
  }
}

inline void AddNoiseToPoints(const double slide, const double sigma,
                             double rate, eigen_vector<Vec3d>* points,
                             std::mt19937* engine) {
  std::normal_distribution<double> dist_nor(0.0, sigma);
  std::uniform_real_distribution<double> dist_uni(0.0, 1.0);

  rate                  = std::clamp(rate, 0.0, 1.0);
  const auto num_sample = size_t(std::round(double(points->size()) * rate));

  std::unordered_set<int64_t> sampled_ids;
  std::uniform_int_distribution<int64_t> dist_int(0,
                                                  int64_t(points->size()) - 1);

  while (sampled_ids.size() < num_sample) {
    sampled_ids.emplace(dist_int(*engine));
  }

  for (const auto sampled_id : sampled_ids) {
    Vec3d& p      = points->at(size_t(sampled_id));
    const Vec3d b = SampleInUnitSphere(dist_uni(*engine), dist_uni(*engine),
                                       dist_uni(*engine))
                        .normalized();

    const double d = slide + dist_nor(*engine);
    p += d * b;
  }
}

}  // namespace gsrap::example
