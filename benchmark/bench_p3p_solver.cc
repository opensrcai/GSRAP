#include <cstdint>
#include <ctime>

//
#include <algorithm>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <ios>
#include <iostream>
#include <limits>
#include <numeric>
#include <random>
#include <sstream>
#include <string>
// GSRAP
#include "gsrap/macros.h"
#include "gsrap/p3p_solver.h"

// common
#include "common/camera.h"
#include "common/type.h"
#include "common/utils.h"

GSRAP_IGNORE_STRICT_WARNING_PUSH

// Eigen
#include "Eigen/Core"
#include "Eigen/Dense"
#include "Eigen/Geometry"

// Cerial
#include "cereal/archives/json.hpp"
#include "cereal/types/string.hpp"
#include "cereal/types/vector.hpp"

// cxxopts
#include "cxxopts.hpp"

GSRAP_IGNORE_STRICT_WARNING_POP

using namespace gsrap::example;  // NOLINT

namespace gsrap::benchmark {

struct Result {
  std::vector<std::vector<std::vector<double>>> rotation;     // * x 3 x 3
  std::vector<std::vector<double>> translation;               // * x 3
  std::vector<std::vector<std::vector<double>>> gt_rotation;  // * x 3 x 3
  std::vector<std::vector<double>> gt_translation;            // * x 3
  std::vector<double> fx;
  std::vector<double> fy;
  std::vector<double> time_ns;
  std::vector<double> depth_mean;
  std::vector<double> depth_std_deviation;
  std::vector<uint32_t> number_of_solutions;

  template <class Archive>
  void serialize(Archive& ar) const {
    ar(  // NOLINT
        CEREAL_NVP(rotation), CEREAL_NVP(translation), CEREAL_NVP(gt_rotation),
        CEREAL_NVP(gt_translation), CEREAL_NVP(fx), CEREAL_NVP(fy),
        CEREAL_NVP(time_ns), CEREAL_NVP(depth_mean),
        CEREAL_NVP(depth_std_deviation), CEREAL_NVP(number_of_solutions));
  }
};

struct Output {
  std::string name;
  std::string date;

  std::string desc;

  uint32_t width;
  uint32_t height;

  uint64_t num_bench_itr;
  uint64_t num_found_solutions;

  Result result;

  template <class Archive>
  void serialize(Archive& ar) const {
    ar(  // NOLINT
        CEREAL_NVP(name), CEREAL_NVP(date), CEREAL_NVP(desc), CEREAL_NVP(width),
        CEREAL_NVP(height), CEREAL_NVP(num_bench_itr),
        CEREAL_NVP(num_found_solutions), CEREAL_NVP(result));
  }
};

static std::string FillZero(int64_t value, int32_t n) {
  std::stringstream ss;
  ss << std::setfill('0') << std::right << std::setw(n) << value;

  return ss.str();
}

static std::string GetUtcTime() {
  time_t t         = 0;
  struct tm utc_tm = {};
  time(&t);
  utc_tm = *gmtime(&t);  // [NOTE] No thread safe NOLINT

  const std::string year = FillZero(utc_tm.tm_year + 1900, 4);
  const std::string mon  = FillZero(utc_tm.tm_mon + 1, 2);
  const std::string day  = FillZero(utc_tm.tm_mday, 2);
  const std::string hour = FillZero(utc_tm.tm_hour, 2);
  const std::string min  = FillZero(utc_tm.tm_min, 2);
  const std::string sec  = FillZero(utc_tm.tm_sec, 2);

  return year + "-" + mon + "-" + day + "T" + hour + ":" + min + ":" + sec +
         "+00:00";
}

using Pair     = std::pair<uint32_t, uint32_t>;
using Iterator = std::vector<Pair>::const_iterator;

static Camera CreateCameras(std::mt19937* engine) {
  Camera camera = {};

  const double xm = -100;
  const double xM = 100;

  const double ym = -100;
  const double yM = 100;

  const double zm = -100;
  const double zM = 100;

  std::uniform_real_distribution<> dist(0.0, 1.0);

  const double tx = dist(*engine);
  const double ty = dist(*engine);
  const double tz = dist(*engine);

  const Vec3d trans_cw((1.0 - tx) * xm + tx * xM, (1.0 - ty) * ym + ty * yM,
                       (1.0 - tz) * zm + tz * zM);  // camera to world

  const Vec3d axis =
      SampleInUnitSphere(dist(*engine), dist(*engine), dist(*engine))
          .normalized();

  const double theta = dist(*engine) * (2.0 * M_PI);

  Mat3d rotation_wc = {};  // camera to world
  rotation_wc       = Eigen::AngleAxisd(theta, axis);
#if 0
  std::cout << rotation_wc << std::endl;
  std::cout << "det: " << rotation_wc.determinant() << std::endl;
  std::cout << "cross check: "
            << (rotation_wc.col(0).cross(rotation_wc.col(1)) -
                rotation_wc.col(2))
                   .norm()
            << std::endl;
  std::cout << std::endl;
#endif

  camera.rotation_cw    = rotation_wc.transpose();
  camera.transltaion_cw = -camera.rotation_cw * trans_cw;

  return camera;
}

static std::tuple<eigen_vector<Vec3d> /*Points*/,
                  eigen_vector<Vec2d> /*Keypts*/, double /*fx*/, double /*fy*/>
CreatePointCloud(const size_t num_sample, const uint32_t width,
                 const uint32_t height, const double fov_min /*horizontal*/,
                 const double fov_max /*horizontal*/, const double depth_min,
                 const double depth_max,
                 const Mat3d& rotation_wc /*camera to world */,
                 const Vec3d& translation_wc /*camera to world*/,
                 std::mt19937* engine) {
  std::uniform_real_distribution<double> uniform_rng(0.0, 1.0);

  const double tmp   = uniform_rng(*engine);
  const double fov_x = (1.0 - tmp) * fov_min + tmp * fov_max;
#if 0
  const double w = width;
  const double h = height;
  const double fov_y = 2.0 * atan(h / w * tan(fov_x * 0.5));
#endif

  const double fx = 0.5 * width / tan(fov_x * 0.5);
#if 0
  const double fy = 0.5 * height / tan(fov_y * 0.5);  // = fx
#else
  const double fy = fx;
#endif

  const double inv_fx = 1.0 / fx;
  const double inv_fy = 1.0 / fy;

  eigen_vector<Vec3d> points;
  eigen_vector<Vec2d> keypts;

  for (size_t itr = 0; itr < num_sample; ++itr) {
    // cx == cy == 0
    const double u = (uniform_rng(*engine) - 0.5) * width;
    const double v = (uniform_rng(*engine) - 0.5) * height;

    const double z = [&]() {
      const double t = uniform_rng(*engine);
      return (1.0 - t) * depth_min + t * depth_max;
    }();

    const double x = inv_fx * u * z;
    const double y = inv_fy * v * z;

    const Vec3d pc(x, y, z);                             // camera coords
    const Vec3d pw = rotation_wc * pc + translation_wc;  // world corrds

    points.emplace_back(pw);
    keypts.emplace_back(u, v);
  }

  return {points, keypts, fx, fy};
}

static eigen_vector<Vec3d> ToBearingVectors(const eigen_vector<Vec2d>& keypts,
                                            const double fx, const double fy,
                                            const bool add_noise,
                                            std::mt19937* engine = nullptr,
                                            const double sigma_x = 1.0,
                                            const double sigma_y = 1.0) {
  eigen_vector<Vec3d> bearings;

  const double inv_fx = 1.0 / fx;
  const double inv_fy = 1.0 / fy;

  std::normal_distribution<> distx(0.0, sigma_x);
  std::normal_distribution<> disty(0.0, sigma_y);

  for (const auto& keypt : keypts) {
    double u = keypt.x();
    double v = keypt.y();
    if (add_noise) {
      u += distx(*engine);
      v += disty(*engine);
    }
    const double x = inv_fx * u;
    const double y = inv_fy * v;
    const double z = 1.0;

    bearings.emplace_back(Vec3d(x, y, z).normalized());
  }

  return bearings;
}

static std::vector<std::vector<double>> Mat3d2StdVector(const Mat3d& mat) {
  std::vector<std::vector<double>> ret(3, std::vector<double>(3, 0));

  ret[0][0] = mat(0, 0);
  ret[0][1] = mat(0, 1);
  ret[0][2] = mat(0, 2);

  ret[1][0] = mat(1, 0);
  ret[1][1] = mat(1, 1);
  ret[1][2] = mat(1, 2);

  ret[2][0] = mat(2, 0);
  ret[2][1] = mat(2, 1);
  ret[2][2] = mat(2, 2);

  return ret;
}

static std::vector<double> Vec3d2StdVector(const Vec3d& vec) {
  std::vector<double> ret(3, 0);
  ret[0] = vec.x();
  ret[1] = vec.y();
  ret[2] = vec.z();

  return ret;
}

static void SaveOutputAsJson(const std::string& filepath,
                             const Output& output) {
  std::ofstream os(filepath, std::ios::binary);
  cereal::JSONOutputArchive archive(os);

  output.serialize(archive);  // NOLINT
}

static void BenchP3pSolver(int argc, const char* argv[]) {
  std::string desc;
  try {
    cxxopts::Options options(argv[0], " - example command line options");
    options.positional_help("[optional args]").show_positional_help();

    options.set_width(70)
        .set_tab_expansion()
        .allow_unrecognised_options()
        .add_options()("d,desc", "description", cxxopts::value<std::string>());

    auto result = options.parse(argc, argv);

    if (result.count("desc") != 0u) {
      desc = result["desc"].as<std::string>();
    }
  } catch (const cxxopts::OptionException& e) {
    std::cout << "error parsing options: " << e.what() << std::endl;
    exit(1);  // NOLINT
  }

  const uint32_t num_points = 4;  // Number of points

  const Vec3d center(0.5 / sqrt(2.0), -0.5 / sqrt(2.0),
                     0.5 / sqrt(2.0));  // Center of spherical point cloud.

  std::mt19937 engine(123456789);

#if 0
  {
    auto [pts, kps, fx, fy] = CreatePointCloud(
        1000000, 1920, 1920, 60 * M_PI / 180.0, 60 * M_PI / 180.0, 0.01, 100.0,
        Mat3d::Identity(), Vec3d::Zero(), &engine);
    double max_fov_x     = 0.0;
    double max_fov_x_dot = 0.0;
    double max_fov_y     = 0.0;
    double max_fov_y_dot = 0.0;
    for (const auto& pt : pts) {
      if (pt.z() > 0.0) {
        max_fov_x = std::max(max_fov_x,
                             atan(abs(pt.x() / pt.z())) * 2.0 * 180.0 / M_PI);
        max_fov_x_dot =
            std::max(max_fov_x_dot,
                     acos(Vec3d(0.0, 0.0, 1.0)
                              .dot(Vec3d(pt.x(), 0.0, pt.z()).normalized())) *
                         2.0 * 180.0 / M_PI);
        max_fov_y = std::max(max_fov_y,
                             atan(abs(pt.y() / pt.z())) * 2.0 * 180.0 / M_PI);
        max_fov_y_dot =
            std::max(max_fov_y_dot,
                     acos(Vec3d(0.0, 0.0, 1.0)
                              .dot(Vec3d(0.0, pt.y(), pt.z()).normalized())) *
                         2.0 * 180.0 / M_PI);
      }
    }
    std::cout << "max fov x from tan: " << max_fov_x << std::endl;
    std::cout << "max fov x from dot: " << max_fov_x_dot << std::endl;
    std::cout << "max fov y from tan: " << max_fov_y << std::endl;
    std::cout << "max fov y from dot: " << max_fov_y_dot << std::endl;
  }
#endif

  Output output              = {};
  output.name                = "p3p solver";
  output.date                = GetUtcTime();
  output.desc                = desc;
  output.width               = 1920;
  output.height              = 1080;
  output.num_bench_itr       = 100000;
  output.num_found_solutions = 0;

  std::vector<int64_t> elapsed_times;
  for (size_t bench_itr = 0; bench_itr < output.num_bench_itr; ++bench_itr) {
    const auto camera = CreateCameras(&engine);

    const Mat3d rotation_wc    = camera.rotation_cw.transpose();
    const Vec3d translation_wc = -rotation_wc * camera.transltaion_cw;

    auto [points, keypts, fx, fy] = CreatePointCloud(
        num_points, output.width, output.height, 45.0 * M_PI / 180.0,
        90.0 * M_PI / 180.0, 0.5, 10.0, rotation_wc, translation_wc, &engine);

    const eigen_vector<Vec3d> three_points(points.begin(), points.end() - 1);

    const bool add_noise = true;
    const double sigma_x = 1.0 / 3.0;
    const double sigma_y = 1.0 / 3.0;

    const eigen_vector<Vec3d> bearings =
        ToBearingVectors(keypts, fx, fy, add_noise, &engine, sigma_x, sigma_y);
    const eigen_vector<Vec3d> three_bearings(bearings.begin(),
                                             bearings.end() - 1);

    const auto start  = std::chrono::high_resolution_clock::now();
    const auto result = gsrap::SolveP3pProblemWithoutRansac(
        three_bearings, three_points);  //  TODO(kyawakyaw): Check accuracy
    const auto end = std::chrono::high_resolution_clock::now();

    const std::chrono::nanoseconds elapsed_time =
        std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);

    elapsed_times.emplace_back(elapsed_time.count());

    Mat3d rotation;
    Vec3d translation;
    double min_error = std::numeric_limits<double>::max();

    if (result.empty()) {
      continue;
    }
    output.num_found_solutions++;

    const double u = fx * bearings.back().x() / bearings.back().z();
    const double v = fy * bearings.back().y() / bearings.back().z();
    for (const auto& rt : result) {
      const Vec3d bearing_ =
          (rt.first * points.back() + rt.second).normalized();

      const double u_ = fx * bearing_.x() / bearing_.z();
      const double v_ = fy * bearing_.y() / bearing_.z();

      const double error = (u - u_) * (u - u_) / (sigma_x * sigma_x) +
                           (v - v_) * (v - v_) / (sigma_y * sigma_y);
      if (error < min_error) {
        rotation    = rt.first;
        translation = rt.second;
        min_error   = error;
      }
    }

#if 0
    if (flag) {
      std::cout << "R  :\n" << rotation << std::endl;
      std::cout << "t  :\n" << translation << std::endl;
      std::cout << "gtR:\n" << camera.rotation_cg << std::endl;
      std::cout << "gtt:\n" << camera.transltaion_cg << std::endl;
      std::cout << "\n" << std::endl;
    } else {
      std::cout << "error" << std::endl;
    }
#endif

    // Fill output.result
    output.result.rotation.emplace_back(Mat3d2StdVector(rotation));
    output.result.translation.emplace_back(Vec3d2StdVector(translation));
    output.result.gt_rotation.emplace_back(Mat3d2StdVector(camera.rotation_cw));
    output.result.gt_translation.emplace_back(
        Vec3d2StdVector(camera.transltaion_cw));
    output.result.fx.emplace_back(fx);
    output.result.fy.emplace_back(fy);
    output.result.time_ns.emplace_back(double(elapsed_time.count()));

    {
      std::vector<double> depth;
      for (const auto& pt : three_points) {
        const auto ptc = camera.rotation_cw * pt + camera.transltaion_cw;
        depth.emplace_back(ptc.z());
      }
      const double depth_mean =
          std::accumulate(depth.begin(), depth.end(), 0.0) /
          double(depth.size());

      std::vector<double> tmp;
      tmp.reserve(depth.size());
      for (const auto& d : depth) {
        tmp.emplace_back(pow(d - depth_mean, 2.0));
      }
      const double depth_std_deviation = sqrt(
          std::accumulate(tmp.begin(), tmp.end(), 0.0) / double(tmp.size()));

      output.result.depth_mean.emplace_back(depth_mean);
      output.result.depth_std_deviation.emplace_back(depth_std_deviation);

      output.result.number_of_solutions.emplace_back(uint32_t(result.size()));
    }
  }

  SaveOutputAsJson("out-" + output.date + ".json", output);

  std::sort(elapsed_times.begin(), elapsed_times.end());

  const double mean =
      double(std::accumulate(elapsed_times.begin(), elapsed_times.end(), 0)) /
      double(elapsed_times.size());

  const int64_t median = elapsed_times.at(elapsed_times.size() / 2);

  const int64_t maxi = elapsed_times.back();
  std::cout << "mean: " << mean << " ns, median: " << median
            << " ns, max: " << maxi << " ns" << std::endl;
}

}  // namespace gsrap::benchmark

int main(int argc, const char* argv[]) {
  gsrap::benchmark::BenchP3pSolver(argc, argv);
  return 0;
}
