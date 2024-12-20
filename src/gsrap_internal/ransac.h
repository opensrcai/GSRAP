#ifndef GSRAP_INTERNAL_RANSAC_H
#define GSRAP_INTERNAL_RANSAC_H

#include <assert.h>
#include <math.h>

#include <algorithm>
#include <atomic>
#include <limits>
#include <memory>
#include <mutex>
#include <optional>
#include <random>
#include <thread>
#include <utility>
#include <vector>

#include "gsrap/gsrap-def.h"

namespace gsrap {

/// Struct that is packed RANSAC result
template <typename Iterator, typename Model_t>
struct RansacResult {
  RansacResult(std::unique_ptr<Model_t>&& _model,
               std::vector<Iterator>&& _inliers, const size_t _num_inliers,
               const double& _inlier_ratio)
      : model(std::move(_model)),
        inliers(std::move(_inliers)),
        num_inliers(_num_inliers),
        inlier_ratio(_inlier_ratio) {}

  /// Model
  std::unique_ptr<Model_t> model;

  /// Iterators to inliers
  std::vector<Iterator> inliers;

  /// number of inliers
  size_t num_inliers;
  /// inlier rates
  double inlier_ratio;
};

// This function checks RansacPolicy
inline RansacPolicy CheckRansacPolicy(RansacPolicy ransac_policy) {
  assert(ransac_policy.num_ransac_itr_lower_limit <=
         ransac_policy.num_ransac_itr_upper_limit);
  if (ransac_policy.num_ransac_itr_lower_limit >
      ransac_policy.num_ransac_itr_upper_limit) {
    std::swap(ransac_policy.num_ransac_itr_lower_limit,
              ransac_policy.num_ransac_itr_upper_limit);
  }

  ransac_policy.probability =
      std::clamp(ransac_policy.probability, 0.0, 0.999999999);

  assert(ransac_policy.num_threads >= 1);
  ransac_policy.num_threads = std::max(uint32_t(1), ransac_policy.num_threads);

  return ransac_policy;
}

template <typename T>
T Pow(T a, size_t x) {
  T ret = T(1);
  while (x > 0) {
    if (x & 1) {
      ret *= a;
    }
    a *= a;
    x >>= 1;
  }
  return ret;
}

/**
 * @fn
 *
 * @brief This function execute RANSAC
 * @param ransac_policy RANSAC policy.
 * @param first Interator to data.
 * @param last Interator to data.
 * @param FitModels Function which fits model from sampled data.
 * @param IsInlier Function which determins if data are inliers for the model.
 * @param Sample Function which samples Iterators of data.
 * @return Pair of optional of struct RansacResult and struct RansacReport
 * @details TODO
 */
template <typename Iterator, typename FitModelsFunc, typename IsInlierFunc,
          typename SampleFunc>
// TODO(kywakyawa): refectoring
auto Ransac(RansacPolicy ransac_policy, const Iterator first,  // NOLINT
            const Iterator last, const FitModelsFunc FitModels,
            const IsInlierFunc IsInlier, const SampleFunc Sample)
    -> std::pair<
        std::optional<RansacResult<
            Iterator, typename decltype(FitModels(
                          std::declval<const std::vector<Iterator>&>()))::
                          value_type::element_type>>,
        RansacReport> {
  // Derive the type of the model from the return type of fuction
  // FitModels.
  using Model = typename decltype(FitModels(
      std::declval<const std::vector<Iterator>&>()))::value_type::element_type;

  // Prepare a variable to store the return value.
  std::optional<RansacResult<Iterator, Model>> ret;

  // Calculate number of data from iterator.
  const auto num_data = size_t(std::distance(first, last));

  ransac_policy = CheckRansacPolicy(ransac_policy);

  const bool early_stop =
      (ransac_policy.flags & RANSAC_POLICY_FLAGS_EARLY_STOP);

  RansacReport ransac_report = {};

  // Sets global iteration counter
  std::atomic<uint32_t> itr_cnt{0};
  const size_t itr_lower_limit =
      ransac_policy.num_ransac_itr_lower_limit;  //  immutable
  size_t itr_upperlimit  = ransac_policy.num_ransac_itr_upper_limit;  // mutable
  size_t best_num_inliers = 0;     // variable for inlier check early stopping
  std::mutex mtx_ret;             // mutex for variable ret.
  std::mutex mtx_itr_upperlimit;  // mutex for variable itr_upperlimit

  // Since the calculation becomes unstable when the inlier rate is small, we
  // use a threshold value compared with a power of the inlier rate. To compare
  // with this threshold is equivalent to determining whether the number of
  // times calculated by the inlier rate does not exceed the maximum number of
  // times in the policy.
  const double prob_sample_success_thr =
      1 -
      std::pow(1.0 - ransac_policy.probability, 1.0 / double(itr_upperlimit));

  std::vector<std::thread> threads;

  // Define lambda function.
  // TODO(kyawakyawa): refactoring
  const auto Work = [&](const size_t thread_id) {  // NOLINT
    // TODO(kywakyawa): use better seed
    std::mt19937 engine(thread_id * size_t(12345));
    while (true) {
      // Increment the counter to thread-safe
      const uint32_t _itr_cnt = itr_cnt++;

      // If the number of iterations is greater than or equal to variable
      // itr_lower_limit, we check that it reaches the upper limit.
      if (_itr_cnt >= itr_lower_limit) {
        // If the number of iterations reaches the upper limit, we get out of
        // the loop. If early stopping is enabled, we note that variable
        // itr_upperlimit is mutable.
        if (early_stop) {
          std::lock_guard<std::mutex> lock(mtx_itr_upperlimit);
          if (_itr_cnt >= itr_upperlimit) {
            break;
          }
        } else {
          if (_itr_cnt >= itr_upperlimit) {
            break;
          }
        }
      }

      // Sample data to fit models.
      const std::vector<Iterator> sampled_data = Sample(first, last, &engine);

      // Fit models from sampled data.
      std::vector<std::unique_ptr<Model>> models = FitModels(sampled_data);

      // For each model, we compute inlier rate and update result.
      for (auto& model : models) {
        std::vector<Iterator> inliers;
        inliers.reserve(models.size());

        bool inlier_check_early_stop =
            false;  // variable for inlier check early stopping
        size_t num_data_no_checked =
            num_data;  // variable for inlier check early stopping
        // Correct iterators to inlier data.
        for (Iterator itr = first; itr != last; ++itr) {
          // Inlier check early stopping
          if (best_num_inliers > inliers.size() + num_data_no_checked) {
            inlier_check_early_stop = true;
            break;
          }
          if (IsInlier(*itr, *(model.get()))) {
            inliers.emplace_back(itr);
          }
          num_data_no_checked--;
        }

        // If we get a better inlier rate, we update the results.
        if (!inlier_check_early_stop) {
          const size_t num_inliers = inliers.size();
          std::lock_guard<std::mutex> lock(mtx_ret);
          if (!ret.has_value() || num_inliers > ret.value().num_inliers) {
            const double inlier_ratio = double(num_inliers) / double(num_data);
            ret.emplace(std::move(model), std::move(inliers), num_inliers,
                        inlier_ratio);
            best_num_inliers = num_inliers;  // updates best_num_inliers for
                                            // inlier check early stopping

            // If early stoppping is enabled, we update variable itr_upperlimit
            // thread-safely based on probability.
            if (early_stop) {
              double prob_sample_success = 1;
              if (ransac_policy.flags &
                  RANSAC_POLICY_FLAGS_USE_PROBABILITY_WITHOUT_DUPLICATION_SAMPLE) {
                // Probability that all will be inliers when sampled without
                // allowing for duplication.
                for (int64_t i = 0; size_t(i) < sampled_data.size(); ++i) {
                  prob_sample_success *=
                      std::max<double>(0, double(int64_t(num_inliers) - i)) /
                      double(int64_t(num_data) - i);
                }
              } else {
                // This is a commonly used probability.
                prob_sample_success =
                    Pow(inlier_ratio, sampled_data.size());
              }

              if (prob_sample_success > prob_sample_success_thr) {
                std::lock_guard<std::mutex> lock_itr_upperlimit(
                    mtx_itr_upperlimit);
                const size_t N =
                    size_t(log(1.0 - ransac_policy.probability) /
                           log(std::max(std::numeric_limits<double>::epsilon(),
                                        1.0 - prob_sample_success)));
                itr_upperlimit = std::min(itr_upperlimit, N);
              }
            }
          }
        }
      }
    }
  };

  // Lauche threads
  for (size_t thread_id = 0; thread_id < ransac_policy.num_threads - 1;
       ++thread_id) {
    threads.emplace_back(std::bind(Work, thread_id));
  }
  // Main thread
  Work(ransac_policy.num_threads - 1);

  for (auto& thr : threads) {
    thr.join();
  }

  // Update ransac_report
  ransac_report.num_iteration = itr_cnt - ransac_policy.num_threads;
  ransac_report.ransac_termination_info =
      (itr_upperlimit < ransac_policy.num_ransac_itr_upper_limit
           ? RansacTerminationInfo::EARLY_STOPED
           : RansacTerminationInfo::REACHED_UPPER_LIMIT);
  return std::make_pair(std::move(ret), ransac_report);
}

}  // namespace gsrap

#endif  // GSRAP_INTERNAL_RANSAC_H
