#ifndef GSRAP_INTERNAL_SAMPLE_H
#define GSRAP_INTERNAL_SAMPLE_H

#include <assert.h>
#include <stddef.h>
#include <stdint.h>

#include <iterator>
#include <random>
#include <unordered_set>
#include <vector>

namespace gsrap {

// Function object to sample `N` pairs of bearing vector.
template <typename Iterator, size_t N>
struct SampleN {
  std::vector<Iterator> operator()(const size_t num_data, Iterator first,
                                   [[maybe_unused]] Iterator last,
                                   std::mt19937* engine) const {
    assert(num_data == size_t(std::distance(first, last)));
    if (num_data < N) {
      // TODO(kyawakyawa): error handling
      assert(false);
      return std::vector<Iterator>();
    }

    std::unordered_set<int64_t> sampled_ids;
    std::uniform_int_distribution<int64_t> dist(0, int64_t(num_data) - 1);
    while (sampled_ids.size() < N) {
      sampled_ids.emplace(dist(*engine));
    }

    std::vector<Iterator> ret;
    ret.reserve(sampled_ids.size());
    for (const int64_t sampled_id : sampled_ids) {
      // NOTE: If Iterator is std::unordered_map<>::const_iterator, this process
      // is slow.
      ret.emplace_back(std::next(first, sampled_id));
    }

    return ret;
  }
};

// Function object to sample `num_sample` pairs of bearing vector.
template <typename Iterator>
struct Sample {
  std::vector<Iterator> operator()(const size_t num_data,
                                   const size_t num_sample, Iterator first,
                                   [[maybe_unused]] Iterator last,
                                   std::mt19937* engine) const {
    assert(num_data == size_t(std::distance(first, last)));
    if (num_data < num_sample) {
      // TODO(kyawakyawa): error handling
      assert(false);
      return std::vector<Iterator>();
    }

    std::unordered_set<int64_t> sampled_ids;
    std::uniform_int_distribution<int64_t> dist(0, int64_t(num_data) - 1);
    while (sampled_ids.size() < num_sample) {
      sampled_ids.emplace(dist(*engine));
    }

    std::vector<Iterator> ret;
    ret.reserve(sampled_ids.size());
    for (const int64_t sampled_id : sampled_ids) {
      // NOTE: If Iterator is std::unordered_map<uint32_t,
      // uint32_t>::const_iterator, this process is slow.
      ret.emplace_back(std::next(first, sampled_id));
    }
    return ret;
  }
};

}  // namespace gsrap
#endif  // GSRAP_INTERNAL_SAMPLE_H
