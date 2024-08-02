#include <cmath>
#include <cstddef>
#include <vector>
#include <benchmark/benchmark.h>

#include <libff_liop/algebra/fields/binary/gf64.hpp>
#include <libff_liop/common/utils.hpp>
#include "libiop/algebra/utils.hpp"

namespace libiop {

static void BM_all_gf64_subset_sums(benchmark::State &state)
{
    const size_t sz = state.range(0);
    const size_t log_sz = libff_liop::log2(sz);

    const std::vector<libff_liop::gf64> basis = random_vector<libff_liop::gf64>(log_sz);
    const libff_liop::gf64 shift = libff_liop::gf64::random_element();

    for (auto _ : state)
    {
        all_subset_sums<libff_liop::gf64>(basis, shift);
    }

    state.SetItemsProcessed(state.iterations() * sz);
}

BENCHMARK(BM_all_gf64_subset_sums)->Range(1ull<<4, 1ull<<20)->Unit(benchmark::kMicrosecond);

static void BM_random_gf64_vector(benchmark::State &state)
{
    const size_t sz = state.range(0);

    for (auto _ : state)
    {
        const std::vector<libff_liop::gf64> vec = random_vector<libff_liop::gf64>(sz);
    }

    state.SetItemsProcessed(state.iterations() * sz);
}

BENCHMARK(BM_random_gf64_vector)->Range(1ull<<4, 1ull<<20)->Unit(benchmark::kMicrosecond);

}

BENCHMARK_MAIN();
