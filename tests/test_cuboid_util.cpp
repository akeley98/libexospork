#include "test_cuboid_util.hpp"

#include <stdio.h>
#include <stdlib.h>
#include <utility>
#include <vector>

#include "../lib/util/cuboid_util.hpp"

namespace camspork {

struct TestCuboidInterval
{
    uint32_t lo;
    uint32_t hi;
    bool operator==(TestCuboidInterval other) const
    {
        return lo == other.lo && hi == other.hi;
    }
};

class TestCuboidCallback
{
    int line;
    std::vector<TestCuboidInterval> actual;
    std::vector<TestCuboidInterval> expected;

  public:
    TestCuboidCallback(int line, std::vector<TestCuboidInterval> expected)
      : line(line)
      , actual({})
      , expected(std::move(expected))
    {
    }

    void check()
    {
        if (actual == expected) {
            printf("Test passed, line %i\n", line);
            return;
        }
        printf("Test failed, line %i\n", line);
        printf("  actual:");
        for (TestCuboidInterval t: actual) {
            printf(" (%u, %u)", t.lo, t.hi);
        }
        printf("\n  expected:");
        for (TestCuboidInterval t: expected) {
            printf(" (%u, %u)", t.lo, t.hi);
        }
        printf("\n");
    }

    void operator() (uint32_t lo, uint32_t hi)
    {
        actual.push_back({lo, hi});
    }
};

static void cuboid_util_test_case(
    int line,
    std::vector<uint32_t> outer,
    std::vector<uint32_t> offset,
    std::vector<uint32_t> inner,
    std::vector<TestCuboidInterval> expected)
{
    TestCuboidCallback callback(line, std::move(expected));
    cuboid_to_intervals<uint32_t>(
        outer.begin(), outer.end(), offset.begin(), offset.end(), inner.begin(), inner.end(), callback);
    callback.check();
}

void test_cuboid_util()
{
    cuboid_util_test_case(
        __LINE__,
        {20, 10, 40}, {0, 0, 0}, {20, 10, 40}, {{0, 8000}}
    );
    cuboid_util_test_case(
        __LINE__,
        {20, 10, 40}, {5, 0, 0}, {15, 10, 40}, {{2000, 8000}}
    );
    cuboid_util_test_case(
        __LINE__,
        {20, 10, 40}, {5, 0, 0}, {2, 2, 30}, {{2000, 2030}, {2040, 2070}, {2400, 2430}, {2440, 2470}}
    );
    cuboid_util_test_case(
        __LINE__,
        {20, 10, 40}, {5, 0, 0}, {2, 2, 0}, {}
    );
    cuboid_util_test_case(
        __LINE__,
        {20, 10, 40}, {5, 0, 0}, {0, 2, 40}, {}
    );
    cuboid_util_test_case(
        __LINE__,
        {10, 10, 100}, {5, 4, 3}, {3, 2, 2}, {
            {5403, 5405}, {5503, 5505},
            {6403, 6405}, {6503, 6505},
            {7403, 7405}, {7503, 7505},
        }
    );
    cuboid_util_test_case(
        __LINE__,
        {10, 10, 100}, {5, 4, 0}, {3, 2, 100}, {
            {5400, 5600},
            {6400, 6600},
            {7400, 7600},
        }
    );
    cuboid_util_test_case(
        __LINE__,
        {10, 10, 100}, {5, 4, 0}, {3, 1, 50}, {
            {5400, 5450},
            {6400, 6450},
            {7400, 7450},
        }
    );
    cuboid_util_test_case(
        __LINE__,
        {}, {}, {}, {{0, 1}}
    );
    cuboid_util_test_case(
        __LINE__,
        {200}, {5}, {15}, {{5, 20}}
    );
    cuboid_util_test_case(
        __LINE__,
        {5, 100}, {0, 18}, {5, 10}, {{18, 28}, {118, 128}, {218, 228}, {318, 328}, {418, 428}}
    );
}

}
