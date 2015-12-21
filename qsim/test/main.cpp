#include <stdio.h>
#include "gtest/gtest.h"
#include "test_QSimMath.h"

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    int err = RUN_ALL_TESTS();
    char str = ' ';
    std::cin >> str;
    return err;
}