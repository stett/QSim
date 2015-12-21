
# GTest
add_subdirectory(googletest)

# Include everything
include_directories(${QSim_SOURCE_DIR}/include
                    ${QSim_SOURCE_DIR}/test)
include_directories(${gtest_SOURCE_DIR}/include)

# Define the executable test program
add_executable(QSimTest main.cpp)

# Link everything to the FrameTest executable
target_link_libraries(QSimTest gtest gtest_main)
