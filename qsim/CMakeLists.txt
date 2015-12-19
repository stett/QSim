
# Set the version of CMake
cmake_minimum_required (VERSION 3.0.1)

# Set project title
project(QSim)

# Set compiler options
if (UNIX)
    add_definitions(-std=c++0x -pthread -Wall -Wextra)
elseif(WIN32)
    add_definitions(/Wall)
endif()

# Include external and internal headers
include_directories(
    ${QSim_SOURCE_DIR}/extern/include
    ${QSim_SOURCE_DIR}/include
    ${QSim_SOURCE_DIR}/source
)

# Look for static libs in the correct directory
link_directories(${QSim_SOURCE_DIR}/extern/lib)

# Mark all source files in the main QSim directory for compilation
file(GLOB_RECURSE QSimSourceFiles "source/*.cpp")

# Create an executable from the source files
add_executable(QSim ${QSimSourceFiles})

# Add some preprocessor definitions, necessary for statically linking SFML
target_compile_definitions(QSim PRIVATE SFML_STATIC)
target_compile_definitions(QSim PRIVATE SFGUI_STATIC)
target_compile_definitions(QSim PRIVATE GSL_DLL)

# Statically link SFML to QSim
if(CMAKE_BUILD_TYPE MATCHES Release)
    target_link_libraries(QSim sfgui-s sfml-graphics-s sfml-window-s sfml-system-s opengl32 freetype jpeg winmm gdi32)
elseif(CMAKE_BUILD_TYPE MATCHES Debug)
    target_link_libraries(QSim sfgui-s-d sfml-graphics-s-d sfml-window-s-d sfml-system-s-d opengl32 freetype jpeg winmm gdi32)
endif()

# Copy dll's and other necessary binaries to the build directory
add_custom_command(
    TARGET QSim
    PRE_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
    "${QSim_SOURCE_DIR}/extern/bin"
    $<TARGET_FILE_DIR:QSim>)