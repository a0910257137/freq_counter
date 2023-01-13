project(CCounter)
cmake_minimum_required(VERSION 3.16)
set (CMAKE_CXX_STANDARD 11)
if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
    message(FATAL_ERROR "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there. You may need to remove CMakeCache.txt. ")
endif()

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)
file(GLOB SRC_LIST "${PROJECT_SOURCE_DIR}/src/*.c*")
find_package(OpenCV REQUIRED)
find_package(Threads REQUIRED)

find_path(PYTHON3_NUMPY_INCLUDE_DIR 
        /home2/anders/miniconda3/envs/cpp/lib/python3.5/site-packages/numpy/core/include)
find_package(Python3 COMPONENTS Interpreter Development NumPy REQUIRED)

find_package(PythonLibs 3.0 REQUIRED)
include_directories(${PYTHON3_INCLUDE_DIRS} ${NumPy_INCLUDE_DIRS})
include_directories( ${OpenCV_INCLUDE_DIRS})

# populate matplotlib repository
include(FetchContent)
FetchContent_Declare(
    matplotlib
    GIT_REPOSITORY https://github.com/lava/matplotlib-cpp.git
    GIT_TAG        ef0383f1315d32e0156335e10b82e90b334f6d9f
)
FetchContent_GetProperties(matplotlib)
if(NOT matplotlib_POPULATED)
    FetchContent_Populate(matplotlib)
endif()
include_directories(SYSTEM ${matplotlib_SOURCE_DIR})

add_executable(infer.exe  ${SRC_LIST})
link_directories( ${OpenCV_LIBRARY_DIRS})
target_link_libraries(infer.exe
        PRIVATE
                ${PYTHON_LIBRARIES}
                Python3::NumPy 
                ${OpenCV_LIBS})
