cmake_minimum_required(VERSION 3.0.0)
project(learnLqr VERSION 0.1.0)

include(CTest)
enable_testing()

find_package(Eigen3 REQUIRED)
find_package(OsqpEigen REQUIRED)
find_package(OpenCV REQUIRED)

include_directories(${OPENCV_INCLUDE_DIRS})
include_directories(${EIGEN3_INCLUDE_DIRS})

add_executable(learnLqr main.cpp)

target_link_libraries( learnLqr OsqpEigen::OsqpEigen ${OpenCV_LIBS} ${EIGEN2_LIBS})

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)