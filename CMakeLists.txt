cmake_minimum_required(VERSION 3.16)
project(pluto_hi_monitor LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(PkgConfig REQUIRED)

# Dépendances IIO et FFTW
find_library(IIO_LIB iio REQUIRED)
find_library(AD9361_LIB ad9361 REQUIRED)
find_library(FFTW3_LIB fftw3 REQUIRED)
find_library(M_LIB m REQUIRED)

# Dépendance sigutils via pkg-config
pkg_check_modules(SIGUTILS REQUIRED sigutils)

add_executable(spectre_integration src/main.cpp)

target_include_directories(spectre_integration PRIVATE
    ${SIGUTILS_INCLUDE_DIRS}
)

target_link_libraries(spectre_integration PRIVATE
    ${IIO_LIB}
    ${AD9361_LIB}
    ${FFTW3_LIB}
    ${M_LIB}
    ${SIGUTILS_LIBRARIES}
)

target_compile_options(spectre_integration PRIVATE -O3 -march=native)
