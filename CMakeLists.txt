cmake_minimum_required(VERSION 3.27)
project(multifloat LANGUAGES CXX Fortran)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
enable_testing()

include_directories(src)
find_program(FYPP_EXECUTABLE fypp REQUIRED)

set(MULTIFLOATS_GENERATED_DIR ${CMAKE_CURRENT_BINARY_DIR}/generated)
set(MULTIFLOATS_GENERATED_SRC ${MULTIFLOATS_GENERATED_DIR}/multifloats.f90)

add_custom_command(
    OUTPUT ${MULTIFLOATS_GENERATED_SRC}
    COMMAND ${CMAKE_COMMAND} -E make_directory ${MULTIFLOATS_GENERATED_DIR}
    COMMAND ${FYPP_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/fsrc/multifloats.fypp ${MULTIFLOATS_GENERATED_SRC}
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/fsrc/multifloats.fypp
    VERBATIM)

add_library(multifloat-fortran ${MULTIFLOATS_GENERATED_SRC})
set_target_properties(multifloat-fortran PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)
target_include_directories(multifloat-fortran PUBLIC
    ${CMAKE_BINARY_DIR}/modules)

add_executable(test-fortran ./fsrc/test.f90)
target_link_libraries(test-fortran PRIVATE multifloat-fortran)
add_test(NAME test-fortran COMMAND test-fortran)

add_library(blas-multifloat
    ./blas/wgemm.f90
    ./blas/wtrsm.f90)
target_link_libraries(blas-multifloat PRIVATE multifloat-fortran)

add_executable(multifloat_fuzz ./fsrc/multifloat_fuzz.f90)
target_link_libraries(multifloat_fuzz PRIVATE multifloat-fortran)
add_executable(multifloat_test ./fsrc/multifloat_test.f90)
target_link_libraries(multifloat_test PRIVATE multifloat-fortran)

# C++ tests for src/multifloats.hh.
#
# These tests need __float128 + libquadmath as a high-precision reference,
# which Apple Clang does not provide. We locate a compatible g++ (e.g.
# Homebrew g++-15) and compile the test executables with it directly via
# add_custom_command. Override with -DMULTIFLOATS_QUADMATH_CXX=/path/to/g++.
find_program(MULTIFLOATS_QUADMATH_CXX
    NAMES g++-15 g++-14 g++-13 g++
    PATHS /opt/homebrew/bin /usr/local/bin /usr/bin
    DOC "C++ compiler used to build the C++ tests against __float128")

if(MULTIFLOATS_QUADMATH_CXX)
    set(MF_TEST_INCLUDES -I${CMAKE_CURRENT_SOURCE_DIR}/src)
    set(MF_TEST_FLAGS -std=c++17 -O2 -Wall -Wextra -fext-numeric-literals)
    set(MF_TEST_LIBS -lquadmath)

    foreach(name multifloats_test multifloats_fuzz)
        set(_src ${CMAKE_CURRENT_SOURCE_DIR}/test/${name}.cc)
        set(_bin ${CMAKE_CURRENT_BINARY_DIR}/${name})
        add_custom_command(
            OUTPUT ${_bin}
            COMMAND ${MULTIFLOATS_QUADMATH_CXX}
                    ${MF_TEST_FLAGS}
                    ${MF_TEST_INCLUDES}
                    ${_src}
                    -o ${_bin}
                    ${MF_TEST_LIBS}
            DEPENDS ${_src} ${CMAKE_CURRENT_SOURCE_DIR}/src/multifloats.hh
            COMMENT "Building ${name} with ${MULTIFLOATS_QUADMATH_CXX}"
            VERBATIM)
        add_custom_target(build_${name} ALL DEPENDS ${_bin})
        add_test(NAME ${name} COMMAND ${_bin})
    endforeach()
else()
    message(WARNING "No quadmath-capable g++ found; C++ tests will not be built.")
endif()
