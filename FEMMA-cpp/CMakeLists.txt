cmake_minimum_required(VERSION 3.16)
project(femma_cpp VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

if(MSVC)
  # MSVC: Release defaults to /O2; /DNOMINMAX avoids collisions with Eigen
  add_compile_options(/W3 /DNOMINMAX)
else()
  set(CMAKE_CXX_FLAGS_RELEASE "-O2")
  add_compile_options(-Wall -Wextra)
endif()

option(FEMMA_BUILD_TESTS    "Build the test suite"            ON)
option(FEMMA_BUILD_EXAMPLES "Build the example programmes"    ON)
option(FEMMA_BUILD_TOOLS    "Build the MATLAB cross-validation harnesses" OFF)

# Eigen: all three paths create an Eigen3::Eigen CMake target so libraries
# can always link it without caring about the install location.
find_package(Eigen3 3.3 QUIET NO_MODULE)
if(TARGET Eigen3::Eigen)
  message(STATUS "FEMMA: using Eigen3::Eigen target (system or vcpkg)")
elseif(EXISTS ${CMAKE_SOURCE_DIR}/third_party/Eigen)
  add_library(Eigen3_vendored INTERFACE)
  target_include_directories(Eigen3_vendored SYSTEM INTERFACE
      ${CMAKE_SOURCE_DIR}/third_party)
  add_library(Eigen3::Eigen ALIAS Eigen3_vendored)
  message(STATUS "FEMMA: using vendored Eigen under third_party/")
else()
  include(FetchContent)
  FetchContent_Declare(eigen
    GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
    GIT_TAG 3.4.0
    GIT_SHALLOW TRUE)
  FetchContent_MakeAvailable(eigen)
  message(STATUS "FEMMA: fetched Eigen 3.4.0")
endif()

# --- optimiser library (GPL, translated from Svanberg's MMA/GCMMA) ----------
add_library(femma_mma
  src/mma/subsolv.cpp
  src/mma/mmasub.cpp
  src/mma/kktcheck.cpp
  src/mma/asymp.cpp
  src/mma/gcmmasub.cpp
  src/mma/concheck.cpp
  src/mma/raaupdate.cpp)
target_include_directories(femma_mma PUBLIC include)
target_link_libraries(femma_mma PUBLIC Eigen3::Eigen)

# --- core FEMMA library (original translation of the FEMMA MATLAB) ----------
add_library(femma_core
  src/core/operator_create.cpp
  src/core/state_create.cpp
  src/core/mesh_fem.cpp
  src/core/ele_matrix_linear.cpp
  src/core/ele_matrix_grad_linear.cpp
  src/core/ele_matrix_quad.cpp
  src/core/ele_matrix_grad_quad.cpp
  src/core/ele_matrix_hermite.cpp
  src/core/ele_matrix_grad_hermite.cpp
  src/core/global_stiff.cpp
  src/core/femsolve.cpp
  src/core/helm_matrix.cpp
  src/core/helm_solve.cpp
  src/core/ensemble_spin.cpp
  src/core/combine_grad.cpp
  src/core/sigmoid.cpp
  src/core/power_constr.cpp
  src/core/spin_solve_phase.cpp
  src/core/spin_solve_xy.cpp
  src/core/bloch.cpp
  src/core/bloch_solve_phase.cpp
  src/core/mma_init_phase.cpp
  src/core/mma_init_xy.cpp
  src/core/fem_mma.cpp)
target_include_directories(femma_core PUBLIC include)
target_link_libraries(femma_core PUBLIC Eigen3::Eigen)

# the driver calls the optimiser, so core depends on mma; this is the only
# edge across the GPL boundary
target_link_libraries(femma_core PUBLIC femma_mma)

# OpenMP for the embarrassingly parallel ensemble loop
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
  target_link_libraries(femma_core PUBLIC OpenMP::OpenMP_CXX)
endif()

# convenience alias so downstream targets can link a single name
add_library(femma_fem INTERFACE)
target_link_libraries(femma_fem INTERFACE femma_core)

if(FEMMA_BUILD_TESTS)
  enable_testing()
  macro(femma_test name)
    add_executable(test_${name} tests/test_${name}.cpp)
    target_link_libraries(test_${name} PRIVATE femma_fem)
    add_test(NAME ${name} COMMAND test_${name})
  endmacro()
  femma_test(operators)
  femma_test(assembly)
  femma_test(solve)
  femma_test(gradient)
  femma_test(helm)
  femma_test(ensemble)
  femma_test(spinsolve)
  femma_test(mma)
  femma_test(optimise)
  femma_test(spinsolve_xy)
  femma_test(optimise_xy)
  femma_test(gcmma)
  femma_test(quad)
  femma_test(optimise_quad)
  femma_test(hermite)
  femma_test(optimise_hermite)
  femma_test(bloch)
  femma_test(optimise_bloch)
  femma_test(bloch_xy_grad)
endif()

if(FEMMA_BUILD_EXAMPLES)
  macro(femma_example name)
    add_executable(${name} examples/${name}.cpp)
    target_link_libraries(${name} PRIVATE femma_fem)
  endmacro()
  femma_example(ex_femma2)
  femma_example(lvn_phase)
  femma_example(element_orders)
  femma_example(bloch_excitation)
endif()

# cross-validation harnesses that diff against MATLAB dumps in /tmp; off by
# default because they need the reference data files at run time
if(FEMMA_BUILD_TOOLS)
  foreach(tool dump_case dump_spinsolve bv_iter1 bv_traj bv_wave bl_xval)
    add_executable(${tool} tools/${tool}.cpp)
    target_link_libraries(${tool} PRIVATE femma_fem)
  endforeach()
endif()
