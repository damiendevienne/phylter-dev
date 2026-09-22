cmake_minimum_required(VERSION 3.21)

project(phylter2 VERSION 0.1.0 LANGUAGES CXX)

option(PHYLTER2_BUILD_TESTS "Build the C++ test suite" ON)
option(PHYLTER2_OPENMP "Enable controlled parallel RV computation" ON)
option(PHYLTER2_BUNDLE_RUNTIME "Bundle Linux numerical and compiler runtime libraries" OFF)

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

add_library(phylter2_core
  src/matrix.cpp
  src/lapack.cpp
  src/distatis.cpp
  src/newick.cpp
  src/preparation.cpp
  src/statistics.cpp
  src/medcouple.cpp
  src/analysis.cpp
  src/rv.cpp
  src/iterative_eigen.cpp
)

target_include_directories(phylter2_core
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)
target_compile_features(phylter2_core PUBLIC cxx_std_20)
target_link_libraries(phylter2_core PUBLIC LAPACK::LAPACK BLAS::BLAS)
set(PHYLTER2_HAVE_OPENMP OFF)
if(PHYLTER2_OPENMP)
  find_package(OpenMP)
  if(OpenMP_CXX_FOUND)
    set(PHYLTER2_HAVE_OPENMP ON)
    target_link_libraries(phylter2_core PUBLIC OpenMP::OpenMP_CXX)
  endif()
endif()

if(MSVC)
  target_compile_options(phylter2_core PRIVATE /W4 /permissive-)
else()
  target_compile_options(phylter2_core PRIVATE -Wall -Wextra -Wpedantic -Wconversion)
endif()

add_executable(phylter src/main.cpp)
target_link_libraries(phylter PRIVATE phylter2_core)
target_compile_features(phylter PRIVATE cxx_std_20)

include(GNUInstallDirs)
install(TARGETS phylter
  RUNTIME_DEPENDENCY_SET phylterRuntime
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT Runtime
)
install(TARGETS phylter2_core
  EXPORT phylter2Targets
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Development
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Development
)
install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} COMPONENT Development)
install(FILES LICENSE NOTICE README.md VALIDATION.md RV.md BENCHMARK-MATRIX-FREE.md
  DESTINATION ${CMAKE_INSTALL_DATADIR}/phylter COMPONENT Runtime)
install(FILES packaging/README.md DESTINATION ${CMAKE_INSTALL_DATADIR}/phylter/packaging COMPONENT Runtime)
install(FILES tests/fixtures/carnivora.nwk tests/fixtures/carnivora.outliers.tsv
  tests/fixtures/carnivora.scores.tsv tests/fixtures/PROVENANCE.txt
  DESTINATION ${CMAKE_INSTALL_DATADIR}/phylter/examples COMPONENT Runtime)
if(PHYLTER2_BUNDLE_RUNTIME)
  if(NOT CMAKE_SYSTEM_NAME STREQUAL "Linux")
    message(FATAL_ERROR "Runtime bundling is currently validated only on Linux")
  endif()
  set_target_properties(phylter PROPERTIES INSTALL_RPATH "$ORIGIN/../${CMAKE_INSTALL_LIBDIR}")
  # RPATH is inherited by LAPACK's transitive dependencies (RUNPATH is not).
  target_link_options(phylter PRIVATE "-Wl,--disable-new-dtags")
  install(RUNTIME_DEPENDENCY_SET phylterRuntime
    PRE_EXCLUDE_REGEXES "^ld-linux" "^libc\\.so" "^libm\\.so" "^libmvec\\.so" "^libdl\\.so"
                       "^libpthread\\.so" "^librt\\.so"
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Runtime)
endif()

set(CPACK_PACKAGE_NAME phylter2)
set(CPACK_PACKAGE_VERSION ${PROJECT_VERSION})
set(CPACK_PACKAGE_FILE_NAME "phylter2-${PROJECT_VERSION}-${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR}")
set(CPACK_GENERATOR TGZ)
set(CPACK_COMPONENTS_ALL Runtime)
set(CPACK_ARCHIVE_COMPONENT_INSTALL ON)
set(CPACK_COMPONENTS_GROUPING ALL_COMPONENTS_IN_ONE)
set(CPACK_INCLUDE_TOPLEVEL_DIRECTORY ON)
include(CPack)

if(PHYLTER2_BUILD_TESTS)
  include(CTest)
  if(BUILD_TESTING)
    add_executable(phylter2-core-tests tests/core_tests.cpp)
    target_link_libraries(phylter2-core-tests PRIVATE phylter2_core)
    add_test(NAME phylter2-core-tests COMMAND phylter2-core-tests)
    set_tests_properties(phylter2-core-tests PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1")
    add_executable(phylter2-rv-tests tests/rv_tests.cpp)
    target_link_libraries(phylter2-rv-tests PRIVATE phylter2_core)
    add_test(NAME phylter2-rv-tests COMMAND phylter2-rv-tests)
    set_tests_properties(phylter2-rv-tests PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1")
    add_executable(phylter2-regression-tests tests/regression_tests.cpp)
    target_link_libraries(phylter2-regression-tests PRIVATE phylter2_core)
    add_test(NAME phylter2-r-reference COMMAND phylter2-regression-tests
      ${CMAKE_CURRENT_SOURCE_DIR}/tests/fixtures)
    add_test(NAME phylter2-cli-check
      COMMAND phylter check --trees ${CMAKE_CURRENT_SOURCE_DIR}/tests/fixtures/simple-tree.nwk)
    add_test(NAME phylter2-cli-analysis COMMAND ${CMAKE_COMMAND}
      -DPHYLTER=$<TARGET_FILE:phylter>
      -DFIXTURES=${CMAKE_CURRENT_SOURCE_DIR}/tests/fixtures
      -DOUTPUT_DIR=${CMAKE_CURRENT_BINARY_DIR}/cli-tests
      -DPARALLEL=${PHYLTER2_HAVE_OPENMP}
      -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/cli-tests.cmake)
  endif()
endif()
