cmake_minimum_required(VERSION 3.21)

project(phylter VERSION 0.1.0 LANGUAGES CXX)

option(PHYLTER2_BUILD_TESTS "Build the C++ test suite" ON)
option(PHYLTER2_OPENMP "Enable controlled parallel RV computation" ON)
option(PHYLTER2_BUNDLE_RUNTIME "Bundle numerical and compiler runtime libraries" OFF)

find_package(Threads REQUIRED)

set(PHYLTER2_BLAS_LIBRARY "" CACHE FILEPATH "Explicit BLAS library (for toolchains where FindBLAS fails)")
set(PHYLTER2_LAPACK_LIBRARY "" CACHE FILEPATH "Explicit LAPACK library (for toolchains where FindLAPACK fails)")
set(PHYLTER2_NUMERICAL_LIBRARIES)
if(PHYLTER2_BLAS_LIBRARY OR PHYLTER2_LAPACK_LIBRARY)
  if(NOT EXISTS "${PHYLTER2_BLAS_LIBRARY}" OR NOT EXISTS "${PHYLTER2_LAPACK_LIBRARY}")
    message(FATAL_ERROR "Both PHYLTER2_BLAS_LIBRARY and PHYLTER2_LAPACK_LIBRARY must name existing files")
  endif()
  list(APPEND PHYLTER2_NUMERICAL_LIBRARIES
    "${PHYLTER2_LAPACK_LIBRARY}" "${PHYLTER2_BLAS_LIBRARY}")
elseif(WIN32 AND DEFINED ENV{CONDA_PREFIX})
  # Conda-forge puts its MSVC import libraries below Library/lib. CMake's
  # FindBLAS does not reliably search that layout on Windows, even when the
  # environment prefix is present in CMAKE_PREFIX_PATH.
  find_library(PHYLTER2_CONDA_BLAS NAMES blas
    PATHS "$ENV{CONDA_PREFIX}/Library/lib" NO_DEFAULT_PATH)
  find_library(PHYLTER2_CONDA_LAPACK NAMES lapack
    PATHS "$ENV{CONDA_PREFIX}/Library/lib" NO_DEFAULT_PATH)
  if(PHYLTER2_CONDA_BLAS AND PHYLTER2_CONDA_LAPACK)
    message(STATUS "Using BLAS/LAPACK from Conda: $ENV{CONDA_PREFIX}")
    list(APPEND PHYLTER2_NUMERICAL_LIBRARIES
      "${PHYLTER2_CONDA_LAPACK}" "${PHYLTER2_CONDA_BLAS}")
  endif()
endif()
if(NOT PHYLTER2_NUMERICAL_LIBRARIES)
  find_package(BLAS REQUIRED)
  find_package(LAPACK REQUIRED)
  list(APPEND PHYLTER2_NUMERICAL_LIBRARIES LAPACK::LAPACK BLAS::BLAS)
endif()

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
target_link_libraries(phylter2_core PUBLIC ${PHYLTER2_NUMERICAL_LIBRARIES})
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
target_link_libraries(phylter PRIVATE phylter2_core Threads::Threads)
if(WIN32)
  target_link_libraries(phylter PRIVATE psapi)
endif()
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
  if(WIN32)
    install(RUNTIME_DEPENDENCY_SET phylterRuntime
      DIRECTORIES "$ENV{CONDA_PREFIX}/Library/bin"
      PRE_EXCLUDE_REGEXES "^api-ms-" "^ext-ms-"
      POST_EXCLUDE_REGEXES ".*[Ww]indows[/\\\\][Ss]ystem32.*"
      RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT Runtime
      LIBRARY DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT Runtime)
  elseif(APPLE)
    set_target_properties(phylter PROPERTIES INSTALL_RPATH "@loader_path/../${CMAKE_INSTALL_LIBDIR}")
    install(RUNTIME_DEPENDENCY_SET phylterRuntime
      DIRECTORIES "$ENV{CONDA_PREFIX}/lib"
      POST_EXCLUDE_REGEXES "^/System/Library/" "^/usr/lib/"
      LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Runtime)
  else()
    set_target_properties(phylter PROPERTIES INSTALL_RPATH "$ORIGIN/../${CMAKE_INSTALL_LIBDIR}")
    # RPATH is inherited by LAPACK's transitive dependencies (RUNPATH is not).
    target_link_options(phylter PRIVATE "-Wl,--disable-new-dtags")
    install(RUNTIME_DEPENDENCY_SET phylterRuntime
      DIRECTORIES "$ENV{CONDA_PREFIX}/lib"
      PRE_EXCLUDE_REGEXES "^ld-linux" "^libc\\.so" "^libm\\.so" "^libmvec\\.so" "^libdl\\.so"
                         "^libpthread\\.so" "^librt\\.so"
      LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Runtime)
  endif()
endif()

set(CPACK_PACKAGE_NAME phylter)
set(CPACK_PACKAGE_VERSION ${PROJECT_VERSION})
set(CPACK_PACKAGE_FILE_NAME "phylter-${PROJECT_VERSION}-${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR}")
if(WIN32)
  set(CPACK_GENERATOR ZIP)
else()
  set(CPACK_GENERATOR TGZ)
endif()
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
      COMMAND phylter --check --trees ${CMAKE_CURRENT_SOURCE_DIR}/tests/fixtures/simple-tree.nwk)
    add_test(NAME phylter2-cli-analysis COMMAND ${CMAKE_COMMAND}
      -DPHYLTER=$<TARGET_FILE:phylter>
      -DFIXTURES=${CMAKE_CURRENT_SOURCE_DIR}/tests/fixtures
      -DOUTPUT_DIR=${CMAKE_CURRENT_BINARY_DIR}/cli-tests
      -DPARALLEL=${PHYLTER2_HAVE_OPENMP}
      -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/cli-tests.cmake)
  endif()
endif()
