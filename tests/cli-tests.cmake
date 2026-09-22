if(NOT DEFINED PHYLTER OR NOT DEFINED FIXTURES OR NOT DEFINED OUTPUT_DIR)
  message(FATAL_ERROR "Missing CLI test paths")
endif()
file(MAKE_DIRECTORY "${OUTPUT_DIR}")
execute_process(COMMAND "${PHYLTER}" --check --trees "${FIXTURES}/carnivora.nwk"
  RESULT_VARIABLE status OUTPUT_VARIABLE stdout ERROR_VARIABLE stderr)
if(NOT status EQUAL 0 OR NOT stdout MATCHES "Genes:[ ]+125" OR
   NOT stdout MATCHES "Phylter 0.1.0-dev" OR
   NOT stdout MATCHES "Tips per gene:" OR
   NOT stdout MATCHES "Patristic distance:[ ]+available" OR
   NOT stdout MATCHES "carnivora:1, carnivora:2" OR
   NOT stdout MATCHES "tree dataset is well formed")
  message(FATAL_ERROR "Input check failed: ${stdout}${stderr}")
endif()
execute_process(COMMAND "${PHYLTER}" --check --trees "${FIXTURES}/carnivora.nwk" --distance nodal
  RESULT_VARIABLE status OUTPUT_QUIET ERROR_VARIABLE stderr)
if(status EQUAL 0 OR NOT stderr MATCHES "--check only accepts --trees")
  message(FATAL_ERROR "Input check accepted an analysis option: ${stderr}")
endif()
set(branchless_tree "${OUTPUT_DIR}/branchless-input.nwk")
file(WRITE "${branchless_tree}" "((A,B),C);\n")
execute_process(COMMAND "${PHYLTER}" --check --trees "${branchless_tree}"
  RESULT_VARIABLE status OUTPUT_VARIABLE stdout ERROR_VARIABLE stderr)
if(NOT status EQUAL 0 OR NOT stdout MATCHES "Patristic distance: unavailable" OR
   NOT stdout MATCHES "Nodal distance:[ ]+available" OR
   NOT stdout MATCHES "branchless-input")
  message(FATAL_ERROR "Branchless-tree capabilities are incorrect: ${stdout}${stderr}")
endif()
execute_process(COMMAND "${PHYLTER}" run --trees "${FIXTURES}/carnivora.nwk"
  --out "${OUTPUT_DIR}/obsolete-command"
  RESULT_VARIABLE status OUTPUT_QUIET ERROR_VARIABLE stderr)
if(status EQUAL 0 OR NOT stderr MATCHES "unknown option: run")
  message(FATAL_ERROR "Obsolete run command was not rejected clearly: ${stderr}")
endif()
set(invalid_tree "${OUTPUT_DIR}/invalid-input.nwk")
file(WRITE "${invalid_tree}" "(A:1,B:broken,C:1);\n")
execute_process(COMMAND "${PHYLTER}" --check --trees "${invalid_tree}"
  RESULT_VARIABLE status OUTPUT_QUIET ERROR_VARIABLE stderr)
if(status EQUAL 0 OR NOT stderr MATCHES "invalid-input\\.nwk: invalid Newick")
  message(FATAL_ERROR "Invalid input error does not name its file: ${stderr}")
endif()
execute_process(COMMAND "${PHYLTER}" --trees "${FIXTURES}/carnivora.nwk"
  --out "${OUTPUT_DIR}/serial" --force
  RESULT_VARIABLE status OUTPUT_VARIABLE stdout ERROR_VARIABLE stderr)
if(NOT status EQUAL 0)
  message(FATAL_ERROR "CLI analysis failed: ${stderr}")
endif()
foreach(section "Dataset" "run phylter --help for all options"
                "Distance scaling:" "Acceptance rule:" "Estimated memory:" "Optimization"
                "candidate cell outlier"
                "Results" "Elapsed time:" "Peak memory:" "Output files"
                "Note: rerun with --diagnostics DIR")
  if(NOT stdout MATCHES "${section}")
    message(FATAL_ERROR "CLI output is missing section '${section}': ${stdout}")
  endif()
endforeach()
if(NOT stdout MATCHES "Data loss:[ ]+1.42% \\(94/6625 outlier cells removed\\)")
  message(FATAL_ERROR "CLI result summary is incomplete: ${stdout}")
endif()
file(READ "${OUTPUT_DIR}/serial.outliers.tsv" actual)
file(READ "${FIXTURES}/carnivora.outliers.tsv" expected)
string(REPLACE "\r\n" "\n" actual "${actual}")
string(REPLACE "\r\n" "\n" expected "${expected}")
if(NOT actual STREQUAL "gene\tspecies\n${expected}")
  message(FATAL_ERROR "CLI output does not match frozen R outliers")
endif()
execute_process(COMMAND "${PHYLTER}" --trees "${FIXTURES}/carnivora.nwk"
  --out "${OUTPUT_DIR}/dense" --rv-method dense --force
  RESULT_VARIABLE status OUTPUT_QUIET ERROR_VARIABLE stderr)
if(NOT status EQUAL 0)
  message(FATAL_ERROR "Dense reference method failed: ${stderr}")
endif()
execute_process(COMMAND "${CMAKE_COMMAND}" -E compare_files
  "${OUTPUT_DIR}/serial.outliers.tsv" "${OUTPUT_DIR}/dense.outliers.tsv"
  RESULT_VARIABLE status)
if(NOT status EQUAL 0)
  message(FATAL_ERROR "Dense and matrix-free outliers differ")
endif()
execute_process(COMMAND "${PHYLTER}" --trees "${FIXTURES}/carnivora.nwk"
  --out "${OUTPUT_DIR}/invalid-method" --rv-method unknown
  RESULT_VARIABLE status OUTPUT_QUIET ERROR_QUIET)
if(status EQUAL 0)
  message(FATAL_ERROR "Unknown RV method was accepted")
endif()
execute_process(COMMAND "${PHYLTER}" --trees "${FIXTURES}/carnivora.nwk"
  --out "${OUTPUT_DIR}/serial" RESULT_VARIABLE status OUTPUT_QUIET ERROR_QUIET)
if(status EQUAL 0)
  message(FATAL_ERROR "CLI overwrote outputs without --force")
endif()
execute_process(COMMAND "${PHYLTER}" --trees "${FIXTURES}/carnivora.nwk"
  --out "${OUTPUT_DIR}/invalid" --k nan RESULT_VARIABLE status OUTPUT_QUIET ERROR_QUIET)
if(status EQUAL 0)
  message(FATAL_ERROR "CLI accepted NaN threshold")
endif()
if(PARALLEL)
  execute_process(COMMAND "${PHYLTER}" --trees "${FIXTURES}/carnivora.nwk"
    --out "${OUTPUT_DIR}/parallel" --threads 4 --force
    RESULT_VARIABLE status OUTPUT_QUIET ERROR_VARIABLE stderr)
  if(NOT status EQUAL 0)
    message(FATAL_ERROR "Parallel analysis failed: ${stderr}")
  endif()
  foreach(suffix outliers.tsv scores.tsv discarded.tsv)
    execute_process(COMMAND "${CMAKE_COMMAND}" -E compare_files
      "${OUTPUT_DIR}/serial.${suffix}" "${OUTPUT_DIR}/parallel.${suffix}"
      RESULT_VARIABLE status)
    if(NOT status EQUAL 0)
      message(FATAL_ERROR "Parallel output differs: ${suffix}")
    endif()
  endforeach()
endif()
