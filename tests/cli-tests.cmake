if(NOT DEFINED PHYLTER OR NOT DEFINED FIXTURES OR NOT DEFINED OUTPUT_DIR)
  message(FATAL_ERROR "Missing CLI test paths")
endif()
file(MAKE_DIRECTORY "${OUTPUT_DIR}")
execute_process(COMMAND "${PHYLTER}" run --trees "${FIXTURES}/carnivora.nwk"
  --out "${OUTPUT_DIR}/serial" --force
  RESULT_VARIABLE status OUTPUT_VARIABLE stdout ERROR_VARIABLE stderr)
if(NOT status EQUAL 0)
  message(FATAL_ERROR "CLI analysis failed: ${stderr}")
endif()
file(READ "${OUTPUT_DIR}/serial.outliers.tsv" actual)
file(READ "${FIXTURES}/carnivora.outliers.tsv" expected)
string(REPLACE "\r\n" "\n" actual "${actual}")
string(REPLACE "\r\n" "\n" expected "${expected}")
if(NOT actual STREQUAL "gene\tspecies\n${expected}")
  message(FATAL_ERROR "CLI output does not match frozen R outliers")
endif()
execute_process(COMMAND "${PHYLTER}" run --trees "${FIXTURES}/carnivora.nwk"
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
execute_process(COMMAND "${PHYLTER}" run --trees "${FIXTURES}/carnivora.nwk"
  --out "${OUTPUT_DIR}/invalid-method" --rv-method unknown
  RESULT_VARIABLE status OUTPUT_QUIET ERROR_QUIET)
if(status EQUAL 0)
  message(FATAL_ERROR "Unknown RV method was accepted")
endif()
execute_process(COMMAND "${PHYLTER}" run --trees "${FIXTURES}/carnivora.nwk"
  --out "${OUTPUT_DIR}/serial" RESULT_VARIABLE status OUTPUT_QUIET ERROR_QUIET)
if(status EQUAL 0)
  message(FATAL_ERROR "CLI overwrote outputs without --force")
endif()
execute_process(COMMAND "${PHYLTER}" run --trees "${FIXTURES}/carnivora.nwk"
  --out "${OUTPUT_DIR}/invalid" --k nan RESULT_VARIABLE status OUTPUT_QUIET ERROR_QUIET)
if(status EQUAL 0)
  message(FATAL_ERROR "CLI accepted NaN threshold")
endif()
if(PARALLEL)
  execute_process(COMMAND "${PHYLTER}" run --trees "${FIXTURES}/carnivora.nwk"
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
