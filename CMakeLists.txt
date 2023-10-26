cmake_minimum_required(VERSION 3.15)

project(Pacose VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

add_subdirectory(${PROJECT_SOURCE_DIR}/solver-proxy)
add_custom_command(
  OUTPUT ${PROJECT_SOURCE_DIR}/maxpre2/src/lib/libmaxpre.a
  COMMAND "make lib"
  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/maxpre2)

add_library(maxpre-static STATIC IMPORTED)
set_target_properties(maxpre-static PROPERTIES IMPORTED_LOCATION ${PROJECT_SOURCE_DIR}/maxpre2/src/lib/libmaxpre.a)
#target_include_directories(maxpre-static INTERFACE
#  ${PROJECT_SOURCE_DIR}/maxpre/src
#)

# Include directories
include_directories(
    #${PROJECT_SOURCE_DIR}
    ${PROJECT_SOURCE_DIR}/Helper
    ${PROJECT_SOURCE_DIR}/maxSAT
    ${PROJECT_SOURCE_DIR}/maxSAT/DGPW
    ${PROJECT_SOURCE_DIR}/VeriPB_ProofLogger
)

file(GLOB SOURCE_FILES 
  "${PROJECT_SOURCE_DIR}/Main.cpp"
  "${PROJECT_SOURCE_DIR}/maxSAT/*.cpp"
  "${PROJECT_SOURCE_DIR}/maxSAT/DGPW/*.cpp"
  "${PROJECT_SOURCE_DIR}/VeriPB_Prooflogger/VeriPBProoflogger_Pacose.cc"
  "${PROJECT_SOURCE_DIR}/VeriPB_Prooflogger/MaxSATProoflogger_Pacose.cc"
  "${PROJECT_SOURCE_DIR}/VeriPB_Prooflogger/cadicalprooftracer.cpp"
)
add_executable(Pacose ${SOURCE_FILES})
target_link_libraries(Pacose PRIVATE 
  solver-proxy-static
  maxpre-static
  z
)

set(CMAKE_CTEST_ARGUMENTS -V --output-on-failure)
enable_testing()
add_test(NAME RegressionSuite COMMAND ${PROJECT_SOURCE_DIR}/MaxSATRegressionSuite/run_tests.sh -s ${PROJECT_SOURCE_DIR}/bin/Pacose)
#add_test(NAME RegressionSuite COMMAND ${PROJECT_SOURCE_DIR}/MaxSATRegressionSuite/run_tests.sh ${PROJECT_SOURCE_DIR}/bin/Pacose) ## long test run
