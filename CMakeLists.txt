cmake_minimum_required(VERSION 3.15)

project(Pacose VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)


if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Debug" CACHE STRING "Choose the type of build." FORCE)
endif()

set_property(GLOBAL PROPERTY USE_FOLDERS ON)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|AppleClang|GNU")
  add_compile_options(-Wall -Wuninitialized) # -Wextra -Wpedantic
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
  add_compile_options(/W4 /WX /EHsc)
endif()

file(GLOB_RECURSE SOURCE_FILES 
  ${PROJECT_SOURCE_DIR}/Main.cpp
  ${PROJECT_SOURCE_DIR}/maxSAT/*.cpp
)
list(REMOVE_ITEM SOURCE_FILES "${PROJECT_SOURCE_DIR}/solver-proxy/Main.cpp")

include_directories(
    #${PROJECT_SOURCE_DIR}
    ${PROJECT_SOURCE_DIR}/Helper
    ${PROJECT_SOURCE_DIR}/maxSAT
    ${PROJECT_SOURCE_DIR}/maxSAT/DGPW
)

find_library(MAXPRE_LIB NAMES maxpre HINTS ${CMAKE_SOURCE_DIR}/maxpre2/src/lib/)
if(NOT MAXPRE_LIB)
    message(FATAL_ERROR "maxpre library not found")
endif()

find_library(SOLVER_PROXY_LIB NAMES SATSolverProxy HINTS ${CMAKE_SOURCE_DIR}/solver-proxy/bin)
if(NOT MAXPRE_LIB)
    message(FATAL_ERROR "maxpre library not found")
endif()

add_executable(${PROJECT_NAME} ${SOURCE_FILES})
target_link_libraries(${PROJECT_NAME} PRIVATE 
  ${MAXPRE_LIB} 
  ${SOLVER_PROXY_LIB}
  z
)

set_target_properties(${PROJECT_NAME} PROPERTIES LINK_SEARCH_START_STATIC 1)
set_target_properties(${PROJECT_NAME} PROPERTIES LINK_SEARCH_END_STATIC 1)

if (CMAKE_SOURCE_DIR STREQUAL PROJECT_SOURCE_DIR)
  set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
  set(BUILD_SHARED_LIBS OFF)
  set(CMAKE_EXE_LINKER_FLAGS "-static")
endif()

set(CMAKE_CTEST_ARGUMENTS -V --output-on-failure)
enable_testing()
add_test(NAME RegressionSuite COMMAND ${PROJECT_SOURCE_DIR}/MaxSATRegressionSuite/run_tests.sh -s ${PROJECT_SOURCE_DIR}/bin/Pacose)
#add_test(NAME RegressionSuite COMMAND ${PROJECT_SOURCE_DIR}/MaxSATRegressionSuite/run_tests.sh ${PROJECT_SOURCE_DIR}/bin/Pacose) ## long test run
