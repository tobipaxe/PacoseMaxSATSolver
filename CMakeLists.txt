cmake_minimum_required(VERSION 3.15)

project(Pacose VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build." FORCE)
  #  set(CMAKE_BUILD_TYPE "Debug" CACHE STRING "Choose the type of build." FORCE)
endif()

# Use folders in the resulting project files.
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

# Set the output folder where your program will be created
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

# Add compiler warnings
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|AppleClang|GNU")
    if(CMAKE_BUILD_TYPE STREQUAL "Debug") # -Wextra -Wpedantic
      add_compile_options(-g -Wall -Wuninitialized)
    elseif(CMAKE_BUILD_TYPE STREQUAL "Release")
      add_compile_options(-O3)
    endif()
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
  add_compile_options(/W4 /WX /EHsc)
endif()

# Add your sources here
file(GLOB_RECURSE SOURCE_FILES 
  ${PROJECT_SOURCE_DIR}/Main.cpp
  ${PROJECT_SOURCE_DIR}/maxSAT/*.cpp
)
list(REMOVE_ITEM SOURCE_FILES "${PROJECT_SOURCE_DIR}/solver-proxy/Main.cpp")

execute_process(
    COMMAND make -C "${PROJECT_SOURCE_DIR}/solver-proxy/build"
)

# Include directories
include_directories(
    #${PROJECT_SOURCE_DIR}
    ${PROJECT_SOURCE_DIR}/Helper
    ${PROJECT_SOURCE_DIR}/maxSAT
    ${PROJECT_SOURCE_DIR}/maxSAT/DGPW
)

# Libraries
find_library(MAXPRE_LIB NAMES maxpre HINTS ${CMAKE_SOURCE_DIR}/maxpre2/src/lib/)
if(NOT MAXPRE_LIB)
    message(FATAL_ERROR "maxpre library not found")
endif()

# find_library(SOLVER_PROXY_LIB NAMES SATSolverProxy HINTS ${CMAKE_SOURCE_DIR}/solver-proxy/bin)
# if(NOT MAXPRE_LIB)
#     message(FATAL_ERROR "maxpre library not found")
# endif()

add_executable(${PROJECT_NAME} ${SOURCE_FILES})
target_link_libraries(${PROJECT_NAME} PRIVATE 
  # ${MAXPRE_LIB} 
  ${SOLVER_PROXY_LIB}
  z
)

# To ensure -static flag works as expected
set_target_properties(${PROJECT_NAME} PROPERTIES LINK_SEARCH_START_STATIC 1)
set_target_properties(${PROJECT_NAME} PROPERTIES LINK_SEARCH_END_STATIC 1)

# If we are building as a part of a larger project, we may not want to force these settings
if (CMAKE_SOURCE_DIR STREQUAL PROJECT_SOURCE_DIR)
  # Static build
  set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
  set(BUILD_SHARED_LIBS OFF)
  set(CMAKE_EXE_LINKER_FLAGS "-static")
endif()

set(CMAKE_CTEST_ARGUMENTS -V --output-on-failure)
enable_testing()
#add_test(NAME RegressionSuite COMMAND ${PROJECT_SOURCE_DIR}/MaxSATRegressionSuite/run_tests.sh -s ${PROJECT_SOURCE_DIR}/bin/Pacose)
add_test(NAME RegressionSuite COMMAND ${PROJECT_SOURCE_DIR}/MaxSATRegressionSuite/run_tests.sh ${PROJECT_SOURCE_DIR}/bin/Pacose) ## long test run
