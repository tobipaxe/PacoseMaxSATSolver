cmake_minimum_required(VERSION 3.00)                 # CMake version check
project("SATSolverProxy")                                      # Create project "simple_example"

set(CMAKE_CXX_STANDARD 11)

set (EXTRA_OPTIONS "-O3")
#set (EXTRA_OPTIONS "-g -Wall -Wextra")

#set(CMAKE_BUILD_TYPE "Debug")                          #Possible values are empty, Debug, Release, RelWithDebInfo and MinSizeRel.
set(CMAKE_BUILD_TYPE Release)

# Set the output folder where your program will be created
set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/bin)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR})

include_directories("$(PROJECT_SOURCE_DIR)")

set(SOURCE_FILES
    Main.cpp
    SATSolverProxy.cpp
    Glucose421SolverProxy.cpp
    #Glucose3SolverProxy.cpp
    # CryptoMiniSat5SolverProxy.cpp
    #AntomSolverProxy.cpp
    CadicalSolverProxy.cpp
    ) 

add_library(lib${PROJECT_NAME} STATIC ${SOURCE_FILES})

add_library(glucose-4.2.1 STATIC IMPORTED)
set_target_properties(glucose-4.2.1 PROPERTIES IMPORTED_LOCATION ${PROJECT_SOURCE_DIR}/glucose-4.2.1/build/libglucose.a)
#add_library(glucose-3.0 STATIC IMPORTED)
#set_target_properties(glucose-3.0 PROPERTIES IMPORTED_LOCATION ${PROJECT_SOURCE_DIR}/glucose-3.0/core/libglucose-3.0_release.a)
add_library(cadical STATIC IMPORTED)
set_target_properties(cadical PROPERTIES IMPORTED_LOCATION ${PROJECT_SOURCE_DIR}/cadical/build/libcadical.a)
#add_library(antom STATIC IMPORTED)
#set_target_properties(antom PROPERTIES IMPORTED_LOCATION ${PROJECT_SOURCE_DIR}/antom/implementation/trunk/build/liblibantom.a)
#add_library(cms STATIC IMPORTED)
#set_target_properties(cms PROPERTIES IMPORTED_LOCATION ${PROJECT_SOURCE_DIR}/CryptoMiniSAT/build/lib/libcryptominisat5.so)
add_custom_target(combined
  COMMAND ar -x $<TARGET_FILE:glucose-4.2.1>
#  COMMAND ar -x $<TARGET_FILE:glucose-3.0>
  COMMAND ar -x $<TARGET_FILE:cadical>
  #COMMAND ar -x $<TARGET_FILE:antom>
  #COMMAND ar -x $<TARGET_FILE:cms>
  COMMAND ar -x $<TARGET_FILE:lib${PROJECT_NAME}>
  #COMMAND echo "ar -x $<TARGET_FILE:lib${PROJECT_NAME}>"
  # COMMAND ar -qcs lib${PROJECT_NAME}.a *.o *.or
  COMMAND ar -qcs lib${PROJECT_NAME}.a *.o
  # COMMAND rm *.o *.or
  COMMAND rm *.o
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
  DEPENDS glucose-4.2.1 lib${PROJECT_NAME}
 # DEPENDS glucose-3.0 lib${PROJECT_NAME}
  DEPENDS cadical lib${PROJECT_NAME}
  #DEPENDS antom lib${PROJECT_NAME}
 # DEPENDS cms lib${PROJECT_NAME}
)

add_library(SSolverProxy STATIC IMPORTED GLOBAL)
add_dependencies(SSolverProxy combined)
set_target_properties(SSolverProxy PROPERTIES IMPORTED_LOCATION ${CMAKE_BINARY_DIR}/lib${PROJECT_NAME}.a)

add_executable(${PROJECT_NAME} ${SOURCE_FILES})
# target_link_libraries(${PROJECT_NAME} SSolverProxy /usr/lib/x86_64-linux-gnu/libsqlite3.so.0)
target_link_libraries(${PROJECT_NAME} SSolverProxy)
# target_link_libraries(${PROJECT_NAME} SSolverProxy /usr/lib/x86_64-linux-gnu/libsqlite3.so.0 /usr/lib/x86_64-linux-gnu/libm4ri.a)

