cmake_minimum_required(VERSION 3.15)

project(SATSolverProxy VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

include_directories(
    ${PROJECT_SOURCE_DIR}
    ${PROJECT_SOURCE_DIR}/../VeriPB_ProofLogger)

set(LIB_SOURCE_FILES
    SATSolverProxy.cpp
    Glucose421SolverProxy.cpp
    #Glucose3SolverProxy.cpp
    # CryptoMiniSat5SolverProxy.cpp
    #AntomSolverProxy.cpp
    CadicalSolverProxy.cpp
    "${PROJECT_SOURCE_DIR}/../VeriPB_Prooflogger/VeriPBProoflogger_Pacose.cc"
    "${PROJECT_SOURCE_DIR}/../VeriPB_Prooflogger/MaxSATProoflogger_Pacose.cc"
    "${PROJECT_SOURCE_DIR}/../VeriPB_Prooflogger/cadicalprooftracer.cpp"
)

add_library(solver-proxy-static STATIC ${LIB_SOURCE_FILES})
add_library(solver-proxy-shared SHARED ${LIB_SOURCE_FILES})

set_target_properties(solver-proxy-static PROPERTIES OUTPUT_NAME "SATSolverProxy")
set_target_properties(solver-proxy-shared PROPERTIES OUTPUT_NAME "SATSolverProxy")

add_custom_target(glucose-4.2.1-build
    COMMAND mkdir -p build && cmake -S . -B build && make -C build -j 
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/glucose-4.2.1)
add_library(glucose-4.2.1-static STATIC IMPORTED)
set_target_properties(glucose-4.2.1-static PROPERTIES IMPORTED_LOCATION ${PROJECT_SOURCE_DIR}/glucose-4.2.1/build/libglucose.a)
add_dependencies(glucose-4.2.1-static glucose-4.2.1-build)

#target_include_directories(glucose-4.2.1-static INTERFACE
#  ${PROJECT_SOURCE_DIR}/glucose-4.2.1/core
#  ${PROJECT_SOURCE_DIR}/glucose-4.2.1/mtl
#  ${PROJECT_SOURCE_DIR}/glucose-4.2.1/parallel
#  ${PROJECT_SOURCE_DIR}/glucose-4.2.1/simp
#  ${PROJECT_SOURCE_DIR}/glucose-4.2.1/utils
#)
#add_library(glucose-3.0 STATIC IMPORTED)
#set_target_properties(glucose-3.0 PROPERTIES IMPORTED_LOCATION ${PROJECT_SOURCE_DIR}/glucose-3.0/core/libglucose-3.0_release.a)
add_custom_target(cadical-build
   # COMMAND mkdir -p build/cadical && cd build/cadical && ../../configure && make -j
    COMMAND make -j
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/cadical)
add_library(cadical-static STATIC IMPORTED)
set_target_properties(cadical-static PROPERTIES IMPORTED_LOCATION ${PROJECT_SOURCE_DIR}/cadical/build/libcadical.a)
add_dependencies(cadical-static cadical-build)
#target_include_directories(cadical-static INTERFACE
#  ${PROJECT_SOURCE_DIR}/cadical/src
#)
#add_library(antom STATIC IMPORTED)
#set_target_properties(antom PROPERTIES IMPORTED_LOCATION ${PROJECT_SOURCE_DIR}/antom/implementation/trunk/build/liblibantom.a)
#add_library(cms STATIC IMPORTED)
#set_target_properties(cms PROPERTIES IMPORTED_LOCATION ${PROJECT_SOURCE_DIR}/CryptoMiniSAT/build/lib/libcryptominisat5.so)

target_link_libraries(solver-proxy-static PRIVATE glucose-4.2.1-static cadical-static)
IF(APPLE)
    target_link_libraries(solver-proxy-shared PUBLIC glucose-4.2.1-static cadical-static)
ENDIF(APPLE)

add_executable(SATSolverProxy Main.cpp)
target_link_libraries(SATSolverProxy solver-proxy-static)

