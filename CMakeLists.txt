cmake_minimum_required(VERSION 3.14)
project(cppsolver VERSION "0.3.0")

configure_file(${PROJECT_SOURCE_DIR}/config.hpp.in ${PROJECT_BINARY_DIR}/config.hpp)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
        MESSAGE("CODE COVERAGE")
        INCLUDE(CodeCoverage)
        APPEND_COVERAGE_COMPILER_FLAGS()
    endif()
    SET(CMAKE_VERBOSE_MAKEFILE on)
    add_definitions(-pedantic -Wall)
endif()

set(SOLVER_LIB solverlib)

include_directories(libraries/eigen)

set(Boost_USE_MULTITHREADED ON)
FIND_PACKAGE( Boost 1.63 COMPONENTS program_options
        system filesystem REQUIRED )

INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIR})
SET(requiredlibs ${requiredlibs} ${Boost_LIBRARIES})
link_directories(${Boost_LIBRARY_DIRS})

# include so we find the config.hpp file
include_directories(${PROJECT_BINARY_DIR})

include_directories(src)
add_subdirectory(src)

if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
        add_subdirectory(test)
    endif()
endif()

set(SOURCE_FILES main.cpp)
add_executable(cppsolver ${SOURCE_FILES})
target_link_libraries(cppsolver ${SOLVER_LIB})