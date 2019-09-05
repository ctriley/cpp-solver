project(cppsolver-test)

set(SOURCE_FILES simplextest.cpp)

add_executable(src_test ${SOURCE_FILES})

target_link_libraries(src_test ${SOLVER_LIB})
target_link_libraries(src_test gtest gtest_main)

#add_test(NAME all COMMAND state-test)
install(TARGETS src_test DESTINATION bin)

if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
    SETUP_TARGET_FOR_COVERAGE(
        src_coverage  # Name for custom target.
        src-test         # Name of the test driver executable that runs the tests.
        # NOTE! This should always have a ZERO as exit code
        # otherwise the coverage generation will not complete.
        src-coverage            # Name of output directory.
    )
endif()
