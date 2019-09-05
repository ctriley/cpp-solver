project(cppsolver)

SET(SOURCE_FILES_SRC
        simplex.cpp)

add_library(${SOLVER_LIB} SHARED ${SOURCE_FILES_SRC})
target_link_libraries(${SOLVER_LIB} ${requiredlibs})