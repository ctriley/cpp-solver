# Source: https://github.com/libigl/CoMISo/blob/master/cmake/FindGUROBI.cmake

# - Try to find GUROBI

#  GUROBI_BASE - The libraries needed to use Gurobi

# Once done this will define
#  GUROBI_FOUND - System has Gurobi
#  GUROBI_INCLUDE_DIRS - The Gurobi include directories
#  GUROBI_LIBRARIES - The libraries needed to use Gurobi

set (GUROBI_ENABLE ON CACHE BOOL "Enable gurobi?")

if ( GUROBI_ENABLE )
    set (GUROBI_BASE "c:" CACHE PATH "Base path of your gurobi installation")

    if (GUROBI_INCLUDE_DIR)
        # in cache already
        set(GUROBI_FOUND TRUE)
        set(GUROBI_INCLUDE_DIRS "${GUROBI_INCLUDE_DIR}" )
        set(GUROBI_LIBRARIES "${GUROBI_CXX_LIBRARY};${GUROBI_LIBRARY}" )
    else (GUROBI_INCLUDE_DIR)



        find_path(GUROBI_INCLUDE_DIR
                NAMES gurobi_c++.h
                PATHS "$ENV{GUROBI_HOME}/include"
                "/usr/local/pacerepov1/gurobi/8.1.0/linux64/include"
                "/usr/local/pacerepov1/gurobi/7.5.1/linux64/include"
		        "/sw/arc/centos7/gurobi/gurobi751/linux64/include"
                "/sw/arc/centos7/gurobi/gurobi652/linux64/include"
                "/opt/gurobi811/linux64/include"
               	"/opt/gurobi/8.1.0/include/"
                "/opt/gurobi800/linux64/include"
                "/opt/gurobi752/linux64/include"
                "/Library/gurobi801/mac64/include"
                "/Library/gurobi800/mac64/include"
                "/Library/gurobi752/mac64/include"
                "/Library/gurobi702/mac64/include"
                "/Library/gurobi652/mac64/include"
                "/Library/gurobi562/mac64/include"
                "/Library/gurobi502/mac64/include"
                "C:\\libs\\gurobi502\\include"
                "C:\\libs\\gurobi562\\include"
                "C:\\libs\\gurobi652\\include"
                "C:\\libs\\gurobi702\\include"
                "C:\\libs\\gurobi752\\include"
                "${GUROBI_BASE}/include"
                )
	MESSAGE(${GUROBI_INCLUDE_DIR})
        find_library( GUROBI_LIBRARY
                NAMES gurobi
                gurobi81
                gurobi80
                gurobi75
                gurobi70
                gurobi65
                gurobi60
                gurobi56
                gurobi55
                gurobi51
                gurobi50
                gurobi46
                gurobi45

                PATHS "$ENV{GUROBI_HOME}/lib"
                "/nv/hpvh1/criley38/data/software/lib"
                "/home/ctriley/didi-simulator/lib"
                "/sw/arc/centos7/gurobi/gurobi751/linux64/lib"
                "/sw/arc/centos7/gurobi/gurobi652/linux64/lib"
                "/opt/gurobi811/linux64/lib"
                "/opt/gurobi/8.1.0/linux64/lib"
                "/opt/gurobi752/linux64/lib"
                "/Library/gurobi801/mac64/lib"
                "/Library/gurobi800/mac64/lib"
                "/Library/gurobi752/mac64/lib"
                "/Library/gurobi702/mac64/lib"
                "/Library/gurobi652/mac64/lib"
                "/Library/gurobi562/mac64/lib"
                "/Library/gurobi502/mac64/lib"
                "C:\\libs\\gurobi752\\lib"
                "C:\\libs\\gurobi702\\lib"
                "C:\\libs\\gurobi652\\lib"
                "C:\\libs\\gurobi562\\lib"
                "C:\\libs\\gurobi502\\lib"
                "${GUROBI_BASE}/lib"
                )

        if ( CMAKE_GENERATOR MATCHES "^Visual Studio 12.*Win64" )
            SET(GUROBI_LIB_NAME "gurobi_c++md2013")
        endif()

        find_library( GUROBI_CXX_LIBRARY
                NAMES gurobi_c++
                ${GUROBI_LIB_NAME}
                PATHS "$ENV{GUROBI_HOME}/lib"
                "/nv/hpvh1/criley38/data/software/lib"
		        "/home/ctriley/didi-simulator/lib"
                "/sw/arc/centos7/gurobi/gurobi751/linux64/lib"
                "/sw/arc/centos7/gurobi/gurobi652/linux64/lib"
                "/opt/gurobi811/linux64/lib"
                "/opt/gurobi/8.1.0/linux64/lib"
                "/opt/gurobi752/linux64/lib"
                "/Library/gurobi801/mac64/lib"
                "/Library/gurobi800/mac64/lib"
                "/Library/gurobi752/mac64/lib"
                "/Library/gurobi702/mac64/lib"
                "/Library/gurobi652/mac64/lib"
                "/Library/gurobi562/mac64/lib"
                "/Library/gurobi502/mac64/lib"
                "C:\\libs\\gurobi752\\lib"
                "C:\\libs\\gurobi702\\lib"
                "C:\\libs\\gurobi652\\lib"
                "C:\\libs\\gurobi562\\lib"
                "C:\\libs\\gurobi502\\lib"
                "${GUROBI_BASE}/lib"
                )

        # Binary dir for DLLs
        find_path(GUROBI_BIN_DIR
                NAMES "gurobi60.dll"
                PATHS "${GUROBI_INCLUDE_DIR}/../bin"
                "${GUROBI_BASE}/bin"
                DOC "Directory containing the GUROBI DLLs"
                )

        set(GUROBI_INCLUDE_DIRS "${GUROBI_INCLUDE_DIR}" )
        set(GUROBI_LIBRARIES "${GUROBI_CXX_LIBRARY};${GUROBI_LIBRARY}" )

        # use c++ headers as default
        #set(GUROBI_COMPILER_FLAGS "-DIL_STD" CACHE STRING "Gurobi Compiler Flags")

        include(FindPackageHandleStandardArgs)
        # handle the QUIETLY and REQUIRED arguments and set LIBCPLEX_FOUND to TRUE
        # if all listed variables are TRUE
        find_package_handle_standard_args(GUROBI  DEFAULT_MSG
                GUROBI_CXX_LIBRARY GUROBI_LIBRARY GUROBI_INCLUDE_DIR)

        mark_as_advanced(GUROBI_INCLUDE_DIR GUROBI_LIBRARY GUROBI_CXX_LIBRARY GUROBI_BIN_DIR )

    endif(GUROBI_INCLUDE_DIR)

endif()
