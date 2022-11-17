# Copyright (c) 2017, Mate Soos
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

include_directories( ${PROJECT_SOURCE_DIR} )
include_directories( ${BOSPHORUS_INCLUDE_DIRS} )

if (NOT WIN32)
    add_cxx_flag_if_supported("-Wno-bitfield-constant-conversion")
    #add_cxx_flag_if_supported("-Wduplicated-cond")
    #add_cxx_flag_if_supported("-Wduplicated-branches")
    add_cxx_flag_if_supported("-Wlogical-op")
    add_cxx_flag_if_supported("-Wrestrict")
    add_cxx_flag_if_supported("-Wnull-dereference")
    add_cxx_flag_if_supported("-Wjump-misses-init")
    add_cxx_flag_if_supported("-Wdouble-promotion")
    add_cxx_flag_if_supported("-Wshadow")
    add_cxx_flag_if_supported("-Wformat=2")
    add_cxx_flag_if_supported("-Wextra-semi")
    add_cxx_flag_if_supported("-pedantic")
    #add_cxx_flag_if_supported("-Wdeprecated")
endif()

include_directories(${CMAKE_CURRENT_BINARY_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

include_directories(${Boost_INCLUDE_DIRS})
link_directories(Boost_LIBRARY_DIRS)

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/GitSHA1.cpp.in" "${CMAKE_CURRENT_BINARY_DIR}/GitSHA1.cpp" @ONLY)
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/cryptominisat.h.in" "${CMAKE_CURRENT_BINARY_DIR}/cmsgen/cryptominisat.h" @ONLY)
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/solvertypesmini.h.in" "${CMAKE_CURRENT_BINARY_DIR}/cmsgen/solvertypesmini.h" @ONLY)
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/cryptominisat_c.h.in" "${CMAKE_CURRENT_BINARY_DIR}/cmsgen/cryptominisat_c.h" @ONLY)

set(cryptoms_lib_files
    cnf.cpp
    propengine.cpp
    varreplacer.cpp
    clausecleaner.cpp
    clauseusagestats.cpp
    prober.cpp
    occsimplifier.cpp
    subsumestrengthen.cpp
    clauseallocator.cpp
    sccfinder.cpp
    solverconf.cpp
    distillerlong.cpp
    distillerlongwithimpl.cpp
    str_impl_w_impl_stamp.cpp
    solutionextender.cpp
    completedetachreattacher.cpp
    searcher.cpp
    solver.cpp
    implcache.cpp
    stamp.cpp
    hyperengine.cpp
    subsumeimplicit.cpp
    datasync.cpp
    reducedb.cpp
    intree.cpp
    searchstats.cpp
    xorfinder.cpp
    cryptominisat_c.cpp
    walksat.cpp
    sls.cpp
    gaussian.cpp
    matrixfinder.cpp
    packedrow.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/GitSHA1.cpp
)

set(cryptoms_lib_link_libs "")

add_library(cmsgen
    ${cryptoms_lib_files}
    cryptominisat.cpp
)

GENERATE_EXPORT_HEADER(cmsgen
         BASE_NAME cmsgen
         #EXPORT_MACRO_NAME cmsgen_EXPORT
         #EXPORT_FILE_NAME MyLibrary_Export.h
         #STATIC_DEFINE MyLibrary_BUILT_AS_STATIC
)

# indicate that we depend on pthread, and compile in the actual library
target_link_libraries(cmsgen
    LINK_PUBLIC ${cryptoms_lib_link_libs}
    LINK_PUBLIC ${CMAKE_THREAD_LIBS_INIT}
)

if (NOT WIN32)
    set_target_properties(cmsgen PROPERTIES
        PUBLIC_HEADER "${cmsgen_public_headers}"
        VERSION ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}
        SOVERSION ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}
    )
else()
    set_target_properties(cmsgen PROPERTIES
        OUTPUT_NAME cmsgenwin
        PUBLIC_HEADER "${cmsgen_public_headers}"
        VERSION ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}
        SOVERSION ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}
    )
endif()

# -----------------------------------------------------------------------------
# Copy public headers into build directory include directory.
# The cmsgenConfig.cmake we generate in the build directory depends on
# this.
# -----------------------------------------------------------------------------
set(HEADER_DEST "${PROJECT_BINARY_DIR}/include/cmsgen")
add_custom_target(CopyPublicHeaders ALL)
get_target_property(cmsgen_public_headers cmsgen PUBLIC_HEADER)
foreach(public_header ${cmsgen_public_headers})
    get_filename_component(HEADER_NAME ${public_header} NAME)
    add_custom_command(TARGET CopyPublicHeaders PRE_BUILD
                       COMMAND ${CMAKE_COMMAND} -E make_directory
                               "${HEADER_DEST}"
                       COMMAND ${CMAKE_COMMAND} -E echo
                       "Copying ${HEADER_NAME} to ${HEADER_DEST}"
                       COMMAND ${CMAKE_COMMAND} -E
                           copy_if_different
                           ${public_header}
                           "${HEADER_DEST}"
                      )
endforeach()

install(TARGETS cmsgen
    EXPORT ${CRYPTOMINISAT5_EXPORT_NAME}
    LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
    ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
    RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
    PUBLIC_HEADER DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/cmsgen"
)

if (NOT ONLY_SIMPLE)
    add_executable(cmsgen-bin
        main.cpp
        main_common.cpp
        main_exe.cpp
        signalcode.cpp
    )
endif()

set(cryptoms_exec_link_libs
    cmsgen
)

IF (ZLIB_FOUND)
    SET(cryptoms_exec_link_libs ${cryptoms_exec_link_libs} ${ZLIB_LIBRARY})
ENDIF()

if (NOT ONLY_SIMPLE)
    set_target_properties(cmsgen-bin PROPERTIES
        OUTPUT_NAME cmsgen
        RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}
        INSTALL_RPATH_USE_LINK_PATH TRUE)
    target_link_libraries(cmsgen-bin
        ${cryptoms_exec_link_libs}
        "${Boost_LIBRARIES}"
    )
    install(TARGETS cmsgen-bin
        # EXPORT ${CRYPTOMINISAT5_EXPORT_NAME}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    )
    SET(CPACK_PACKAGE_EXECUTABLES "cmsgen")
endif()
